import gzip
import os
import argparse
import numpy as np
from collections import defaultdict
import itertools
from multiprocessing import Pool


def batched(iterable, n, *, strict=False):
    if n < 1:
        raise ValueError('n must be at least one')
    iterator = iter(iterable)
    while batch := tuple(itertools.islice(iterator, n)):
        if strict and len(batch) != n:
            raise ValueError('batched(): incomplete batch')
        yield batch


def pairwise(iterable):
    iterator = iter(iterable)
    a = next(iterator, None)
    for b in iterator:
        yield a, b
        a = b


def cmd_int(cmd):
    with os.popen(cmd) as h1:
        return int(list(h1)[0])


def uint16(data):
    return np.frombuffer(data, dtype=np.uint16)[0]


def get_header_size(bam):
    with gzip.open(bam, 'rt') as h1:
        return sum(len(i) for i in itertools.takewhile(lambda x: x.startswith('@'), h1))


def index_bgzip(file):
    total = bytearray()
    with open(file, 'rb') as h1:
        while tmp1 := h1.read(18)[-2:]:
            block_size = int(uint16(tmp1)) - 17
            tmp2 = h1.read(block_size)[-4:]
            total.extend(tmp1 + tmp2)
    total = np.frombuffer(total, dtype=[('zip', np.uint16), ('unzip', np.uint32)]).copy()
    total['zip'] += 1
    if total['zip'].sum() != os.path.getsize(file):
        raise ValueError(f'{file} file format error.')
    bam_index = np.full(len(total) * 2 + 1, len(total), dtype=np.uint64)
    bam_index[1::2] = np.cumsum(total['zip'], dtype=np.uint64)
    bam_index[2::2] = np.cumsum(total['unzip'], dtype=np.uint64)
    index_name = f'{file}.gzi'
    bam_index.tofile(index_name)


def chunk_sam(sam_file, *, chunk_num):
    start_pos = get_header_size(sam_file)
    total = int(np.fromfile(f'{sam_file}.gzi', dtype=np.uint64)[-1])
    chunk_size = int((total - start_pos) / chunk_num) + 1
    keep_pos = [start_pos]
    while start_pos < total:
        start_pos += chunk_size
        if start_pos >= total:
            keep_pos.append(total)
            break
        start_pos += cmd_int(f'bgzip -I {sam_file}.gzi -b {start_pos} {sam_file} | head -1 | wc -c')
        keep_pos.append(start_pos)
    return [(num, start, end - start) for num, (start, end) in enumerate(pairwise(keep_pos))]


def get_sub(*, bam, cb_group, nt, out_name, tag_name='CB', snp_site, ref, chunk, hap):
    tmp_dir = f'{out_name}-dir'
    if os.path.exists(tmp_dir):
        raise FileExistsError(f'{tmp_dir} exits.')
    os.system(f'mkdir {tmp_dir}')
    
    with open(cb_group) as h1:
        data = [i.split() for i in h1]
        if any(len(i) != 2 for i in data):
            raise RuntimeError(f'{cb_group} requires with two columns.')
        group_dict = defaultdict(list)
        for row in data:
            group_dict[row[-1]].append(row[0])
    cb_files = []
    for g_name, val in group_dict.items():
        cb_file = f'{tmp_dir}/{g_name}-cb.txt'
        cb_files.append([g_name, cb_file])
        with open(cb_file, 'wt') as h1:
            for item in val:
                print(item, file=h1)
    header = f'{tmp_dir}/pub-header.txt'
    os.system(rf'''samtools view -H {bam} | perl -ne 'print if !/^\@RG/' > {header} ''')
    keep_bam = []
    with Pool(nt) as p1:
        for g_name, cb_file in cb_files:
            sam_name = f'{tmp_dir}/{g_name}-sub.sam.gz'
            sam_header = f'{sam_name}.header'
            os.system(rf'''echo {g_name} | perl -nale 'print "\@RG\tID:$_\tSM:$_"' | cat {header} - | \
            bgzip > {sam_header} ''')
            
            os.system(rf'''samtools view -@ {nt} -O SAM --tag-file {tag_name}:{cb_file} -o {sam_name} {bam}''')
            index_bgzip(sam_name)
            regions = chunk_sam(sam_name, chunk_num=nt * 2)
            cmds = [rf'''bgzip -b {start} -s {q_len} {sam_name} |  \
                    mawk -F '\t' 'BEGIN{{OFS="\t"}}
                {{
                    if($5==255){{$5=60}}
                    rg = -1
                    for (i=12; i<= NF; i++){{
                        if(substr($i,1,3)=="RG:"){{
                            rg = i
                        }}
                    }}
                    if(rg==-1){{$(NF+1) = "RG:Z:{g_name}"}}else{{$rg="RG:Z:{g_name}"}}
                    print
                }}' | bgzip > {sam_name}-mod-{num}.gz''' for num, start, q_len in regions]
            
            p1.map(os.system, cmds, chunksize=1)
            sub_mod_files = [rf'''{sam_name}-mod-{num}.gz''' for num, *_ in regions]
            os.system(rf'''cat {sam_header} {' '.join(sub_mod_files)} | \
            samtools view --write-index -@ {nt} -O BAM -o {tmp_dir}/{g_name}-sub-mod.bam
            rm {' '.join(sub_mod_files)} {sam_header} {sam_name} {sam_name}.gzi {cb_file}''')
            keep_bam.append(f'{tmp_dir}/{g_name}-sub-mod.bam')
    os.system(f'rm {header}')
    
    # calling the PL
    with os.popen(f'tabix -l {snp_site}') as h1, Pool(nt) as p1:
        want_chrome = [i.strip() for i in h1]
        cmds = [rf'''tabix {snp_site} {chrome} | mawk '!a[$1" "$2]++' | bgzip > {tmp_dir}/site-{chrome}.gz
        tabix -s1 -b2 -e2 -C {tmp_dir}/site-{chrome}.gz''' for chrome in want_chrome]
        p1.map(os.system, cmds, chunksize=1)
        keep_size = {chrome: f'{tmp_dir}/site-{chrome}.gz' for chrome in want_chrome}
        
        keep_bcf = defaultdict(list)
        call_cmds = []
        for sub_bam in keep_bam:
            for chrome, site_file in keep_size.items():
                cmd1 = rf'''samtools view -O BAM {sub_bam} {chrome} | \
                bcftools mpileup -f {ref} -I -E -T {site_file} -Ob - | \
                bcftools call --write-index -Aim -C alleles -T {site_file} -Ob -o {sub_bam}-{chrome}.bcf'''
                call_cmds.append(cmd1)
                keep_bcf[chrome].append(f'{sub_bam}-{chrome}.bcf')
        p1.map(os.system, call_cmds, chunksize=1)
        
        merge_cmds = [rf'''bcftools merge -Ob -o {tmp_dir}/merge-call-{chrome}-pl.bcf --write-index {' '.join(bcf_list)}
        ''' for chrome, bcf_list in keep_bcf.items()]
        final_bcf_list = [f'{tmp_dir}/merge-call-{chrome}-pl.bcf' for chrome in keep_bcf]
        
        p1.map(os.system, merge_cmds, chunksize=1)
        
        os.system(rf'''rm {' '.join(f'{i} {i}.csi' for i in itertools.chain.from_iterable(keep_bcf.values()))}''')
        os.system(rf'''rm {' '.join(f'{i} {i}.csi' for i in keep_size.values())}''')
        final_bcf = f'{tmp_dir}/merge-all-pl.bcf'
        os.system(rf'''bcftools concat --write-index --threads {nt} -Ob \
            -o {final_bcf} {' '.join(final_bcf_list)}
            rm {' '.join(f'{i} {i}.csi' for i in final_bcf_list)}''')
    
    # run imputation
    with open(chunk) as h1:
        chunk_data = [i.split() for i in h1]
    imp_cmds = [rf'''GLIMPSE_phase --input {final_bcf} --reference {hap} \
        --input-region {row[2]} --output-region {row[3]} --output {tmp_dir}/chunk-{num}-imp.bcf
        bcftools index {tmp_dir}/chunk-{num}-imp.bcf''' for num, row in enumerate(chunk_data)]
    snp_num = [int(row[-1]) for row in chunk_data]
    order_cmd = [cmd for num, cmd in sorted(zip(snp_num, imp_cmds), reverse=True)]
    with Pool(nt) as p1:
        p1.map(os.system, order_cmd, chunksize=1)
        imp_bcf = [f'{tmp_dir}/chunk-{num}-imp.bcf' for num, *_ in enumerate(chunk_data)]
        fina_imp = f'{out_name}-imp-raw.bcf'
        os.system(rf'''bcftools concat --ligate --write-index \
            --threads {nt} -Ob -o {fina_imp} {' '.join(imp_bcf)}
            rm -R {tmp_dir}''')
    return fina_imp


def main(args=None):
    """Main entry point. If args is None, parse from command line."""
    parser = argparse.ArgumentParser(description='Imputed haplotype using groups of cell.')
    # 必需参数
    parser.add_argument('--bam', required=True, help='BAM file, should be sorted and index by samtools')
    parser.add_argument('--out-name', required=True, help='out file name')
    parser.add_argument('--ref', required=True, help='reference genome file')
    parser.add_argument('--snp-site', required=True, help='SNPs site file of the haplotype panel, indexed with tabix')
    parser.add_argument('--hap', required=True, help='haplotype panel file')
    parser.add_argument('--chunk', required=True, help='chunked region file')
    parser.add_argument('--cb-group', required=True, help='cell barcode group file')
    # 可选参数
    parser.add_argument('--nt', '--threads', type=int, default=16, help='number of threads to used, default: 16')
    if args is None:
        args = parser.parse_args()
    # 调用主要函数
    return get_sub(
        bam=args.bam,
        cb_group=args.cb_group,
        nt=args.nt,
        out_name=args.out_name,
        ref=args.ref,
        snp_site=args.snp_site,
        hap=args.hap,
        chunk=args.chunk
    )


if __name__ == '__main__':
    main()
