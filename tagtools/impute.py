import gzip
import os
import argparse
import numpy as np
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
    return index_name


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


def get_sub(*, bam, cb_group, nt, out_name, tag_name='CB', snp_site, ref, chunk, hap, tag_type='Z'):
    tmp_dir = f'{out_name}-dir'
    if os.path.exists(tmp_dir):
        raise FileExistsError(f'{tmp_dir} exits.')
    os.system(f'mkdir {tmp_dir}')

    # check the bam index
    if not os.path.exists(f'{bam}.csi') and not os.path.exists(f'{bam}.bai'):
        raise FileNotFoundError(f'{bam}.csi or {bam}.bai not found, {bam} should sorted and indexed.')

    with open(chunk) as h1:
        chunk_data = [i.split() for i in h1]
        snp_num = [int(row[-1]) for row in chunk_data]  # num, chrome, in-region, out-region *_, snp-number
        chrome_set = {i[1] for i in chunk_data}  # keep-set

    with os.popen(rf'samtools idxstats {bam}') as h1:
        order_chrome = []
        for row in h1:
            parts = row.split()
            if parts[0] in chrome_set and parts[0] not in order_chrome:
                order_chrome.append(parts[0])  # unique and order chrome list

        order_chunk_nums = []
        for fix_chrome in order_chrome:
            tmp_chunk_nums = [i[0] for i in chunk_data if i[1] == fix_chrome]
            order_chunk_nums.extend(tmp_chunk_nums)  # imputed chunked bcf files will order as the listed

    with open(cb_group) as h1:
        data = [i.split() for i in h1]
        if any(len(i) != 2 for i in data):
            raise RuntimeError(f'{cb_group} requires with two columns.')

    header = f'{tmp_dir}/pub-header.txt.gz'
    mid_sam = f'{header}.sam.gz'
    os.system(rf'''samtools view -H {bam} | perl -ne 'print if !/^\@RG/' | bgzip > {header}
    perl -nale 'print $F[-1]' {cb_group} | sort | uniq | perl -nale 'print "\@RG\tID:$_\tSM:$_"' | bgzip >> {header}
    perl -nale 'print $F[0]' {cb_group} > {header}.barcode
    samtools view -@ {nt} --tag-file {tag_name}:{header}.barcode -O SAM -o {mid_sam} {bam}
    ''')
    index_bgzip(mid_sam)
    regions = chunk_sam(mid_sam, chunk_num=nt * 2)
    tag_len = len(tag_name)
    pre_len = len(f'{tag_name}:{tag_type}:') + 1
    cmds = [rf'''bgzip -b {start} -s {q_len} {mid_sam} |  \
    mawk -F '\t' 'BEGIN{{OFS="\t"}}
    FNR==NR{{ a[$1]=$2; next }}
    {{
        if($5==255){{ $5=60 }}
        rg = -1
        grp = ""
        for (i=12; i<=NF; i++) {{
            if (substr($i,1,3)=="RG:") {{ rg = i }}
            if (substr($i,1,{tag_len})=="{tag_name}") {{
                cb1 = substr($i,{pre_len})
                if (cb1 in a) {{ grp = a[cb1] }}
            }}
        }}
        if (rg != -1 && grp != "") {{
            $rg = "RG:Z:" grp
            print
        }} else if (rg == -1 && grp != "") {{
            $(NF+1) = "RG:Z:" grp
            print
        }}
    }}' {cb_group} - | bgzip > {mid_sam}-mod-{num}.gz''' for num, start, q_len in regions]
    mid_bam = f'{mid_sam}-mod.bam'
    with Pool(nt) as p1:
        p1.map(os.system, cmds, chunksize=1)
    sub_mod_files = [rf'''{mid_sam}-mod-{num}.gz''' for num, *_ in regions]
    os.system(rf'''cat {header} {' '.join(sub_mod_files)} | samtools view --write-index -@ {nt} -O BAM -o {mid_bam}''')

    # calling PL
    with open(cb_group) as h1:
        num_group = len({i[-1] for i in (line.split() for line in h1)})  # how many group
        call_depth = num_group * 100  # each group set 100x read depth

    call_cmds = [rf'''tabix {snp_site} {out_region} | bgzip > {mid_bam}-chunk{num}.gz
    samtools view -O BAM {mid_bam} {in_region} | \
    bcftools mpileup --max-depth {call_depth} -f {ref} -I -E -T {mid_bam}-chunk{num}.gz -Ob - | \
    bcftools call --write-index -Aim -C alleles -T {mid_bam}-chunk{num}.gz -Ob -o {mid_bam}-chunk{num}.bcf
    ''' for num, chrome, in_region, out_region, *_ in chunk_data]
    call_cmds = [cmd for num, cmd in sorted(zip(snp_num, call_cmds), reverse=True)]  # sort by SNP number
    with Pool(30) as p1:
        p1.map(os.system, call_cmds, chunksize=1)

    want_pl = [rf'''{mid_bam}-chunk{num}.bcf''' for num, *_ in chunk_data]
    os.system(rf'''bcftools concat --threads {nt} --write-index -Ob -o {mid_bam}-pl.bcf {' '.join(want_pl)}''')

    imp_cmds = [rf'''GLIMPSE_phase --input {mid_bam}-pl.bcf --reference {hap} \
        --input-region {in_region} --output-region {out_region} --output {mid_bam}-{num}-imp.bcf
        bcftools index {mid_bam}-{num}-imp.bcf''' for num, chrome, in_region, out_region, *_ in chunk_data]

    imp_cmds = [cmd for num, cmd in sorted(zip(snp_num, imp_cmds), reverse=True)]
    with Pool(nt) as p1:
        p1.map(os.system, imp_cmds, chunksize=1)
        imp_bcf = {num: f'{mid_bam}-{num}-imp.bcf' for num, *_ in chunk_data}
        imp_bcf = [imp_bcf[i] for i in order_chunk_nums]  # re-order the imputed bcf keep chrome as bam order
        fina_imp = f'{out_name}-imp-raw.bcf'
    os.system(rf'''bcftools concat --ligate --write-index --threads {nt} -Ob -o {fina_imp} {' '.join(imp_bcf)}
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
