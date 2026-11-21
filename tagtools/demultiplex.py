import os
import argparse
import numpy as np
import random
import itertools
import sys
import gzip
from collections import defaultdict
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


def cb_acc(cb_count):
    with gzip.open(cb_count, 'rt') as h1:
        m1 = itertools.groupby((i.split() for i in h1), key=lambda x: x[0])
        return [(g_name, sum(int(i[-1]) for i in rows)) for g_name, rows in m1]


def split_by_tag(*, sort_tag_sam, nt, tag_name='CB', tag_type='Z'):
    index_name = index_bgzip(sort_tag_sam)
    regions = chunk_sam(sam_file=sort_tag_sam, chunk_num=nt * 2)
    header = f'{sort_tag_sam}.header'
    os.system(rf'''samtools view -H {sort_tag_sam} > {header}''')
    pre_len = len(f'{tag_name}:{tag_type}:') + 1
    tag_len = len(tag_name)
    cmd1 = [rf'''bgzip -b {start} -s {q_len} {sort_tag_sam} | \
        mawk -F '\t' '{{
            for(i=12;i<=NF;i++){{
                if(substr($i, 1, {tag_len})=="{tag_name}"){{ tag=substr($i, {pre_len}); print tag"\t"length+1; next;}}
            }} }}' | bgzip > {sort_tag_sam}-{num}.wc.gz ''' for num, start, q_len in regions]
    
    header_size = get_header_size(sort_tag_sam)
    with Pool(nt) as p1:
        p1.map(os.system, cmd1, chunksize=1)
        cb_result_files = [f'{sort_tag_sam}-{num}.wc.gz' for num, *_ in regions]
        final = list(itertools.chain.from_iterable(p1.map(cb_acc, cb_result_files, chunksize=1)))
        cell_order = []
        cell_count = defaultdict(int)
        for cell, count in final:
            cell_count[cell] += count
            if cell not in cell_order:
                cell_order.append(cell)
        cell_len = [cell_count[i] for i in cell_order]
        c_index = [(start, end - start) for start, end in pairwise(itertools.accumulate(cell_len, initial=header_size))]
        c_index = [(cell, start, q_len) for cell, (start, q_len) in zip(cell_order, c_index)]
        cmd2 = [rf'''bgzip -b {start} -s {q_len} {sort_tag_sam} | cat {header} - | \
            samtools view -O BAM -o {sort_tag_sam}-{cell}.bam
            samtools index {sort_tag_sam}-{cell}.bam''' for cell, start, q_len in c_index]
        p1.map(os.system, cmd2, chunksize=1)
    os.system(rf'''rm {' '.join(cb_result_files)} {index_name} {header}''')
    return [(cell, f'{sort_tag_sam}-{cell}.bam') for cell in cell_order]


def split_vcf(*, vcf, bam, nt, out_dir, cb_list, tag_name='CB', info=0.4, tag_type='Z'):
    if tag_name == 'RG':  # the RG tag used for group-name in split BAM by group of CB
        raise ValueError(f'{tag_name} should not be "RG".')
    
    if os.path.exists(out_dir):
        raise FileExistsError(f'{out_dir} exist.')
    os.system(f'mkdir {out_dir}')
    out_name = f'{out_dir}/demuxlet'
    
    if not os.path.exists(f'{bam}.bai') and not os.path.exists(f'{bam}.csi'):
        raise FileNotFoundError(f'{bam}.csi or {bam}.bai not found.')
    
    cb_group = f'{out_name}-cb-group.txt'
    cb_d1 = defaultdict(list)
    with open(cb_group, 'wt') as h1, open(cb_list) as h2:
        cb_data = [i.split()[0] for i in h2]
        split_cb_list_size = int(len(cb_data) / (2 * nt)) + 1  # split the CB list into group, the split chunk size
        random.shuffle(cb_data)
        for num, chunk_items in enumerate(batched(cb_data, split_cb_list_size), 1):
            for item in chunk_items:
                g_name = f'chunk{num}'
                print(item, g_name, file=h1, sep='\t')
                cb_d1[g_name].append(item)  # each group-name append cell barcode into list
    
    filter1_vcf = f'{out_name}-maf-info.vcf.gz'
    filter1_pos = f'{filter1_vcf}.pos'
    mid_sam = f'{out_name}.sam.gz'
    header = f'{mid_sam}.header.gz'
    cmd1 = rf'''bcftools view --threads {nt} -e 'INFO/INFO<{info}' -Oz -o {filter1_vcf} -q 0.001:minor {vcf}
    bcftools query -f '%CHROM\t%POS\n' {filter1_vcf}  > {filter1_pos}
    samtools depth -@ {nt} -b {filter1_pos} {bam} | \
    perl -nale 'print "$F[0]\t$F[1]" if $F[-1]>0' > {filter1_pos}-1 && mv {filter1_pos}-1 {filter1_pos}
    gzip -dc {filter1_vcf} | mawk 'FNR==NR{{a[$0]=1; next}}
        /^#/{{print; next}} $1"\t"$2 in a' {filter1_vcf}.pos - | \
        bgzip > {filter1_vcf}-1 && mv {filter1_vcf}-1 {filter1_vcf}
    samtools view -@ {nt} -L {filter1_pos} --tag-file {tag_name}:{cb_list} -O SAM -o {mid_sam} {bam}
    samtools view -H {bam} | bgzip > {header}
    '''
    os.system(cmd1)
    
    # add CB to RG tag and sort by RG tag to sam
    index_file = index_bgzip(mid_sam)
    regions = chunk_sam(mid_sam, chunk_num=nt * 2)
    tag_len = len(tag_name)
    pre_len = len(f'{tag_name}:{tag_type}:') + 1
    add_rg_cmds = [rf'''bgzip -b {start} -s {q_len} {mid_sam} |  \
    mawk -F '\t' 'BEGIN{{OFS="\t"}}
    FNR==NR{{ a[$1]=$2; next }}
    {{
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
    
    with Pool(nt) as p1:
        p1.map(os.system, add_rg_cmds, chunksize=1)
    keep_chunk_files = [rf'''{mid_sam}-mod-{num}.gz''' for num, *_ in regions]
    sorted_sam = f'{mid_sam}-sort.sam.gz'
    os.system(rf'''cat {header} {' '.join(keep_chunk_files)} | \
        samtools sort -@ {nt} -O SAM -o {sorted_sam} -T {sorted_sam}-tmp -t RG && rm {' '.join(keep_chunk_files)}''')
    
    cell_bam = split_by_tag(sort_tag_sam=sorted_sam, nt=nt, tag_type='Z', tag_name='RG')  # split by RG
    
    cb_group_files = {}
    for g_name, vals in cb_d1.items():
        tmp_name = f'{mid_sam}-{g_name}.txt'
        cb_group_files[g_name] = tmp_name
        with open(tmp_name, 'wt') as h1:
            for item in vals:
                print(item, file=h1)
    
    with Pool(nt) as p1:
        dem_cmds = [rf'''
        samtools depth -b {filter1_pos} {sub_bam} | perl -nale 'print "$F[0]\t$F[1]" if $F[-1]>0' > {sub_bam}.pos
        gzip -dc {filter1_vcf} | mawk 'FNR==NR{{a[$0]=1; next}} /^#/{{print; next}} $1"\t"$2 in a' {sub_bam}.pos - | \
        bgzip  > {sub_bam}.vcf.gz && bcftools index -t {sub_bam}.vcf.gz
        demuxlet --sam {sub_bam} --vcf {sub_bam}.vcf.gz --group-list {cb_group_files[g_name]} \
            --out {sub_bam}-{g_name}-predict
        rm {sub_bam}.pos {sub_bam}.vcf.gz {sub_bam}.vcf.gz.tbi {sub_bam} {sub_bam}.bai'''
                    for g_name, sub_bam in cell_bam]
        p1.map(os.system, dem_cmds, chunksize=1)
        
        # remove the split CB group files, and the filter pos file
        os.system(rf'''rm {' '.join(cb_group_files.values())} {filter1_pos}''')
        
        # remove the mid same (and index file), sorted sam, header, cell barcode group file
        os.system(rf'''rm {mid_sam} {index_file} {sorted_sam} {header} {cb_group}''')
        
        want_best = [f'{sub_bam}-{g_name}-predict.best' for g_name, sub_bam in cell_bam if
                     os.path.exists(f'{sub_bam}-{g_name}-predict.best')]
        want_single = [f'{sub_bam}-{g_name}-predict.single' for g_name, sub_bam in cell_bam if
                       os.path.exists(f'{sub_bam}-{g_name}-predict.single')]
        want_sing2 = [f'{sub_bam}-{g_name}-predict.sing2' for g_name, sub_bam in cell_bam if
                      os.path.exists(f'{sub_bam}-{g_name}-predict.sing2')]
    os.system(rf'''cat {' '.join(want_best)} | perl -nale 'print && next if $.==1;
        print if $F[0] ne "BARCODE"' > {out_name}.best && rm {' '.join(want_best)}''')
    os.system(rf'''cat {' '.join(want_single)} | perl -nale 'print && next if $.==1;
        print if $F[0] ne "BARCODE"' > {out_name}.single  && rm {' '.join(want_single)}''')
    os.system(rf'''cat {' '.join(want_sing2)} | perl -nale 'print && next if $.==1;
            print if $F[0] ne "BARCODE"' > {out_name}.sing2  && rm {' '.join(want_sing2)}''')


def main(args=None):
    """Main entry point. If args is None, parse from command line."""
    parser = argparse.ArgumentParser(description='Demultiplexing using demuxlet.')
    # 必需参数
    parser.add_argument('--vcf', required=True, help='imputed VCF/BCF file, indexed with bcftools')
    parser.add_argument('--bam', required=True, help='BAM file, sorted and index with samtools')
    parser.add_argument('--cb-list', required=True, help='cell barcode list file')
    parser.add_argument('--out-dir', required=True, help='out directory name')
    # 可选参数
    parser.add_argument('--nt', '--threads', type=int, default=16, help='number of threads to used, default: 16')
    parser.add_argument('--info', type=float, default=0.4, help='INFO score to filter imputed genotype, default: 0.4')
    if args is None:
        args = parser.parse_args()
    # 调用主要函数
    return split_vcf(
        vcf=args.vcf,
        bam=args.bam,
        out_dir=args.out_dir,
        nt=args.nt,
        cb_list=args.cb_list,
        info=args.info
    )


if __name__ == '__main__':
    main()
