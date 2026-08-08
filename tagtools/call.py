#!/usr/bin/env python3
"""
TagTools Call Module - Call variants and perform genotyping

Usage: tagtools call [options]
"""

import os
import sys
import gzip
import random
import itertools
import argparse
from multiprocessing import Pool
from collections import defaultdict


def batched(iterable, n, *, strict=False):
    """Batch an iterable into chunks of size n"""
    if n < 1:
        raise ValueError('n must be at least one')
    iterator = iter(iterable)
    while batch := tuple(itertools.islice(iterator, n)):
        if strict and len(batch) != n:
            raise ValueError('batched(): incomplete batch')
        yield batch


def pairwise(iterable):
    """Return successive overlapping pairs from an iterable"""
    iterator = iter(iterable)
    a = next(iterator, None)
    for b in iterator:
        yield a, b
        a = b


def cmd_int(cmd):
    """Execute command and return integer output"""
    with os.popen(cmd) as h1:
        return int(list(h1)[0])


def uint16(data):
    """Convert bytes to uint16"""
    return np.frombuffer(data, dtype=np.uint16)[0]


def sam_header_size(sam):
    """Get size of SAM header in bytes"""
    with gzip.open(sam, 'rt') as h1:
        return sum(len(i) for i in itertools.takewhile(lambda x: x.startswith('@'), h1))


def index_bgzip(file, index_name=None):
    """
    Create index for bgzip compressed file
    
    Args:
        file: Path to bgzip compressed file
        index_name: Output path for index file (default: file.gzi)
    
    Returns:
        str: Path to the created index file
    """
    total = bytearray()
    with open(file, 'rb') as h1:
        while compressed_size := h1.read(18)[-2:]:
            fwd_bytes = int(uint16(compressed_size)) - 17
            uncompressed_size = h1.read(fwd_bytes)[-4:]
            total.extend(compressed_size + uncompressed_size)
    
    total = np.frombuffer(total, dtype=[('zip', np.uint16), ('unzip', np.uint32)]).copy()
    total['zip'] += 1
    
    if total['zip'].sum() != os.path.getsize(file):
        raise ValueError(f'{file} file format error.')
    
    bam_index = np.full(len(total) * 2 + 1, len(total), dtype=np.uint64)
    bam_index[1::2] = np.cumsum(total['zip'], dtype=np.uint64)
    bam_index[2::2] = np.cumsum(total['unzip'], dtype=np.uint64)
    
    if index_name is None:
        index_name = f'{file}.gzi'
    
    bam_index.tofile(index_name)
    return index_name


def chunk_bgz_txt(sam_file, *, chunk_num, start_pos=0):
    """Split bgzip file into chunks for parallel processing"""
    if not os.path.exists(f'{sam_file}.gzi'):
        raise FileNotFoundError(f'Bgzip index file {sam_file}.gzi not found.')
    
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
    """Accumulate cell barcode counts from chunk file"""
    with gzip.open(cb_count, 'rt') as h1:
        m1 = itertools.groupby((i.split() for i in h1), key=lambda x: x[0])
        return [(g_name, sum(int(i[-1]) for i in rows)) for g_name, rows in m1]


def tag_index_split(*, bam, tag_name, out_dir, nt, mod_mpq=False):
    """Index and split BAM file by tag"""
    if not os.path.exists(out_dir):
        raise FileNotFoundError(f'{out_dir} not found.')
    
    if mod_mpq:
        mod_mpq = rf'''mawk -F '\t' 'BEGIN{{OFS="\t"}} {{ if($5==255){{$5=60;}} print;}}' |'''
    else:
        mod_mpq = ''
    
    mid_sam = f'{out_dir}/tmp.sam.gz'
    os.system(rf'''samtools sort -@ {nt} -t {tag_name} -O SAM -T {mid_sam}-tmp -o {mid_sam} {bam}''')
    index_bgzip(mid_sam)
    header_size = sam_header_size(mid_sam)
    regions = chunk_bgz_txt(sam_file=mid_sam, chunk_num=nt * 2, start_pos=header_size)
    
    tag_pref_len = len(f'{tag_name}:Z:') + 1
    tag_len = len(tag_name)
    
    chunk_cal = [rf'''bgzip -b {start} -s {q_len} {mid_sam} | \
        mawk -F '\t' '{{
            for(i=12;i<=NF;i++){{
                if(substr($i, 1, {tag_len})=="{tag_name}"){{ tag=substr($i, {tag_pref_len});
                    print tag"\t"length+1; next;}}
            }} }}' | bgzip > {mid_sam}-{num}.wc.gz ''' for num, start, q_len in regions]
    
    cb_result_files = [f'{mid_sam}-{num}.wc.gz' for num, *_ in regions]
    
    with Pool(nt) as p1:
        p1.map(os.system, chunk_cal, chunksize=1)
        m1 = itertools.chain.from_iterable(p1.map(cb_acc, cb_result_files, chunksize=1))
    
    final = [(g_name, sum(int(i[-1]) for i in rows)) for g_name, rows in itertools.groupby(m1, key=lambda x: x[0])]
    cell_order = [i[0] for i in final]
    cell_count = {a: b for a, b in final}
    cell_len = [cell_count[i] for i in cell_order]
    c_index = [(start, end - start) for start, end in pairwise(itertools.accumulate(cell_len, initial=header_size))]
    c_index = [(tag, start, q_len) for tag, (start, q_len) in zip(cell_order, c_index)]
    
    header = f'{mid_sam}.header'
    os.system(rf'''samtools view -H {bam} > {header}''')
    
    cmd2 = [rf'''bgzip -b {start} -s {q_len} {mid_sam} | {mod_mpq}  cat {header} - | \
        samtools view -O BAM -o {out_dir}/{tag}.bam --write-index''' for tag, start, q_len in c_index]
    
    with Pool(nt) as p1:
        p1.map(os.system, cmd2, chunksize=1)
    
    os.system(rf'''rm {' '.join(cb_result_files)} {mid_sam} {mid_sam}.gzi {header}''')
    return [(tag, f'{out_dir}/{tag}.bam') for tag in cell_order]


def call_pl(*, bam, nt, out_name, ref, tag_name='CB', tag_list, site, chunk_size=20):
    """
    Call variants and perform genotyping
    
    Args:
        bam: Input BAM file path
        nt: Number of threads
        out_name: Output file name
        ref: Reference genome file
        tag_name: Tag name (default: CB)
        tag_list: Tag list file
        site: Site file
        chunk_size: Chunk size (default: 20)
    
    Returns:
        str: Output file path
    """
    out_dir = f'{out_name}-dir'
    if os.path.exists(out_dir):
        raise FileExistsError(f'{out_dir} exists.')
    
    mid_bam = f'{out_dir}/split.bam'
    new_size = f'{mid_bam}.site.gz'
    
    os.system(rf'''mkdir {out_dir}
        samtools view -O BAM -o {mid_bam} -@ {nt} --tag-file {tag_name}:{tag_list} -L {site} {bam}
        samtools depth -@ {nt} -b {site} {mid_bam} | mawk '$NF>0{{print $1"\t"$2}}' > {new_size}.pos
        zcat {site} | mawk 'FNR==NR{{a[$1" "$2]++; next}} $1" "$2 in a' {new_size}.pos  - | \
        mawk '!a[$1" "$2]++'| sort -k1,1 -k2,2n | bgzip > {new_size} && tabix -s1 -b2 -e2 -C {new_size}
        gzip -dc {new_size} | mawk '{{print NR"\t"$0}}' | bgzip > {new_size}-num.gz
        tabix -s1 -b2 -e2 -C {new_size}-num.gz''')
    
    cell_bam = tag_index_split(bam=mid_bam, tag_name=tag_name, out_dir=out_dir, nt=nt, mod_mpq=True)
    random.shuffle(cell_bam)
    keep_cell = list(batched([i[0] for i in cell_bam], chunk_size))
    keep_bam = list(batched([i[1] for i in cell_bam], chunk_size))
    
    call_pl_cmds = [rf'''samtools depth -b {new_size} {' '.join(bam_list)} | \
        mawk '{{for (i=3; i<=NF; i++) if($i > 0){{print $1"\t"$2; break}} }}' | tabix -R - {new_size} | \
        bgzip > {mid_bam}-chunk{num}.site.gz && tabix -s1 -b2 -e2 -C {mid_bam}-chunk{num}.site.gz
        bcftools mpileup --ignore-RG -f {ref} -I -E -T {mid_bam}-chunk{num}.site.gz -Ob {' '.join(bam_list)} | \
        bcftools call -Am -C alleles -T {mid_bam}-chunk{num}.site.gz -Ob | \
        bcftools view -e 'INFO/DP<1' -Ob | bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\t[%GT:%PL\t]\n" - | \
        mawk '
            BEGIN {{
                sample_ids_str = "{' '.join(cell_list)}"
                split(sample_ids_str, sample_ids, " ")
                sample_count = length(sample_ids)
            }}
            {{
                for (i = 1; i <= sample_count; i++) {{
                    pl = $(i+4)
                    if(substr(pl,1,1) !="."){{
                        split(substr(pl,5), nums, ",")
                        print sample_ids[i]"\t"$1"\t"$2"\t"$3"\t"$4"\t"nums[1]"\t"nums[2]"\t"nums[3]
                    }}
                }}
            }}' | sort -k1,1 -k2,2 -k3,3n | mawk '{{print NR"\t"$0}}' | bgzip > {mid_bam}-chunk{num}.col.pl.gz
        ''' for num, (cell_list, bam_list) in enumerate(zip(keep_cell, keep_bam))]
    
    keep_pl = [rf'''{mid_bam}-chunk{num}.col.pl.gz''' for num, _ in enumerate(keep_cell)]
    
    with Pool(nt) as p1:
        p1.map(os.system, call_pl_cmds, chunksize=1)
    
    os.system(rf'''cat {' '.join(keep_pl)}  > {out_name} && tabix -s2 -b1 -e1 -C {out_name}''')
    os.system(f'rm -R {out_dir}')
    return out_name


def main():
    """Command-line entry point for call module"""
    parser = argparse.ArgumentParser(
        prog='tagtools call',
        description='Call variants and perform genotyping'
    )
    parser.add_argument('--bam', required=True, help='Input BAM file path')
    parser.add_argument('--nt', type=int, default=4, help='Number of threads (default: 4)')
    parser.add_argument('--out-name', required=True, help='Output file name')
    parser.add_argument('--ref', required=True, help='Reference genome file')
    parser.add_argument('--tag-name', default='CB', help='Tag name (default: CB)')
    parser.add_argument('--tag-list', required=True, help='Tag list file')
    parser.add_argument('--site', required=True, help='Site file')
    parser.add_argument('--chunk-size', type=int, default=20, help='Chunk size (default: 20)')
    
    args = parser.parse_args()
    
    try:
        result = call_pl(
            bam=args.bam,
            nt=args.nt,
            out_name=args.out_name,
            ref=args.ref,
            tag_name=args.tag_name,
            tag_list=args.tag_list,
            site=args.site,
            chunk_size=args.chunk_size
        )
        print(f'Variant calling complete. Output file: {result}')
    except Exception as e:
        print(f'Error: {e}', file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
