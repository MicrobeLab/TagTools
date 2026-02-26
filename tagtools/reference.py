import os
import numpy as np
import argparse
from multiprocessing import Pool


def pairwise(iterable):
    iterator = iter(iterable)
    a = next(iterator, None)
    for b in iterator:
        yield a, b
        a = b


def cmd_int(cmd):
    with os.popen(cmd) as h1:
        return int(list(h1)[0])


def index_bgzip(file):
    total = bytearray()
    with open(file, 'rb') as h1:
        while tmp1 := h1.read(18)[-2:]:
            block_size = int(np.frombuffer(tmp1, dtype=np.uint16)[0]) - 17
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


def chunk_bgzip(file, *, chunk_num, start_pos=0):
    total = int(np.fromfile(f'{file}.gzi', dtype=np.uint64)[-1])
    chunk_size = int((total - start_pos) / chunk_num) + 1
    pos = [start_pos]
    while start_pos < total:
        start_pos += chunk_size
        if start_pos >= total:
            pos.append(total)
            break
        start_pos += cmd_int(f'bgzip -b {start_pos} {file} | head -1 |wc -c')
        pos.append(start_pos)
    return [(num, start, end - start) for num, (start, end) in enumerate(pairwise(pos))]


def create_hap_panel(*, vcf_list, out_name, nt, add_chr=False, rm_chr=False):
    """
    :param vcf_list: str, path to a text file containing a list of bgzip-compressed haplotype VCF files,
                      with one file path per line.
    :param out_name: str, prefix for the output file names.
    :param nt: int, number of threads to use for processing.
    :param add_chr: bool, whether to add the "chr" prefix to the haplotype panel.
    :param rm_chr: bool, whether to remove the "chr" prefix from the haplotype panel.
    :return: list, containing:
        - str: path to the generated haplotype panel BCF file.
        - str: path to the SNP site file.
        - str: path to the chunked genomic region file.
    """
    if not add_chr and not rm_chr:
        mod_chr = ''
    elif add_chr:
        mod_chr = rf''' | mawk '{{print "chr"$0}}' '''
    elif rm_chr:
        mod_chr = rf'''| mawk '{{tmp=substr($0,4); print tmp}}' '''
    else:
        raise RuntimeError(f'add_chr and rm_chr both True.')
    
    tmp_name = f'{out_name}-tmp.gz'
    with open(vcf_list) as h1:
        data = [i.split()[0] for i in h1]
    os.system(rf'''cat {' '.join(data)} > {tmp_name}''')
    index_bgzip(tmp_name)
    regions = chunk_bgzip(tmp_name, chunk_num=nt * 2)
    header = f'{tmp_name}.header'
    os.system(rf'''bcftools view -h {data[0]} | perl -nale 'print if !/^##contig=<ID/' > {header} ''')
    
    make_hap_cmds = [rf'''bgzip -b {start} -s {q_len} {tmp_name} |mawk -F '\t' '
        BEGIN {{
            a["G"] = a["C"] = a["T"] = a["A"] = 1;
            OFS = "\t";
        }}
        (substr($0, 1, 1) != "#") && ($4 in a) && ($5 in a) {{
            if ($3 == ".") {{
                $3 = $1 "-" $2 ":" $4 ":" $5;
            }}
            print;
        }}' {mod_chr} | cat {header} - | bgzip > {tmp_name}-chunk{num}.vcf.gz && tabix -C {tmp_name}-chunk{num}.vcf.gz
    bcftools view -Ob -o {tmp_name}-chunk{num}.bcf {tmp_name}-chunk{num}.vcf.gz
    bcftools view -G -Ob -o {tmp_name}-chunk{num}.site.bcf {tmp_name}-chunk{num}.bcf
    bcftools view -q 0.05:minor -Ob {tmp_name}-chunk{num}.site.bcf | \
    bcftools query -f "%CHROM\t%POS\t%REF,%ALT\n" | bgzip > {tmp_name}-chunk{num}.MAF0.05.site.gz
    bcftools query -f "%CHROM\t%POS\t%REF,%ALT\n" {tmp_name}-chunk{num}.site.bcf | \
    bgzip > {tmp_name}-chunk{num}.site.gz
    rm {tmp_name}-chunk{num}.vcf.gz {tmp_name}-chunk{num}.vcf.gz.csi
    ''' for num, start, q_len in regions]
    keep_bcf = [f'{tmp_name}-chunk{num}.bcf' for num, *_ in regions]
    site_bcf = [f'{tmp_name}-chunk{num}.site.bcf' for num, *_ in regions]
    site_tsv = [f'{tmp_name}-chunk{num}.site.gz' for num, *_ in regions]
    common_site = [f'{tmp_name}-chunk{num}.MAF0.05.site.gz' for num, *_ in regions]
    
    with Pool(nt) as p1:
        p1.map(os.system, make_hap_cmds, chunksize=1)
    
    os.system(rf'''bcftools concat --threads {nt} -Ob -o {out_name}.bcf --write-index {' '.join(keep_bcf)}
    bcftools concat --threads {nt} -Ob -o {out_name}.site.bcf --write-index {' '.join(site_bcf)}
    rm {' '.join(site_bcf)} {' '.join(keep_bcf)}

    cat {' '.join(site_tsv)}  > {out_name}.site.gz
    cat {' '.join(common_site)}  > {out_name}.MAF0.05.site.gz
    tabix -s1 -b2 -e2 -C {out_name}.MAF0.05.site.gz
    tabix -s1 -b2 -e2 -C {out_name}.site.gz && rm {' '.join(site_tsv + common_site)}
    rm {tmp_name} {tmp_name}.gzi {header}''')
    
    # make the chunk file for GLIMPSE imputation
    with os.popen(f'tabix -l {out_name}.site.gz') as h1:
        keep_chr = [i.strip() for i in h1]
    chunk_cmds = [rf'''GLIMPSE_chunk --input {out_name}.site.bcf --region {chrome} --window-size 2000000 \
    --buffer-size 200000 --output {out_name}-chunk{chrome}-region.txt''' for chrome in keep_chr]
    keep_chunk = [f'{out_name}-chunk{chrome}-region.txt' for chrome in keep_chr]
    
    with Pool(nt) as p:
        p.map(os.system, chunk_cmds, chunksize=1)
    os.system(rf'''cat {' '.join(keep_chunk)} > {out_name}.region && rm {' '.join(keep_chunk)}
    rm {out_name}.site.bcf {out_name}.site.bcf.csi''')
    return f'{out_name}.bcf {out_name}.site.gz {out_name}.region'.split()


def main(args=None):
    """Main entry point. If args is None, parse from command line."""
    parser = argparse.ArgumentParser(description="Create haplotype panel from VCF files.")
    parser.add_argument("--vcf-list", type=str, required=True,
                        help="Path to a text file containing a list of bgzip-compressed haplotype VCF files.")
    parser.add_argument("--out-name", type=str, required=True, help="Prefix for the output file names.")
    parser.add_argument("--nt", type=int, required=True, help="Number of threads to use for processing.")
    parser.add_argument("--add-chr", action="store_true",
                        help="Whether to add the 'chr' prefix to the haplotype panel.")
    parser.add_argument("--rm-chr", action="store_true",
                        help="Whether to remove the 'chr' prefix from the haplotype panel.")
    if args is None:
        args = parser.parse_args()
    return create_hap_panel(vcf_list=args.vcf_list, out_name=args.out_name, nt=args.nt, add_chr=args.add_chr,
                     rm_chr=args.rm_chr)


if __name__ == '__main__':
    main()
