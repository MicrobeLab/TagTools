import os
import argparse
import random
import itertools
import sys
from multiprocessing import Pool


def batched(iterable, n, *, strict=False):
    if n < 1:
        raise ValueError('n must be at least one')
    iterator = iter(iterable)
    while batch := tuple(itertools.islice(iterator, n)):
        if strict and len(batch) != n:
            raise ValueError('batched(): incomplete batch')
        yield batch


def split_vcf(*, vcf, bam, nt, out_dir, cb_list, tag_name='CB', info=0.4):
    if os.path.exists(out_dir):
        raise FileExistsError(f'{out_dir} exist.')
    os.system(f'mkdir {out_dir}')
    out_name = f'{out_dir}/demuxlet'
    os.system(rf'''samtools idxstats {bam} > {out_name}-index.txt
    tabix -l {vcf} | mawk 'FNR==NR{{a[$1]=1; next}} $1 in a' - {out_name}-index.txt | \
    perl -nale 'print $F[0]' > {out_name}-index.txt1 && mv {out_name}-index.txt1 {out_name}-index.txt''')
    with open(f'{out_name}-index.txt') as h1:
        want_chrome = [i.strip() for i in h1]
        if len(want_chrome) < 1:
            raise RuntimeError(f'No share chromosomes with {vcf}  and {bam}.')
    
    for chrome in want_chrome:
        cmd1 = rf'''bcftools view --threads {nt} -e 'INFO/INFO<{info}' -q 0.001:minor {vcf} {chrome} | \
        perl -ne 'next if /^##contig=/; print' | bgzip -@ {nt} > {out_name}-{chrome}.vcf.gz
        tabix -s1 -b2 -e2 -C {out_name}-{chrome}.vcf.gz
        bcftools view --threads {nt} -Ob --write-index -o {out_name}-{chrome}.bcf {out_name}-{chrome}.vcf.gz
        rm {out_name}-{chrome}.vcf.gz {out_name}-{chrome}.vcf.gz.csi'''
        os.system(cmd1)
    want_vcf = [rf'''{out_name}-{chrome}.bcf''' for chrome in want_chrome]
    os.system(rf'''bcftools concat -Ob -o {out_name}-merge.bcf --write-index {' '.join(want_vcf)}
    bcftools view -Oz -o {out_name}-merge.vcf.gz --threads {nt} {out_name}-merge.bcf
    rm {' '.join(f'{i} {i}.csi' for i in want_vcf)}
    rm {out_name}-index.txt
    bcftools query -f '%CHROM\t%POS\t%REF,%ALT\n' {out_name}-merge.bcf | bgzip > {out_name}-merge.site.gz
    rm {out_name}-merge.bcf {out_name}-merge.bcf.csi
    tabix -s1 -b2 -e2 -C {out_name}-merge.site.gz
    samtools depth -@ {nt} -b {out_name}-merge.site.gz {bam} | \
    perl -nale 'print "$F[0]\t$F[1]" if $F[-1]>0' | tabix -R - {out_name}-merge.site.gz | \
    bgzip > {out_name}-merge.site.gz-1
    mv {out_name}-merge.site.gz-1 {out_name}-merge.site.gz
    tabix -f -s1 -b2 -e2 -C {out_name}-merge.site.gz
    gzip -dc {out_name}-merge.site.gz |perl -nale 'print "$F[0]\t$F[1]"' > {out_name}-merge.site.gz.pos
    gzip -dc {out_name}-merge.vcf.gz | mawk 'FNR==NR{{a[$0]=1; next}}
        /^#/{{print; next}} $1"\t"$2 in a' {out_name}-merge.site.gz.pos - | bgzip > {out_name}-merge.vcf.gz-1
    mv {out_name}-merge.vcf.gz-1 {out_name}-merge.vcf.gz && rm {out_name}-merge.site.gz.pos''')
    
    filter_site = f'{out_name}-merge.site.gz'
    filter_vcf = f'{out_name}-merge.vcf.gz'
    sub_bam_list = []
    for chrome in want_chrome:
        cmd1 = rf'''tabix {filter_site} {chrome} | bgzip > {filter_site}-{chrome}.gz
        samtools view --write-index -@ {nt} -O BAM -L {filter_site}-{chrome}.gz -o {out_name}-{chrome}-sub.bam {bam}
        rm {filter_site}-{chrome}.gz'''
        os.system(cmd1)
        sub_bam_list.append(f'{out_name}-{chrome}-sub.bam')
    filter_bam = f'{out_name}-filter-region.bam'
    os.system(rf'''samtools cat -@ {nt} -o {filter_bam} {' '.join(sub_bam_list)}
    rm {' '.join(f'{i} {i}.csi' for i in sub_bam_list)}
    samtools index -@ {nt} {filter_bam}''')
    
    # make the demuxlet input files
    with open(cb_list) as h1:
        all_cb = [i.split()[0] for i in h1]
        random.shuffle(all_cb)
    chunk_size = int(len(all_cb) / nt) + 1
    keep_args_demuxlet = []
    for num, chunk_item in enumerate(batched(all_cb, chunk_size)):
        sub_cb = f'{out_name}-chunk{num}-cb.txt'
        with open(sub_cb, 'wt') as h1:
            for item in chunk_item:
                print(item, file=h1)
        sub_bam = f'{out_name}-chunk{num}-cb.bam'
        os.system(rf'''samtools view -O BAM -@ {nt} --tag-file {tag_name}:{sub_cb} -o {sub_bam} {filter_bam}
        samtools index -@ {nt} {sub_bam}
        samtools depth -@ {nt} -b {filter_site} {sub_bam} | \
        perl -nale 'print "$F[0]\t$F[1]" if $F[-1]>0' > {sub_bam}.pos
        gzip -dc {filter_vcf} | mawk 'FNR==NR{{a[$0]=1; next}} /^#/{{print; next}} $1"\t"$2 in a' {sub_bam}.pos - | \
        bgzip -@ {nt} > {sub_bam}.vcf.gz && bcftools index --threads {nt} -t {sub_bam}.vcf.gz
        rm {sub_bam}.pos''')
        keep_args_demuxlet.append([sub_cb, sub_bam, f'{sub_bam}.vcf.gz'])
    
    with Pool(nt) as p1:
        dem_cmds = [rf'''demuxlet --sam {sub_bam} --vcf {sub_vcf} --group-list {sub_cb} --out {sub_bam}-predict
            rm {sub_bam} {sub_bam}.bai {sub_cb} {sub_vcf} {sub_vcf}.tbi'''
                    for sub_cb, sub_bam, sub_vcf in keep_args_demuxlet]
        p1.map(os.system, dem_cmds, chunksize=1)
        want_best = [f'{sub_bam}-predict.best' for _, sub_bam, _ in keep_args_demuxlet if
                     os.path.exists(f'{sub_bam}-predict.best')]
        want_single = [f'{sub_bam}-predict.single' for _, sub_bam, _ in keep_args_demuxlet if
                       os.path.exists(f'{sub_bam}-predict.single')]
        want_sing2 = [f'{sub_bam}-predict.sing2' for _, sub_bam, _ in keep_args_demuxlet if
                      os.path.exists(f'{sub_bam}-predict.sing2')]
    os.system(rf'''cat {' '.join(want_best)} | perl -nale 'print && next if $.==1;
        print if $F[0] ne "BARCODE"' > {out_name}.best && rm {' '.join(want_best)}
        rm {filter_bam} {filter_bam}.bai {filter_site} {filter_site}.csi''')
    os.system(rf'''cat {' '.join(want_single)} | perl -nale 'print && next if $.==1;
        print if $F[0] ne "BARCODE"' > {out_name}.single  && rm {' '.join(want_single)}''')
    os.system(rf'''cat {' '.join(want_sing2)} | perl -nale 'print && next if $.==1;
            print if $F[0] ne "BARCODE"' > {out_name}.sing2  && rm {' '.join(want_sing2)}''')


def main():
    parser = argparse.ArgumentParser(description='Demultiplexing using demuxlet.')
    # 必需参数
    parser.add_argument('--vcf', required=True, help='imputed VCF/BCF file, indexed with bcftools')
    parser.add_argument('--bam', required=True, help='BAM file, sorted and index with samtools')
    parser.add_argument('--cb-list', required=True, help='cell barcode list file')
    parser.add_argument('--out-dir', required=True, help='out directory name')
    # 可选参数
    parser.add_argument('--nt', '--threads', type=int, default=16, help='number of threads to used, default: 16')
    parser.add_argument('--info', type=float, default=0.4, help='INFO score to filter imputed genotype, default: 0.4')
    args = parser.parse_args()
    # 调用主要函数
    split_vcf(
        vcf=args.vcf,
        bam=args.bam,
        out_dir=args.out_dir,
        nt=args.nt,
        cb_list=args.cb_list,
        info=args.info
    )


if __name__ == '__main__':
    main()
