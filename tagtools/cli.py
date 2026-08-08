#!/usr/bin/env python3
"""
TagTools CLI - Main entry point for the TagTools package.

Usage: tagtools <command> [options]

Commands:
    reference      Create haplotype panel from VCF files (Step 1)
    cluster        Cluster prediction from SNP data (Step 2)
    impute         Imputed haplotype using groups of cells (Step 3)
    demultiplex    Demultiplexing using demuxlet (Step 4)
    index          Create index for bgzip compressed file
    split          Split BAM file by tag
    call           Call variants and perform genotyping
"""

import argparse
import sys
from . import reference, cluster, impute, demultiplex, index, split, call


def main():
    # Create main parser
    parser = argparse.ArgumentParser(
        prog='tagtools',
        description='TagTools: A single-cell demultiplexing tool without prior genotype knowledge.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    tagtools reference --vcf-list vcf_files.txt --out-name panel --nt 48
    tagtools cluster --cb-list barcodes.txt --bam input.bam --snp-site sites.gz --sample-num 3 --ref genome.fa --out-name cluster
    tagtools impute --bam input.bam --cb-group groups.txt --snp-site sites.gz --hap panel.bcf --chunk regions.txt --ref genome.fa --out-name impute
    tagtools demultiplex --vcf imputed.bcf --bam input.bam --cb-list barcodes.txt --out-dir demux --nt 48
    tagtools index --file data.gz
    tagtools split --bam sample.bam --tag-list tags.txt --out-dir output/
    tagtools call --bam sample.bam --ref genome.fa --site sites.txt --out-name results
        """
    )

    # Create subparsers for each command
    subparsers = parser.add_subparsers(
        title='Commands',
        dest='command',
        help='Available commands',
        required=True
    )

    # === Command 1: reference (make_ref.py) ===
    ref_parser = subparsers.add_parser(
        'reference',
        help='Create haplotype panel from VCF files (Step 1)',
        description='Create a haplotype panel from VCF files for downstream analysis.'
    )
    ref_parser.add_argument('--vcf-list', type=str, required=True,
                          help='Path to a text file containing a list of bgzip-compressed haplotype VCF files.')
    ref_parser.add_argument('--out-name', type=str, required=True,
                          help='Prefix for the output file names.')
    ref_parser.add_argument('--nt', type=int, required=True,
                          help='Number of threads to use for processing.')
    ref_parser.add_argument('--add-chr', action='store_true',
                          help='Whether to add the "chr" prefix to the haplotype panel.')
    ref_parser.add_argument('--rm-chr', action='store_true',
                          help='Whether to remove the "chr" prefix from the haplotype panel.')
    ref_parser.set_defaults(func=reference.main)

    # === Command 2: cluster (cluster_predict.py) ===
    cluster_parser = subparsers.add_parser(
        'cluster',
        help='Cluster prediction from SNP data (Step 2)',
        description='Perform genotype similarity detection and clustering for single cells.'
    )
    cluster_parser.add_argument('--bam', required=True,
                              help='BAM file, should be sorted and indexed by samtools.')
    cluster_parser.add_argument('--out-name', required=True,
                              help='Output file name prefix.')
    cluster_parser.add_argument('--cb-list', required=True,
                              help='Cell barcode list file.')
    cluster_parser.add_argument('--snp-site', required=True,
                              help='Common SNPs site file, indexed with tabix.')
    cluster_parser.add_argument('--ref', required=True,
                              help='Reference genome file.')
    cluster_parser.add_argument('--sample-num', type=int, required=True,
                              help='Number of samples multiplexed.')
    cluster_parser.add_argument('--nt', '--threads', type=int, default=16,
                              help='Number of threads to use, default: 16.')
    cluster_parser.add_argument('--plot-log', action='store_true',
                              help='Plot the log file, default: false.')
    cluster_parser.set_defaults(func=cluster.main)

    # === Command 3: impute (imputed.py) ===
    impute_parser = subparsers.add_parser(
        'impute',
        help='Imputed haplotype using groups of cells (Step 3)',
        description='Perform low-depth sequencing imputation using sequencing data from different droplet groups.'
    )
    impute_parser.add_argument('--bam', required=True,
                             help='BAM file, should be sorted and indexed by samtools.')
    impute_parser.add_argument('--out-name', required=True,
                             help='Output file name prefix.')
    impute_parser.add_argument('--ref', required=True,
                             help='Reference genome file.')
    impute_parser.add_argument('--snp-site', required=True,
                             help='SNPs site file of the haplotype panel, indexed with tabix.')
    impute_parser.add_argument('--hap', required=True,
                             help='Haplotype panel file.')
    impute_parser.add_argument('--chunk', required=True,
                             help='Chunked region file.')
    impute_parser.add_argument('--cb-group', required=True,
                             help='Cell barcode group file.')
    impute_parser.add_argument('--nt', '--threads', type=int, default=16,
                             help='Number of threads to use, default: 16.')
    impute_parser.set_defaults(func=impute.main)

    # === Command 4: demultiplex (demultiplexing.py) ===
    demux_parser = subparsers.add_parser(
        'demultiplex',
        help='Demultiplexing using demuxlet (Step 4)',
        description='Perform sample demultiplexing using demuxlet with imputed genotype file.'
    )
    demux_parser.add_argument('--vcf', required=True,
                            help='Imputed VCF/BCF file, indexed with bcftools.')
    demux_parser.add_argument('--bam', required=True,
                            help='BAM file, sorted and indexed with samtools.')
    demux_parser.add_argument('--cb-list', required=True,
                            help='Cell barcode list file.')
    demux_parser.add_argument('--out-dir', required=True,
                            help='Output directory name.')
    demux_parser.add_argument('--nt', '--threads', type=int, default=16,
                            help='Number of threads to use, default: 16.')
    demux_parser.add_argument('--info', type=float, default=0.4,
                            help='INFO score to filter imputed genotype, default: 0.4.')
    demux_parser.set_defaults(func=demultiplex.main)

    # === Command 5: index (index.py) ===
    index_parser = subparsers.add_parser(
        'index',
        help='Create index for bgzip compressed file',
        description='Create index for bgzip compressed file to enable random access.'
    )
    index_parser.add_argument('--file', required=True,
                            help='Path to bgzip compressed file')
    index_parser.add_argument('-o', '--index-name',
                            help='Output path for index file (default: file.gzi)')
    index_parser.set_defaults(func=index.main)

    # === Command 6: split (split.py) ===
    split_parser = subparsers.add_parser(
        'split',
        help='Split BAM file by tag',
        description='Split a BAM file into separate files based on a specified tag value.'
    )
    split_parser.add_argument('--bam', required=True,
                            help='Input BAM file path')
    split_parser.add_argument('--nt', type=int, default=4,
                            help='Number of threads (default: 4)')
    split_parser.add_argument('--out-dir', required=True,
                            help='Output directory')
    split_parser.add_argument('--tag-name', default='CB',
                            help='Tag name (default: CB)')
    split_parser.add_argument('--tag-list', required=True,
                            help='Tag list file')
    split_parser.add_argument('--region-file',
                            help='Region file (optional)')
    split_parser.set_defaults(func=split.main)

    # === Command 7: call (call.py) ===
    call_parser = subparsers.add_parser(
        'call',
        help='Call variants and perform genotyping',
        description='Call variants from BAM files and compute genotype likelihoods (PL values).'
    )
    call_parser.add_argument('--bam', required=True,
                           help='Input BAM file path')
    call_parser.add_argument('--nt', type=int, default=4,
                           help='Number of threads (default: 4)')
    call_parser.add_argument('--out-name', required=True,
                           help='Output file name')
    call_parser.add_argument('--ref', required=True,
                           help='Reference genome file')
    call_parser.add_argument('--tag-name', default='CB',
                           help='Tag name (default: CB)')
    call_parser.add_argument('--tag-list', required=True,
                           help='Tag list file')
    call_parser.add_argument('--site', required=True,
                           help='Site file')
    call_parser.add_argument('--chunk-size', type=int, default=20,
                           help='Chunk size (default: 20)')
    call_parser.set_defaults(func=call.main)

    # Parse arguments
    args = parser.parse_args()

    # Execute the selected command
    try:
        args.func(args)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
