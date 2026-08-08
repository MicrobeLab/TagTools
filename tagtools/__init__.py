"""
TagTools: A single-cell demultiplexing tool without prior genotype knowledge.

Modules:
    reference: Create haplotype panel from VCF files (reference.py)
    cluster: Cluster prediction from SNP data (cluster.py)
    impute: Imputed haplotype using groups of cells (impute.py)
    demultiplex: Demultiplexing using demuxlet (demultiplex.py)
    index: Create index for bgzip compressed file (index.py)
    split: Split BAM file by tag (split.py)
    call: Call variants and perform genotyping (call.py)
"""

__version__ = "1.0.0"
__author__ = "Guo-wang Lin"
__email__ = "lingw6@alumni.sysu.edu.cn"

from . import reference, cluster, impute, demultiplex, index, split, call
