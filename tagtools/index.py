#!/usr/bin/env python3
"""
TagTools Index Module - Create index for bgzip compressed files

Usage: tagtools index [options]
"""

import os
import sys
import argparse
import numpy as np


def uint16(data):
    """Convert bytes to uint16"""
    return np.frombuffer(data, dtype=np.uint16)[0]


def index_bgzip(file, index_name=None):
    """
    Create index for bgzip compressed file
    
    Args:
        file: Path to bgzip compressed file
        index_name: Output path for index file (default: file.gzi)
    
    Returns:
        str: Path to the created index file
    """
    if not os.path.exists(file):
        raise FileNotFoundError(f'File not found: {file}')
    
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


def main():
    """Command-line entry point for index module"""
    parser = argparse.ArgumentParser(
        prog='tagtools index',
        description='Create index for bgzip compressed file'
    )
    parser.add_argument('--file', required=True, help='Path to bgzip compressed file')
    parser.add_argument('-o', '--index-name', help='Output path for index file (default: file.gzi)')
    
    args = parser.parse_args()
    
    try:
        result = index_bgzip(file=args.file, index_name=args.index_name)
        print(f'Index file created: {result}')
    except Exception as e:
        print(f'Error: {e}', file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
