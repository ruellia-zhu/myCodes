# -*- coding: utf-8 -*-

import os
import csv
import argparse

if __name__ == '__main__':
    print("Generate sample sheet from a dir")
    print("Usage: python3 generate_samplesheet.py -i '~/RNAseq' -o '~/analysis/samplesheet.csv'")

def find_fastq_pairs(input_dir):
    """
    Scan subdirectories of input_dir and find paired FASTQ files.
    Returns a list of tuples: (sample_name, fastq1_abs_path, fastq2_abs_path).
    """
    pairs = []
    for entry in sorted(os.listdir(input_dir)):
        sample_dir = os.path.join(input_dir, entry)
        if not os.path.isdir(sample_dir):
            continue
        files = os.listdir(sample_dir)
        r1 = [f for f in files if f.endswith('_1.fastq') or f.endswith('_1.fastq.gz')]
        r2 = [f for f in files if f.endswith('_2.fastq') or f.endswith('_2.fastq.gz')]
        if not r1 or not r2:
            print(f"Warning: no valid fastq pairs found in {sample_dir}")
            continue
        r1.sort()
        r2.sort()
        for f1, f2 in zip(r1, r2):
            abs_f1 = os.path.abspath(os.path.join(sample_dir, f1))
            abs_f2 = os.path.abspath(os.path.join(sample_dir, f2))
            pairs.append((entry, abs_f1, abs_f2))
    return pairs


def write_samplesheet(pairs, output_csv):
    """
    Write the samplesheet CSV with header and strandedness auto.
    """
    header = ['sample', 'fastq_1', 'fastq_2', 'strandedness']
    with open(output_csv, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(header)
        for sample, f1, f2 in pairs:
            writer.writerow([sample, f1, f2, 'auto'])

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate samplesheet.csv for nf-core RNAseq')
    parser.add_argument('-i', '--input_dir', required=True,
                        help='Path to RNAseq directory containing sample subdirectories')
    parser.add_argument('-o', '--output', default='samplesheet.csv',
                        help='Name of output CSV file')
    args = parser.parse_args()

    pairs = find_fastq_pairs(args.input_dir)
    if not pairs:
        print('No FASTQ pairs found, exiting.')
        exit(1)
    write_samplesheet(pairs, args.output)
    print(f'Samplesheet written to {args.output} with {len(pairs)} entries.')

