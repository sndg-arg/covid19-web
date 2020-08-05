#!/usr/bin/env python3
if __name__ == '__main__':
    import argparse
    from glob import glob

    parser = argparse.ArgumentParser(description='list')
    parser.add_argument('input_dir', help="where samples processed with process_batch.py are")

    args = parser.parse_args()

    for diff_file in glob(f'{args.input_dir}/*/nucdiff_denovo_consensus/*.vcf'):
        results = []
        with open(diff_file) as h:
            for line in h:
                if not line.startswith("#"):
                    results.append(line)
        print(f'processing {diff_file}: {len(results)} differences found' )
        for line in results:
            print(f'{line}')
