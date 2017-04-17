import argparse
import pandas as pd
import ConfigParser
from scenarios import Condition

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-dn', '--dna_normal', type=argparse.FileType('r'), help='Path to DNA normal MAF.')
    parser.add_argument('-dt', '--dna_tumor', type=argparse.FileType('r'), help='Path to DNA tumor MAF.')
    parser.add_argument('-rn', '--rna_normal', type=argparse.FileType('r'), help='Path to RNA normal MAF.')
    parser.add_argument('-rt', '--rna_tumor', type=argparse.FileType('r'), help='Path to RNA tumor MAF.')
    args = parser.parse_args()

    read_csv_args = {'sep': '\t', 'comment': '#', 'skip_blank_lines': True, 'header': 0}
    inputs = {
        'dna_normal': pd.read_csv(args.dna_normal, **read_csv_args),
        'dna_tumor': pd.read_csv(args.dna_tumor, **read_csv_args),
        'rna_normal': pd.read_csv(args.rna_normal, **read_csv_args),
        'rna_tumor': pd.read_csv(args.rna_tumor, **read_csv_args)
    }

    inputs_merge = inputs['dna_normal'].merge(
        inputs['dna_tumor'],
        how='outer',
        on=['Start_position', 'End_position'],
        suffixes=['_dn', '_dt']
    )

    inputs_merge = inputs_merge.merge(
        inputs['rna_normal'],
        how='outer',
        on=['Start_position', 'End_position'],
    )

    inputs_merge = inputs_merge.merge(
        inputs['rna_tumor'],
        how='outer',
        on=['Start_position', 'End_position'],
        suffixes=['_rn', '_rt']
    )

    print inputs_merge

    # Load all input files into pandas DFs.
    # Merge DFs on position into one agg DF.
    # Have factory function that creates a scenario class (RNAed/LOH/etc.) based on a single row from above DF.
    # Write output.
