import argparse
import pandas as pd
from scenarios import Scenario
from collections import OrderedDict

AF_COLUMN = 'AF1'
OUTPUT_FILE = 'output.tsv'

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-dn', '--dna_normal', type=argparse.FileType('r'), help='Path to DNA normal MAF.')
    parser.add_argument('-dt', '--dna_tumor', type=argparse.FileType('r'), help='Path to DNA tumor MAF.')
    parser.add_argument('-rn', '--rna_normal', type=argparse.FileType('r'), help='Path to RNA normal MAF.')
    parser.add_argument('-rt', '--rna_tumor', type=argparse.FileType('r'), help='Path to RNA tumor MAF.')
    args = parser.parse_args()

    # Load all input files into pandas DFs.
    read_csv_args = {'sep': '\t', 'comment': '#', 'skip_blank_lines': True, 'header': 0}
    inputs = {
        'dna_normal': pd.read_csv(args.dna_normal, **read_csv_args),
        'dna_tumor': pd.read_csv(args.dna_tumor, **read_csv_args),
        'rna_normal': pd.read_csv(args.rna_normal, **read_csv_args),
        'rna_tumor': pd.read_csv(args.rna_tumor, **read_csv_args)
    }

    # Merge DFs on position into one agg DF.
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

    # Extract dna_normal/dna_tumor/rna_normal/rna_tumor values from DF.
    output_maf_map = OrderedDict([
        ('Hugo_Symbol', 'Hugo_Symbol_dn'),
        ('Chromosome', 'Chromosome_dn'),
        ('Start_position', 'Start_position'),
        ('End_position', 'End_position'),
        ('Strand', 'Strand_dn'),
        ('Variant_Classification', 'Variant_Classification_dn'),
        ('Variant_Type', 'Variant_Type_dn'),
        ('Reference_Allele', 'Reference_Allele_dn'),
        ('DNA_Normal_Allele1', 'Tumor_Seq_Allele1_dn'),
        ('DNA_Normal_Allele2', 'Tumor_Seq_Allele2_dn'),
        ('DNA_Tumor_Allele1', 'Tumor_Seq_Allele1_dt'),
        ('DNA_Tumor_Allele2', 'Tumor_Seq_Allele2_dt'),
        ('RNA_Normal_Allele1', 'Tumor_Seq_Allele1_rn'),
        ('RNA_Normal_Allele2', 'Tumor_Seq_Allele2_rn'),
        ('RNA_Tumor_Allele1', 'Tumor_Seq_Allele1_rt'),
        ('RNA_Tumor_Allele2', 'Tumor_Seq_Allele2_rt'),
    ])

    output_rows = OrderedDict()
    for output_column in ['scenario'] + output_maf_map.keys():
        output_rows[output_column] = []

    for i, row in inputs_merge.iterrows():
        dna_normal = row['AF1_dn'] if 'AF1_dn' in row.keys() else 0
        dna_tumor = row['AF1_dt'] if 'AF1_dt' in row.keys() else 0
        rna_normal = row['AF1_rn'] if 'AF1_rn' in row.keys() else 0
        rna_tumor = row['AF1_rt'] if 'AF1_rt' in row.keys() else 0

        quad = {
            'dna_normal': dna_normal,
            'dna_tumor': dna_tumor,
            'rna_normal': rna_normal,
            'rna_tumor': rna_tumor,
        }

        try:
            scenario = Scenario(quad)
        except Scenario.NoScenarioException as e:
            continue

        output_rows['scenario'].append(scenario.name)
        for output_column in output_maf_map:
            input_column = output_maf_map[output_column]
            output_rows[output_column].append(row[input_column])

    # Write output.
    pd.DataFrame.from_dict(output_rows).to_csv(OUTPUT_FILE, index=False, header=True, sep='\t')
