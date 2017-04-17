import argparse
import pandas as pd
from scenarios import Scenario
from collections import OrderedDict

AF_COLUMN = 'AF1'

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-dn', '--dna_normal', type=argparse.FileType('r'), help='Path to DNA normal MAF.')
    parser.add_argument('-dt', '--dna_tumor', type=argparse.FileType('r'), help='Path to DNA tumor MAF.')
    parser.add_argument('-rn', '--rna_normal', type=argparse.FileType('r'), help='Path to RNA normal MAF.')
    parser.add_argument('-rt', '--rna_tumor', type=argparse.FileType('r'), help='Path to RNA tumor MAF.')
    parser.add_argument('-o', '--output', type=argparse.FileType('w'), help='Path to output file.')
    args = parser.parse_args()

    # Load all input files into pandas DFs.
    read_csv_args = {'sep': '\t', 'comment': '#', 'skip_blank_lines': True, 'header': 0}
    input_mafs = {
        'dna_normal': pd.read_csv(args.dna_normal, **read_csv_args),
        'dna_tumor': pd.read_csv(args.dna_tumor, **read_csv_args),
        'rna_normal': pd.read_csv(args.rna_normal, **read_csv_args),
        'rna_tumor': pd.read_csv(args.rna_tumor, **read_csv_args)
    }

    for input_type in input_mafs.keys():
        new_columns = []
        for column in input_mafs[input_type].columns:
            new_column = column
            if column not in ['Start_position', 'End_position']:
                new_column += '_' + input_type

            new_columns.append(new_column)

        input_mafs[input_type].columns = new_columns

    # Merge DFs on position into one agg DF.
    inputs_merge = input_mafs['dna_normal'].merge(
        input_mafs['dna_tumor'],
        how='outer',
        on=['Start_position', 'End_position'],
    )

    inputs_merge = inputs_merge.merge(
        input_mafs['rna_normal'],
        how='outer',
        on=['Start_position', 'End_position'],
    )

    inputs_merge = inputs_merge.merge(
        input_mafs['rna_tumor'],
        how='outer',
        on=['Start_position', 'End_position'],
    )

    # Extract dna_normal/dna_tumor/rna_normal/rna_tumor values from DF.
    output_maf_map = OrderedDict([
        ('Hugo_Symbol', 'Hugo_Symbol_dna_normal'),
        ('Chromosome', 'Chromosome_dna_normal'),
        ('Start_position', 'Start_position'),
        ('End_position', 'End_position'),
        ('Strand', 'Strand_dna_normal'),
        ('Variant_Classification', 'Variant_Classification_dna_normal'),
        ('Variant_Type', 'Variant_Type_dna_normal'),
        ('Reference_Allele', 'Reference_Allele_dna_normal'),
        ('DNA_Normal_Allele1', 'Tumor_Seq_Allele1_dna_normal'),
        ('DNA_Normal_Allele2', 'Tumor_Seq_Allele2_dna_normal'),
        ('DNA_Tumor_Allele1', 'Tumor_Seq_Allele1_dna_tumor'),
        ('DNA_Tumor_Allele2', 'Tumor_Seq_Allele2_dna_tumor'),
        ('RNA_Normal_Allele1', 'Tumor_Seq_Allele1_rna_normal'),
        ('RNA_Normal_Allele2', 'Tumor_Seq_Allele2_rna_normal'),
        ('RNA_Tumor_Allele1', 'Tumor_Seq_Allele1_rna_tumor'),
        ('RNA_Tumor_Allele2', 'Tumor_Seq_Allele2_rna_tumor'),
    ])

    output_rows = OrderedDict()
    for output_column in ['scenario'] + output_maf_map.keys():
        output_rows[output_column] = []

    for i, row in inputs_merge.iterrows():
        quad = {}
        for input_type in input_mafs.keys():
            vaf_column = 'AF1_' + input_type
            quad[input_type] = row[vaf_column] if vaf_column in row.keys() else 0

        try:
            scenario = Scenario(quad)
        except Scenario.NoScenarioException as e:
            continue

        output_rows['scenario'].append(scenario.name)
        for output_column in output_maf_map:
            input_column = output_maf_map[output_column]
            output_rows[output_column].append(row[input_column])

    # Write output.
    pd.DataFrame.from_dict(output_rows).to_csv(args.output, index=False, header=True, sep='\t')
