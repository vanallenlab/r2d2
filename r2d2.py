import argparse
import pandas as pd
import ConfigParser
import sys
from scenarios import Scenario
from collections import OrderedDict

CONFIG_FILENAME='r2d2.ini'

if __name__ == "__main__":
    # Arguments may override config settings
    config = ConfigParser.ConfigParser()
    config.read(CONFIG_FILENAME)
    vaf_column_name = config.get('Settings', 'vaf_column_name')
    scenarios_config_filename = config.get('Settings', 'scenarios_config_file')

    parser = argparse.ArgumentParser()
    parser.add_argument('-dn', '--dna_normal', type=argparse.FileType('r'), help='Path to DNA normal MAF.')
    parser.add_argument('-dt', '--dna_tumor', type=argparse.FileType('r'), help='Path to DNA tumor MAF.')
    parser.add_argument('-rn', '--rna_normal', type=argparse.FileType('r'), help='Path to RNA normal MAF.')
    parser.add_argument('-rt', '--rna_tumor', type=argparse.FileType('r'), help='Path to RNA tumor MAF.')
    parser.add_argument('-o', '--output', type=argparse.FileType('w'), help='Path to output file.')
    parser.add_argument('--scenarios_config', type=argparse.FileType('w'),
                        help='Path to scenario configuration file (default=%s).' % scenarios_config_filename)
    parser.add_argument('--vaf_column_name', type=argparse.FileType('w'),
                        help='Column name to read (default=%s).' % vaf_column_name)
    args = parser.parse_args()

    if not args.dna_normal and not args.dna_tumor and not args.rna_normal and not args.rna_tumor:
        print 'Error: No input files provided.'
        sys.exit(2)

    # Load all input files into pandas DFs.
    read_csv_args = {'sep': '\t', 'comment': '#', 'skip_blank_lines': True, 'header': 0}
    input_mafs = {}
    for maf_type in ['dna_normal', 'dna_tumor', 'rna_normal', 'rna_tumor']:
        maf_arg = getattr(args, maf_type)
        if maf_arg:
            input_mafs[maf_type] = pd.read_csv(maf_arg, **read_csv_args)

            # Correct column names (add maf_type suffix)
            new_columns = []
            for column in input_mafs[maf_type].columns:
                new_column = column
                # Don't change start/end position column names (we will merge on these columns)
                if column not in ['Start_position', 'End_position']:
                    new_column += '_' + maf_type

                new_columns.append(new_column)

            input_mafs[maf_type].columns = new_columns

    # Merge DFs on position into one agg DF.
    merge_columns = ['Start_position', 'End_position']
    inputs_merge = input_mafs['dna_normal'].merge(input_mafs['dna_tumor'], how='outer', on=merge_columns)
    inputs_merge = inputs_merge.merge(input_mafs['rna_normal'], how='outer', on=merge_columns)
    inputs_merge = inputs_merge.merge(input_mafs['rna_tumor'], how='outer', on=merge_columns)

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
            scenario = Scenario(quad, scenarios_config_filename)
        except Scenario.NoScenarioException as e:
            continue

        output_rows['scenario'].append(scenario.name)
        for output_column in output_maf_map:
            input_column = output_maf_map[output_column]
            output_rows[output_column].append(row[input_column])

    # Write output.
    pd.DataFrame.from_dict(output_rows).to_csv(args.output, index=False, header=True, sep='\t')
