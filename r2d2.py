import argparse
import pandas as pd
import ConfigParser
import sys
import logging
from scenarios import Scenario
from collections import OrderedDict

CONFIG_FILENAME='r2d2.ini'

if __name__ == "__main__":
    # Arguments may override config settings
    config = ConfigParser.ConfigParser()
    config.read(CONFIG_FILENAME)
    vaf_column_name = config.get('Settings', 'vaf_column_name')
    scenarios_config_filename = config.get('Settings', 'scenarios_config_file')
    logging.basicConfig(level=logging.INFO)

    # Argument parsing
    parser = argparse.ArgumentParser()
    parser.add_argument('-dn', '--dna_normal', type=argparse.FileType('r'), help='Path to DNA normal MAF.')
    parser.add_argument('-dt', '--dna_tumor', type=argparse.FileType('r'), help='Path to DNA tumor MAF.')
    parser.add_argument('-rn', '--rna_normal', type=argparse.FileType('r'), help='Path to RNA normal MAF.')
    parser.add_argument('-rt', '--rna_tumor', type=argparse.FileType('r'), help='Path to RNA tumor MAF.')
    parser.add_argument('-o', '--output', type=argparse.FileType('w'), help='Path to output file.')

    parser.add_argument('--vaf_column_name', type=argparse.FileType('w'),
                        help='Name of column containing variant allele frequency (default=%s).' % vaf_column_name)
    parser.add_argument('--scenarios_config', type=argparse.FileType('w'),
                        help='Path to scenario configuration file (default=%s).' % scenarios_config_filename)

    parser.add_argument('-dnc', '--dna_normal_column', type=str,
                        help='Name of column in DNA normal MAF containing variant allele frequency.'
                             'Overrides --vaf-column-name')
    parser.add_argument('-dtc', '--dna_tumor_column', type=str,
                        help='Name of column in DNA tumor MAF containing variant allele frequency.'
                             'Overrides --vaf-column-name')
    parser.add_argument('-rnc', '--rna_normal_column', type=str,
                        help='Name of column in RNA normal MAF containing variant allele frequency.'
                             'Overrides --vaf-column-name')
    parser.add_argument('-rtc', '--rna_tumor_column', type=str,
                        help='Name of column in RNA tumor MAF containing variant allele frequency.'
                             'Overrides --vaf-column-name')
    args = parser.parse_args()

    # Get input file set and set MAF VAF column names
    maf_types = ['dna_normal', 'dna_tumor', 'rna_normal', 'rna_tumor']
    vaf_columns = {}

    # We cannot modify a list while looping over it, so we use a list comprehension to make a copy
    for maf_type in [maf_type for maf_type in maf_types]:
        if not getattr(args, maf_type):
            maf_types.remove(maf_type)
            continue

        vaf_column = getattr(args, '%s_column' % maf_type)
        if vaf_column:
            vaf_columns[maf_type] = vaf_column
        elif args.vaf_column_name:
            vaf_columns[maf_type] = args.vaf_column_name
        else:
            vaf_columns[maf_type] = vaf_column_name

    # Check that at least one input file was provided.
    if not maf_types:
        logging.error('No input files provided.')
        sys.exit(2)

    # Load all input files into pandas DFs.
    read_csv_args = {'sep': '\t', 'comment': '#', 'skip_blank_lines': True, 'header': 0}
    input_mafs = {}
    merge_columns = ['Hugo_Symbol', 'Chromosome', 'Start_position', 'End_position',
                     'Strand', 'Variant_Classification', 'Variant_Type']
    for maf_type in maf_types:
        maf_arg = getattr(args, maf_type)
        if maf_arg:
            input_mafs[maf_type] = pd.read_csv(maf_arg, **read_csv_args)
            logging.info('Loaded %s data from %s.' % (maf_type, maf_arg.name))

            # Correct column names (add maf_type suffix)
            new_columns = []
            for column in input_mafs[maf_type].columns:
                new_column = column
                # Don't change merge column names
                if column not in merge_columns:
                    new_column += '_' + maf_type

                new_columns.append(new_column)

            input_mafs[maf_type].columns = new_columns

    # Merge DFs on position into one agg DF.
    input_merge = None
    for maf_type in maf_types:
        if input_merge is None:
            input_merge = input_mafs[maf_type]
        else:
            input_merge = input_merge.merge(input_mafs[maf_type], how='outer', on=merge_columns)

    logging.info('Analyzing %s total variants...' % len(input_merge.index))

    # Extract dna_normal/dna_tumor/rna_normal/rna_tumor values from DF.
    output_maf_map = OrderedDict([
        ('Hugo_Symbol', 'Hugo_Symbol'),
        ('Chromosome', 'Chromosome'),
        ('Start_position', 'Start_position'),
        ('End_position', 'End_position'),
        ('Strand', 'Strand'),
        ('Variant_Classification', 'Variant_Classification'),
        ('Variant_Type', 'Variant_Type'),
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

    for i, row in input_merge.iterrows():
        # We only consider single-nucleotide polymorphisms (exclude DNPs, TNPs, etc.)
        if row['Variant_Type'] != 'SNP':
            continue

        quad = {}
        for input_type in input_mafs.keys():
            vaf_column = '%s_%s' % (vaf_columns[input_type], input_type)
            if vaf_column not in row.index:
                raise Exception('Error: VAF column %s not found in %s sample.' %
                                (vaf_columns[input_type], input_type))

            quad[input_type] = row[vaf_column] if vaf_column in row.keys() else 0
            if pd.isnull(quad[input_type]):
                quad[input_type] = 0

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

    logging.info('Wrote %s discovered scenarios to %s.' % (len(output_rows), args.output.name))
