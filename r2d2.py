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
    scenarios_config_filename = config.get('Settings', 'scenarios_config_file')
    logging.basicConfig(level=logging.INFO)

    # Argument parsing
    parser = argparse.ArgumentParser()
    parser.add_argument('-dn', '--dna_normal', type=argparse.FileType('r'), help='Path to DNA normal MAF.')
    parser.add_argument('-dt', '--dna_tumor', type=argparse.FileType('r'), help='Path to DNA tumor MAF.')
    parser.add_argument('-rn', '--rna_normal', type=argparse.FileType('r'), help='Path to RNA normal MAF.')
    parser.add_argument('-rt', '--rna_tumor', type=argparse.FileType('r'), help='Path to RNA tumor MAF.')
    parser.add_argument('-o', '--output', type=argparse.FileType('w'), help='Path to output file.')
    parser.add_argument('-to', '--total_output', type=argparse.FileType('w'), help='Path to output file.')

    parser.add_argument('-dnv', '--dna_normal_vaf', type=str,
                        help='Name of column in DNA normal MAF containing variant allele frequency.'
                             'Overrides --vaf-column-name')
    parser.add_argument('-dtv', '--dna_tumor_vaf', type=str,
                        help='Name of column in DNA tumor MAF containing variant allele frequency.'
                             'Overrides --vaf-column-name')
    parser.add_argument('-rnv', '--rna_normal_vaf', type=str,
                        help='Name of column in RNA normal MAF containing variant allele frequency.'
                             'Overrides --vaf-column-name')
    parser.add_argument('-rtv', '--rna_tumor_vaf', type=str,
                        help='Name of column in RNA tumor MAF containing variant allele frequency.'
                             'Overrides --vaf-column-name')

    parser.add_argument('-dnrc', '--dna_normal_ref_count', type=str,
                        help='Column name containing count of reference reads in DNA normal file.')
    parser.add_argument('-dnac', '--dna_normal_alt_count', type=str,
                        help='Column name containing count of alternate reads in DNA normal file.')
    parser.add_argument('-dtrc', '--dna_tumor_ref_count', type=str,
                        help='Column name containing count of reference reads in DNA tumor file.')
    parser.add_argument('-dtac', '--dna_tumor_alt_count', type=str,
                        help='Column name containing count of alternate reads in DNA tumor file.')
    parser.add_argument('-rnrc', '--rna_normal_ref_count', type=str,
                        help='Column name containing count of reference reads in RNA normal file.')
    parser.add_argument('-rnac', '--rna_normal_alt_count', type=str,
                        help='Column name containing count of alternate reads in RNA normal file.')
    parser.add_argument('-rtrc', '--rna_tumor_ref_count', type=str,
                        help='Column name containing count of reference reads in RNA tumor file.')
    parser.add_argument('-rtac', '--rna_tumor_alt_count', type=str,
                        help='Column name containing count of alternate reads in RNA tumor file.')

    parser.add_argument('-xc', '--extra_columns', type=str,
                        help='Space-separated list of extra columns to include from each input file.')
    parser.add_argument('-dnxc', '--dna_normal_extra_columns', type=str,
                        help='Space-separated list of extra columns to include from the DNA normal file.')
    parser.add_argument('-dtxc', '--dna_tumors_extra_columns', type=str,
                        help='Space-separated list of extra columns to include from the DNA tumor file.')
    parser.add_argument('-rnxc', '--rna_normal_extra_columns', type=str,
                        help='Space-separated list of extra columns to include from the RNA normal file.')
    parser.add_argument('-rtxc', '--rna_tumor_extra_columns', type=str,
                        help='Space-separated list of extra columns to include from the RNA tumor file.')

    args = parser.parse_args()

    # Get input file set and set MAF VAF column names
    maf_types = ['dna_normal', 'dna_tumor', 'rna_normal', 'rna_tumor']

	if args.ref_column:
		ref_count_column_name = args.ref_column
	if args.alt_column:
		alt_count_column_name = args.alt_column

    # We cannot modify a list while looping over it, so we use a list comprehension to make a copy
    for maf_type in [maf_type for maf_type in maf_types]:
        if not getattr(args, maf_type, None):
            maf_types.remove(maf_type)
            continue

        vaf_column = getattr(args, '%s_vaf' % maf_type, None)
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
        maf_arg = getattr(args, maf_type, None)
        if maf_arg:
            input_mafs[maf_type] = pd.read_csv(maf_arg, low_memory=False, **read_csv_args)
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
        ('Allele1_DNA_Normal', 'Tumor_Seq_Allele1_dna_normal'),
        ('Allele2_DNA_Normal', 'Tumor_Seq_Allele2_dna_normal'),
        ('Allele1_DNA_Tumor', 'Tumor_Seq_Allele1_dna_tumor'),
        ('Allele2_DNA_Tumor', 'Tumor_Seq_Allele2_dna_tumor'),
        ('Allele1_RNA_Normal', 'Tumor_Seq_Allele1_rna_normal'),
        ('Allele2_RNA_Normal', 'Tumor_Seq_Allele2_rna_normal'),
        ('Allele1_RNA_Tumor', 'Tumor_Seq_Allele1_rna_tumor'),
        ('Allele2_RNA_Tumor', 'Tumor_Seq_Allele2_rna_tumor'),
        ('AF_DNA_Normal', '%s_%s' % (vaf_columns['dna_normal'], 'dna_normal')),
        ('AF_DNA_Tumor', '%s_%s' % (vaf_columns['dna_tumor'], 'dna_tumor')),
        ('AF_RNA_Normal', '%s_%s' % (vaf_columns['rna_normal'], 'rna_normal')),
        ('AF_RNA_Tumor', '%s_%s' % (vaf_columns['rna_tumor'], 'rna_tumor')),
    ])

    read_count_column_prefix = 'Read_Count'

    maf_type_column_suffixes = {
            'dna_normal': 'DNA_Normal', 'dna_tumor': 'DNA_Tumor',
            'rna_normal': 'RNA_Normal', 'rna_tumor': 'RNA_Tumor'
    }

    # Add column mappings for all extra columns seen in all 4 files
    if args.extra_columns:
        for xc in args.extra_columns.split():
            for column_suffix in maf_type_column_suffixes:
                this_column_name = '%s_%s' % (xc, maf_type_column_suffixes[column_suffix])
                output_maf_map[this_column_name] = '%s_%s' % (xc, column_suffix)

    # Add column mappings for extra columns seen in only one file
    for column_suffix in maf_type_column_suffixes:
        this_suffix_columns = getattr(args, '%s_extra_columns' % column_suffix, None)
        if this_suffix_columns:
            for column in this_suffix_columns.split():
                this_column_name = '%s_%s' % (column, maf_type_column_suffixes[column_suffix])
                output_maf_map[this_column_name] = '%s_%s' % (column, column_suffix)

    # output_rows will associate the output column names with arrays that will contain that column's values per-row:
    output_rows = OrderedDict()
    expected_rows = OrderedDict()
    for output_column in ['scenario'] + output_maf_map.keys():
        output_rows[output_column] = []
        expected_rows[output_column] = []

    # Create read count arrays:
    for maf_type in maf_types:
        output_rows['%s_%s' % (read_count_column_prefix, maf_type)] = []
        expected_rows['%s_%s' % (read_count_column_prefix, maf_type)] = []

    for i, row in input_merge.iterrows():
        # We only consider single-nucleotide polymorphisms (exclude DNPs, TNPs, etc.)
        if row['Variant_Type'] != 'SNP':
            continue

        # Quad is a dictionary mapping each maf type to the corresponding VAFs on this row.
        quad = {}
        for input_type in input_mafs.keys():
            vaf_column = '%s_%s' % (vaf_columns[input_type], input_type)
            if vaf_column not in row.index:
                raise Exception('Error: VAF column %s not found in %s sample.' %
                                (vaf_columns[input_type], input_type))

            quad[input_type] = row[vaf_column] if vaf_column in row.keys() else 0
            if pd.isnull(quad[input_type]):
                quad[input_type] = 0

        # row_dest either points to the output_rows or expected_rows OrderedDicts
        row_dest = None
        scenario_name = None
        try:
            scenario = Scenario(quad, scenarios_config_filename)
            row_dest = output_rows
            scenario_name = scenario.name
        except Scenario.NoScenarioException as e:
            if not args.total_output:
                continue
            else:
                row_dest = expected_rows

        row_dest['scenario'].append(scenario_name)
        for output_column in output_maf_map:
            input_column = output_maf_map[output_column]
            row_dest[output_column].append(row[input_column])

        for maf_type in maf_types:
            ref_count_column = getattr(args, '%s_ref_count' % maf_type, None)
            alt_count_column = getattr(args, '%s_alt_count' % maf_type, None)

            if not ref_count_column:
                ref_count_column = args.ref_count

            if not alt_count_column:
                alt_count_column = args.alt_count

            if ref_count_column and alt_count_column:
                # convert to type column namespace
                ref_count_column = '%s_%s' % (ref_count_column, maf_type)
                alt_count_column = '%s_%s' % (alt_count_column, maf_type)

                if ref_count_column not in row.index:
                    raise Exception('Error: Ref count column %s not found in %s sample.' %
                                    (ref_count_column, input_type))

                if alt_count_column not in row.index:
                    raise Exception('Error: Alt count column %s not found in %s sample.' %
                                    (alt_count_column, input_type))

                ref_count = row[ref_count_column]
                alt_count = row[alt_count_column]

                if pd.isnull(ref_count):
                    ref_count = 0
                if pd.isnull(alt_count):
                    alt_count = 0

                this_column_prefix = '%s_%s' % (read_count_column_prefix, maf_type)
                row_dest[this_column_prefix].append(ref_count + alt_count)


    # Write output.
    output_df = pd.DataFrame.from_dict(output_rows)
    output_df.to_csv(args.output, index=False, header=True, sep='\t')
    logging.info('Wrote %s discovered scenarios to %s.' % (len(output_df.index), args.output.name))

    if args.total_output:
        total_df = pd.concat([output_df, pd.DataFrame.from_dict(expected_rows)], axis=0)
        total_df.to_csv(args.total_output, index=False, header=True, sep='\t')
        logging.info('Wrote %s total variants to %s.' % (len(total_df.index), args.total_output.name))
