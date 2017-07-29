import argparse
import pandas as pd
import ConfigParser
import os
import sys
import logging
from collections import OrderedDict
from scenarios import Scenario

CONFIG_FILENAME = 'r2d2.ini'

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    # Argument parsing
    parser = argparse.ArgumentParser()
    parser.add_argument('-dn', '--dna_normal', type=argparse.FileType('r'), help='Path to DNA normal MAF.')
    parser.add_argument('-dt', '--dna_tumor', type=argparse.FileType('r'), help='Path to DNA tumor MAF.')
    parser.add_argument('-rn', '--rna_normal', type=argparse.FileType('r'), help='Path to RNA normal MAF.')
    parser.add_argument('-rt', '--rna_tumor', type=argparse.FileType('r'), help='Path to RNA tumor MAF.')
    parser.add_argument('-o', '--output', type=argparse.FileType('w'), help='Path to output file.')
    parser.add_argument('-to', '--total_output', type=argparse.FileType('w'), help='Path to output file.')
    parser.add_argument('-cp', '--config_path', type=str, help='Path to directory containing config files.')

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

    parser.add_argument('-id', '--sample_id', type=str,
                        help='Sample ID string to include in each row of output data.')

    args = parser.parse_args()

    # Arguments may override config settings
    config_file_path = CONFIG_FILENAME
    if args.config_path:
        config_file_path = os.path.join(args.config_path, config_file_path)

    config = ConfigParser.ConfigParser()
    config.read(config_file_path)
    scenarios_config_filename = config.get('Settings', 'scenarios_config_file')
    sample_id_header = config.get('Settings', 'sample_id_header')

    maf_types = ['dna_normal', 'dna_tumor', 'rna_normal', 'rna_tumor']
    ref_count_column = 't_ref_count'
    alt_count_column = 't_alt_count'

    # Remove maf_types we weren't given from the maf_types list.
    # We cannot modify a list while looping over it, so we use a list comprehension to make a copy.
    for maf_type in [maf_type for maf_type in maf_types]:
        # (E.g., --dna_normal provides the maf type 'dna_normal' file.
        if not getattr(args, maf_type, None):
            maf_types.remove(maf_type)
            continue

    # Check that at least one input file was provided.
    if not maf_types:
        logging.error('No input files provided.')
        sys.exit(2)

    # Load all input files into pandas DFs.
    read_csv_args = {'sep': '\t', 'comment': '#', 'skip_blank_lines': True, 'header': 0}
    input_mafs = {}
    merge_columns = ['Hugo_Symbol', 'Chromosome', 'Start_position', 'End_position',
                     'Strand', 'Variant_Classification', 'Variant_Type']

    logging.info('Loading alteration data...')
    for maf_type in maf_types:
        maf_arg = getattr(args, maf_type, None)
        if maf_arg:
            input_mafs[maf_type] = pd.read_csv(maf_arg, low_memory=False, **read_csv_args)
            for merge_column in merge_columns:
                if merge_column not in input_mafs[maf_type].columns:
                    logging.error('Merge column %s not found in %s sample (%s).' %
                                  (merge_column, maf_type, maf_arg.name))
                    sys.exit(2)

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
    ])

    ref_count_column_prefix = 'Ref_Read_Count'
    alt_count_column_prefix = 'Alt_Read_Count'
    vaf_column_prefix = 'VAF'

    maf_type_printable = {
            'dna_normal': 'DNA_Normal', 'dna_tumor': 'DNA_Tumor',
            'rna_normal': 'RNA_Normal', 'rna_tumor': 'RNA_Tumor'
    }

    # Add column mappings for all extra columns seen in all 4 files
    if args.extra_columns:
        for xc in args.extra_columns.split():
            for column_suffix in maf_type_printable:
                this_column_name = '%s_%s' % (xc, maf_type_printable[column_suffix])
                output_maf_map[this_column_name] = '%s_%s' % (xc, column_suffix)

    # Add column mappings for extra columns seen in only one file
    for column_suffix in maf_type_printable:
        this_suffix_columns = getattr(args, '%s_extra_columns' % column_suffix, None)
        if this_suffix_columns:
            for column in this_suffix_columns.split():
                this_column_name = '%s_%s' % (column, maf_type_printable[column_suffix])
                output_maf_map[this_column_name] = '%s_%s' % (column, column_suffix)

    # output_rows will associate the output column names with arrays that will contain that column's values per-row:
    output_rows = OrderedDict()
    expected_rows = OrderedDict()
    for output_column in ['scenario'] + output_maf_map.keys():
        output_rows[output_column] = []
        expected_rows[output_column] = []

    # Create read count arrays:
    for maf_type in maf_types:
        output_rows['%s_%s' % (ref_count_column_prefix, maf_type_printable[maf_type])] = []
        output_rows['%s_%s' % (alt_count_column_prefix, maf_type_printable[maf_type])] = []
        output_rows['%s_%s' % (vaf_column_prefix, maf_type_printable[maf_type])] = []
        expected_rows['%s_%s' % (ref_count_column_prefix, maf_type_printable[maf_type])] = []
        expected_rows['%s_%s' % (alt_count_column_prefix, maf_type_printable[maf_type])] = []
        expected_rows['%s_%s' % (vaf_column_prefix, maf_type_printable[maf_type])] = []

    for i, row in input_merge.iterrows():
        # We exclude DNPs, TNPs, etc.
        if row['Variant_Type'] not in ['SNP', 'INS', 'DEL']:
            continue

        # Quad is a dictionary mapping each maf type to the corresponding VAFs on this row.
        counts = {}
        for input_type in input_mafs.keys():
            ref_column = '%s_%s' % (ref_count_column, input_type)
            if ref_column not in row.index:
                raise Exception('Error: Reference read count column %s not found in %s sample.' %
                                (ref_count_column, input_type))

            alt_column = '%s_%s' % (alt_count_column, input_type)
            if alt_column not in row.index:
                raise Exception('Error: Alternate read count column %s not found in %s sample.' %
                                (alt_count_column, input_type))

            # Cast ref and alt counts to float; if unable to do so, set as nan.
            try:
                ref_count = float(row[ref_column]) if not pd.isnull(row[ref_column]) else 0
            except ValueError:
                ref_count = float('nan')

            try:
                alt_count = float(row[alt_column]) if not pd.isnull(row[alt_column]) else 0
            except ValueError:
                alt_count = float('nan')

            if ref_count == 0 and alt_count == 0:
                vaf = 0
            elif ref_count == float('nan') or alt_count == float('nan'):
                vaf = float('nan')
            else:
                vaf = alt_count / (ref_count + alt_count)

            counts[input_type] = {'ref_count': ref_count, 'alt_count': alt_count, 'vaf': vaf}

        # row_dest either points to the output_rows or expected_rows OrderedDicts
        row_dest = None
        scenario_name = None
        try:
            # Construct VAF quad out of counts (pull vaf value out)
            vaf_quad = {}
            for input_type in input_mafs.keys():
                vaf_quad[input_type] = counts[input_type]['vaf']

            scenario = Scenario(vaf_quad, scenarios_config_filename)
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
            ref_output_column = '%s_%s' % (ref_count_column_prefix, maf_type_printable[maf_type])
            alt_output_column = '%s_%s' % (alt_count_column_prefix, maf_type_printable[maf_type])
            vaf_output_column = '%s_%s' % (vaf_column_prefix, maf_type_printable[maf_type])
            row_dest[ref_output_column].append(counts[maf_type]['ref_count'])
            row_dest[alt_output_column].append(counts[maf_type]['alt_count'])
            row_dest[vaf_output_column].append(counts[maf_type]['vaf'])

    # Write output.
    output_df = pd.DataFrame.from_dict(output_rows)

    # Prepend sample_id to each row if present
    sample_id = getattr(args, 'sample_id', None)
    if sample_id:
        output_df.insert(0, sample_id_header, sample_id)

    output_df.to_csv(args.output, index=False, header=True, sep='\t')
    logging.info('Wrote %s discovered scenarios to %s.' % (len(output_df.index), args.output.name))

    if args.total_output:
        expected_df = pd.DataFrame.from_dict(expected_rows)
        if sample_id:
            expected_df.insert(0, sample_id_header, sample_id)

        total_df = pd.concat([output_df, expected_df], axis=0)
        total_df.to_csv(args.total_output, index=False, header=True, sep='\t')
        logging.info('Wrote %s total variants to %s.' % (len(total_df.index), args.total_output.name))
