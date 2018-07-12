import argparse
import pandas as pd
import ConfigParser
import os
import sys
import logging
from collections import OrderedDict
from scenarios import Scenario

CONFIG_FILENAME = 'r2d2.ini'


class MafTypes(object):
    dna_normal = 'dna_normal'
    dna_tumor = 'dna_tumor'
    rna_normal = 'rna_normal'
    rna_tumor = 'rna_tumor'

    @classmethod
    def all(cls):
        return [MafTypes.dna_normal,
                MafTypes.dna_tumor,
                MafTypes.rna_normal,
                MafTypes.rna_tumor]


class AnalysisTypes(object):
    all_inputs = 'all_inputs'
    no_rna_normal = 'no_rna_normal'
    dna_only = 'dna_only'
    normal_only = 'normal_only'
    tumor_only = 'tumor_only'


ANALYSIS_SAMPLES = {
    AnalysisTypes.all_inputs: [MafTypes.dna_normal, MafTypes.dna_tumor, MafTypes.rna_normal, MafTypes.rna_tumor],
    AnalysisTypes.no_rna_normal: [MafTypes.dna_normal, MafTypes.dna_tumor, MafTypes.rna_tumor],
    AnalysisTypes.dna_only: [MafTypes.dna_normal, MafTypes.dna_tumor],
    AnalysisTypes.normal_only: [MafTypes.dna_normal, MafTypes.rna_normal],
    AnalysisTypes.tumor_only: [MafTypes.dna_tumor, MafTypes.rna_tumor]
}


class R2D2(object):
    def __init__(self, analysis_type=None,
                 dna_normal=None, dna_tumor=None, rna_normal=None, rna_tumor=None, output=None, total_output=None,
                 config_path=None,
                 dna_normal_ref_count='',
                 dna_normal_alt_count='',
                 dna_tumor_ref_count='',
                 dna_tumor_alt_count='',
                 rna_normal_ref_count='',
                 rna_normal_alt_count='',
                 rna_tumor_ref_count='',
                 rna_tumor_alt_count='',
                 extra_columns=None,
                 dna_normal_extra_columns=None,
                 dna_tumor_extra_columns=None,
                 rna_normal_extra_columns=None,
                 rna_tumor_extra_columns=None,
                 sample_id=''):

        self.analysis_type = analysis_type
        self.dna_normal = dna_normal
        self.dna_tumor = dna_tumor
        self.rna_normal = rna_normal
        self.rna_tumor = rna_tumor
        self.output = output
        self.total_output = total_output
        self.sample_id = sample_id

        # Arguments may override config settings
        config_file_name = CONFIG_FILENAME
        if config_path is not None:
            self.config_file_path = os.path.join(args.config_path, config_file_name)
            config_parser = ConfigParser.ConfigParser()
            config_parser.read(self.config_file_path)
            logging.info('Loading configuration from {}'.format(self.config_file_path))

            scenarios_config_file_path = config_parser.get('Settings', 'scenarios_config_file')
            scenarios_config_file_path = os.path.join(args.config_path, scenarios_config_file_path)
            logging.info('Loading scenario definitions from {}'.format(scenarios_config_file_path))

            sample_id_header = config_parser.get('Settings', 'sample_id_header')

        logging.info('Analysis type is {}'.format(self.analysis_type))
        samples_for_analysis = ANALYSIS_SAMPLES.get(self.analysis_type)

    maf_types = MafTypes.all()
    ref_count_column = 't_ref_count'
    alt_count_column = 't_alt_count'

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

    maf_type_printable = {
            'dna_normal': 'DNA_Normal', 'dna_tumor': 'DNA_Tumor',
            'rna_normal': 'RNA_Normal', 'rna_tumor': 'RNA_Tumor'
    }

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
    ])

    for maf_type in maf_types:
        allele1_k = 'Allele1_%s' % maf_type_printable[maf_type]
        allele1_v = 'Tumor_Seq_Allele1_%s' % maf_type
        allele2_k = 'Allele2_%s' % maf_type_printable[maf_type]
        allele2_v = 'Tumor_Seq_Allele2_%s' % maf_type
        output_maf_map[allele1_k] = allele1_v
        output_maf_map[allele2_k] = allele2_v

    ref_count_column_prefix = 'Ref_Read_Count'
    alt_count_column_prefix = 'Alt_Read_Count'
    vaf_column_prefix = 'VAF'

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

            scenario = Scenario(vaf_quad, scenarios_config_file_path)
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


def get_analysis_type(dna_normal=None, dna_tumor=None, rna_normal=None, rna_tumor=None):
    num_files = len([f for f in [dna_normal, dna_tumor, rna_normal, rna_tumor] if f is not None])
    if num_files == 4:
        # If all four files are present, simply run the canonical R2D2 setup
        if dna_normal and dna_tumor and rna_normal and rna_tumor:
            return [AnalysisTypes.all_inputs]

    elif 1 < num_files <= 3:
        if dna_normal and dna_tumor and rna_tumor:
            # If all three are provided, we will run the algorithm with all three
            return [AnalysisTypes.no_rna_normal]

        analysis_types_list = []
        # There are 3 files but they are not the canonical 3 dna_normal, dna_tumor, rna_tumor combination, or there are
        # only two files
        if dna_tumor and dna_normal:
            analysis_types_list.append(AnalysisTypes.dna_only)
        if dna_tumor and rna_tumor:
            analysis_types_list.append(AnalysisTypes.tumor_only)
        if dna_normal and rna_normal:
            analysis_types_list.append(AnalysisTypes.normal_only)
        return analysis_types_list

    elif num_files <= 1:
        # No analyses specified for 1 or less files
        return []


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

    analysis_types = get_analysis_type(dna_normal=args.dn, dna_tumor=args.dt, rna_normal=args.rn, rna_tumor=args.rt)
    if len(analysis_types) == 0:
        logging.error('No analysis can be conducted with the {} maf files provided.'.format(len({args.dn,
                                                                                                 args.dt,
                                                                                                 args.rn,
                                                                                                 args.rt})))
        sys.exit(0)

    for analysis_type in analysis_types:
        r2d2 = R2D2(analysis_type=analysis_type,
                    dna_normal=args.dn,
                    dna_tumor=args.dt,
                    rna_normal=args.rn,
                    rna_tumor=args.rt,
                    output=args.o,
                    total_output=args.to,
                    config_path=args.cp,
                    dna_normal_ref_count=args.dnrc,
                    dna_normal_alt_count=args.dnac,
                    dna_tumor_ref_count=args.dtrc,
                    dna_tumor_alt_count=args.dtac,
                    rna_normal_ref_count=args.rnrf,
                    rna_normal_alt_count=args.rnac,
                    rna_tumor_ref_count=args.rtrc,
                    rna_tumor_alt_count=args.rtac,
                    extra_columns=args.xc,
                    dna_normal_extra_columns=args.dnxc,
                    dna_tumor_extra_columns=args.dtxc,
                    rna_normal_extra_columns=args.rnxc,
                    rna_tumor_extra_columns=args.rtxc,
                    sample_id=args.id)