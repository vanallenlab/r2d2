import ConfigParser
import argparse
import logging
import os
import sys
from collections import OrderedDict

import pandas as pd

from helpers.maf_types import AnalysisTypes, ANALYSIS_SAMPLES
from helpers.scenarios import ScenarioCalculator

CONFIG_FILENAME = 'r2d2.ini'


class R2D2ParsingException(Exception):
    def __init__(self, error_args):
        Exception.__init__(self, error_args)


MAF_TYPE_PRINTABLE = {
    'dna_normal': 'DNA_Normal', 'dna_tumor': 'DNA_Tumor',
    'rna_normal': 'RNA_Normal', 'rna_tumor': 'RNA_Tumor'
}

# Extract dna_normal/dna_tumor/rna_normal/rna_tumor values from DF.
OUTPUT_MAF_MAP = OrderedDict([
    ('Hugo_Symbol', 'Hugo_Symbol'),
    ('Chromosome', 'Chromosome'),
    ('Start_position', 'Start_position'),
    ('End_position', 'End_position'),
    ('Strand', 'Strand'),
    ('Variant_Classification', 'Variant_Classification'),
    ('Variant_Type', 'Variant_Type'),
    ('Reference_Allele', 'Reference_Allele_dna_normal'),
])


class R2D2(object):
    def __init__(self, analysis_type=None,
                 dna_normal=None, dna_tumor=None, rna_normal=None, rna_tumor=None, output_location=None,
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
                 sample_id='',
                 load_data=True):

        self.analysis_type = analysis_type
        self.dna_normal = dna_normal
        self.dna_tumor = dna_tumor
        self.rna_normal = rna_normal
        self.rna_tumor = rna_tumor
        self.output_location = output_location
        self.sample_id = sample_id
        self.extra_columns = extra_columns

        self.dna_normal_extra_columns = dna_normal_extra_columns,
        self.dna_tumor_extra_columns = dna_tumor_extra_columns,
        self.rna_normal_extra_columns = rna_normal_extra_columns,
        self.rna_tumor_extra_columns = rna_tumor_extra_columns,

        self.ref_count_column = 't_ref_count'
        self.alt_count_column = 't_alt_count'

        self.read_csv_args = {'sep': '\t', 'comment': '#', 'skip_blank_lines': True, 'header': 0}
        self.merge_columns = ['Hugo_Symbol', 'Chromosome', 'Start_position', 'End_position', 'Strand',
                              'Variant_Classification', 'Variant_Type']
        # Arguments may override config settings
        config_file_name = CONFIG_FILENAME
        if config_path is None:
            config_path = os.path.dirname(os.path.abspath(__file__))

        self.config_file_path = os.path.join(config_path, config_file_name)
        config_parser = ConfigParser.ConfigParser()
        config_parser.read(self.config_file_path)
        logging.info('Loading configuration from {}'.format(self.config_file_path))

        self.scenarios_config_file_path = os.path.join(config_path,
                                                       config_parser.get('Settings', 'scenarios_config_file'))
        logging.info('Loading scenario definitions from {}'.format(self.scenarios_config_file_path))
        self.sample_id_header = config_parser.get('Settings', 'sample_id_header')

        logging.info('Analysis type is {}'.format(self.analysis_type))
        self.maf_types_for_analysis = ANALYSIS_SAMPLES.get(self.analysis_type, [])
        self.maf_sample_paths = {maf_type: self.__getattribute__(maf_type) for maf_type in self.maf_types_for_analysis}
        self.input_merge = None
        if load_data:
            self.load_data()

    def load_data(self):
        """ Load all input files into pandas DFs. """
        input_mafs = {}
        logging.info('Loading alteration data...')
        for maf_type, maf_sample_path in self.maf_sample_paths.items():
            input_mafs[maf_type] = pd.read_csv(maf_sample_path, low_memory=False, **self.read_csv_args)
            for merge_column in self.merge_columns:
                if merge_column not in input_mafs[maf_type].columns:
                    error_message = 'Merge column {} not found in {} sample ({}).'.format(merge_column,
                                                                                          maf_type,
                                                                                          maf_sample_path)
                    logging.error(error_message)
                    raise R2D2ParsingException(error_message)

            logging.info('Loaded {} data from {}.'.format(maf_type, maf_sample_path))

            # Correct column names (add maf_type suffix)
            new_columns = []
            for column in input_mafs[maf_type].columns:
                new_column = column
                # Don't change merge column names
                if column not in self.merge_columns:
                    new_column = '{}_{}'.format(new_column, maf_type)

                new_columns.append(new_column)

            input_mafs[maf_type].columns = new_columns

        # Merge DFs on position into one agg DF.
        for maf_type, loaded_maf in input_mafs.items():
            if self.input_merge is None:
                self.input_merge = loaded_maf
            else:
                self.input_merge = self.input_merge.merge(loaded_maf, how='outer', on=self.merge_columns)

    @staticmethod
    def __cast_to_float_else_nan(value):
        try:
            count = float(value) if not pd.isnull(value) else 0
        except ValueError:
            count = float('nan')
        return count

    @staticmethod
    def __check_column_in_index(column_name, index, input_type):
        if column_name not in index:
            raise Exception('Error: Reference read count column %s not found in %s sample.' %
                            (column_name, input_type))

    @staticmethod
    def __calculate_vaf(ref_count, alt_count):
        if ref_count == 0 and alt_count == 0:
            return 0
        elif ref_count == float('nan') or alt_count == float('nan'):
            return float('nan')
        else:
            return float(alt_count) / float(ref_count + alt_count)

    def analyze(self):
        logging.info('Analyzing {} total variants...'.format(len(self.input_merge)))
        for maf_type in self.maf_types_for_analysis:
            allele1_k = 'Allele1_{}'.format(MAF_TYPE_PRINTABLE[maf_type])
            allele1_v = 'Tumor_Seq_Allele1_{}'.format(maf_type)
            allele2_k = 'Allele2_{}'.format(MAF_TYPE_PRINTABLE[maf_type])
            allele2_v = 'Tumor_Seq_Allele2_{}'.format(maf_type)
            OUTPUT_MAF_MAP[allele1_k] = allele1_v
            OUTPUT_MAF_MAP[allele2_k] = allele2_v

        ref_count_column_prefix = 'Ref_Read_Count'
        alt_count_column_prefix = 'Alt_Read_Count'
        vaf_column_prefix = 'VAF'

        # Add column mappings for all extra columns seen in all 4 files
        if self.extra_columns:
            for xc in self.extra_columns.split():
                for column_suffix in MAF_TYPE_PRINTABLE:
                    this_column_name = '{}_{}'.format(xc, MAF_TYPE_PRINTABLE[column_suffix])
                    OUTPUT_MAF_MAP[this_column_name] = '%s_%s' % (xc, column_suffix)

        # Add column mappings for extra columns seen in only one file
        for column_suffix in MAF_TYPE_PRINTABLE:
            this_suffix_columns = self.__getattribute__('{}_extra_columns'.format(column_suffix))[0]
            if this_suffix_columns:
                for column in this_suffix_columns.split():
                    this_column_name = '{}_{}'.format(column, MAF_TYPE_PRINTABLE[column_suffix])
                    OUTPUT_MAF_MAP[this_column_name] = '{}_{}'.format(column, column_suffix)

        # output_rows will associate the output column names with arrays that will contain that column's values per-row:
        output_rows = OrderedDict()
        expected_rows = OrderedDict()
        for output_column in ['scenario'] + OUTPUT_MAF_MAP.keys():
            output_rows[output_column] = []
            expected_rows[output_column] = []

        # Create read count arrays:
        for maf_type in self.maf_types_for_analysis:
            output_rows['{}_{}'.format(ref_count_column_prefix, MAF_TYPE_PRINTABLE[maf_type])] = []
            output_rows['{}_{}'.format(alt_count_column_prefix, MAF_TYPE_PRINTABLE[maf_type])] = []
            output_rows['{}_{}'.format(vaf_column_prefix, MAF_TYPE_PRINTABLE[maf_type])] = []
            expected_rows['{}_{}'.format(ref_count_column_prefix, MAF_TYPE_PRINTABLE[maf_type])] = []
            expected_rows['{}_{}'.format(alt_count_column_prefix, MAF_TYPE_PRINTABLE[maf_type])] = []
            expected_rows['{}_{}'.format(vaf_column_prefix, MAF_TYPE_PRINTABLE[maf_type])] = []

        scenario_calculator = ScenarioCalculator(self.scenarios_config_file_path)
        for i, row in self.input_merge.iterrows():
            # We exclude DNPs, TNPs, etc.
            if row['Variant_Type'] not in ['SNP', 'INS', 'DEL']:
                continue

            # Quad is a dictionary mapping each maf type to the corresponding VAFs on this row.
            counts = {}
            for input_type in self.maf_types_for_analysis:
                # Get column names
                ref_column = '{}_{}'.format(self.ref_count_column, input_type)
                alt_column = '{}_{}'.format(self.alt_count_column, input_type)
                # Check to make sure column names are in the sample's dataframe index
                self.__check_column_in_index(ref_column, row.index, input_type)
                self.__check_column_in_index(alt_column, row.index, input_type)

                # Cast ref and alt counts to float; if unable to do so, set as nan.
                ref_count = self.__cast_to_float_else_nan(row[ref_column])
                alt_count = self.__cast_to_float_else_nan(row[alt_column])

                # Calculate variant allele fraction
                vaf = self.__calculate_vaf(ref_count, alt_count)

                counts[input_type] = {'ref_count': ref_count, 'alt_count': alt_count, 'vaf': vaf}

            # row_dest either points to the output_rows or expected_rows OrderedDicts
            row_dest = None
            scenario_name = None
            # Construct VAF quad out of counts (pull vaf value out)
            vaf_quad = {}
            for input_type in self.maf_types_for_analysis:
                vaf_quad[input_type] = counts[input_type]['vaf']
            try:
                scenario_name = scenario_calculator.categorize(self.analysis_type, vaf_quad)
                row_dest = output_rows
            except ScenarioCalculator.NoScenarioException as e:
                logging.error("Could not find scenario for VAF values {} and analysis type {}".format(vaf_quad,
                                                                                                      analysis_type))

            row_dest['scenario'].append(scenario_name)
            for output_column in OUTPUT_MAF_MAP:
                input_column = OUTPUT_MAF_MAP[output_column]
                row_dest[output_column].append(row[input_column])

            for maf_type in self.maf_types_for_analysis:
                ref_output_column = '{}_{}'.format(ref_count_column_prefix, MAF_TYPE_PRINTABLE[maf_type])
                alt_output_column = '{}_{}'.format(alt_count_column_prefix, MAF_TYPE_PRINTABLE[maf_type])
                vaf_output_column = '{}_{}'.format(vaf_column_prefix, MAF_TYPE_PRINTABLE[maf_type])
                row_dest[ref_output_column].append(counts[maf_type]['ref_count'])
                row_dest[alt_output_column].append(counts[maf_type]['alt_count'])
                row_dest[vaf_output_column].append(counts[maf_type]['vaf'])

        # Write output.
        output_df = pd.DataFrame.from_dict(output_rows)

        # Prepend sample_id to each row if present
        if self.sample_id:
            output_df.insert(0, self.sample_id_header, self.sample_id)

        return output_df

    def output_results(self, output_df):
        output_df.to_csv(self.output, index=False, header=True, sep='\t')
        logging.info('Wrote {} discovered scenarios to {}.'.format(len(output_df.index), self.output.name))


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
                    output_location=args.o,
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
        output_df = r2d2.analyze()
        r2d2.output_results(output_df)
