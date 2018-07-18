import unittest
import numpy as np
import os
import re
from collections import defaultdict
from scenarios import Condition, ScenarioCalculator
from r2d2 import get_analysis_type, R2D2, R2D2ParsingException
from maf_types import AnalysisTypes


class TestR2D2(unittest.TestCase):
    def test_r2d2_init(self):
        dn = 'a'
        dt = 'b'
        rn = 'c'
        rt = 'd'
        r2d2 = R2D2(analysis_type=AnalysisTypes.all_inputs, dna_normal=dn, dna_tumor=dt, rna_normal=rn, rna_tumor=rt,
                    load_data = False)
        self.assertEqual(r2d2.maf_sample_paths, {
            'dna_normal': 'a',
            'dna_tumor': 'b',
            'rna_normal': 'c',
            'rna_tumor': 'd'
        })

        dn = 'a'
        dt = 'b'
        rt = 'd'
        r2d2 = R2D2(analysis_type=AnalysisTypes.no_rna_normal, dna_normal=dn, dna_tumor=dt, rna_tumor=rt,
                    load_data=False)
        self.assertEqual(r2d2.maf_sample_paths, {
            'dna_normal': 'a',
            'dna_tumor': 'b',
            'rna_tumor': 'd'
        })

        dn = 'a'
        rn = 'c'
        r2d2 = R2D2(analysis_type=AnalysisTypes.normal_only, dna_normal=dn, rna_normal=rn, load_data=False)
        self.assertEqual(r2d2.maf_sample_paths, {
            'dna_normal': 'a',
            'rna_normal': 'c'
        })

        dn = 'a'
        dt = 'b'
        r2d2 = R2D2(analysis_type=AnalysisTypes.dna_only, dna_normal=dn, dna_tumor=dt, load_data=False)
        self.assertEqual(r2d2.maf_sample_paths, {
            'dna_normal': 'a',
            'dna_tumor': 'b'
        })

        dt = 'b'
        rt = 'd'
        r2d2 = R2D2(analysis_type=AnalysisTypes.tumor_only, dna_tumor=dt, rna_tumor=rt, load_data=False)
        self.assertEqual(r2d2.maf_sample_paths, {
            'dna_tumor': 'b',
            'rna_tumor': 'd'
        })

    def test_r2d2_load_data(self):
        dn = '../test_input/scenarios/loh-amp/dna_normal.maf'
        dt = '../test_input/scenarios/loh-amp/dna_tumor.maf'
        rn = '../test_input/scenarios/loh-amp/rna_normal.maf'
        rt = '../test_input/scenarios/loh-amp/rna_tumor.maf'
        r2d2 = R2D2(analysis_type=AnalysisTypes.all_inputs, dna_normal=dn, dna_tumor=dt, rna_normal=rn, rna_tumor=rt)
        r2d2.load_data()

        # Check that there are equal numbers of columns in the dataframe with the dna_normal, dna_tumor, rna_normal, and
        # rna_tumor suffixes
        column_names = r2d2.input_merge.keys()
        dn_cols = [c for c in column_names if c.endswith('dna_normal')]
        dt_cols = [c for c in column_names if c.endswith('dna_tumor')]
        rn_cols = [c for c in column_names if c.endswith('rna_normal')]
        rt_cols = [c for c in column_names if c.endswith('rna_tumor')]

        self.assertEqual(len(dn_cols), len(dt_cols))
        self.assertEqual(len(dn_cols), len(rn_cols))
        self.assertEqual(len(dn_cols), len(rt_cols))

    def test_r2d2_load_data_merge_column_missing(self):
        dn = '../test_input/scenarios/loh-del_dna-normal-missing-merge-column/dna_normal.maf'
        dt = '../test_input/scenarios/loh-del_dna-normal-missing-merge-column/dna_tumor.maf'
        rn = '../test_input/scenarios/loh-del_dna-normal-missing-merge-column/rna_normal.maf'
        rt = '../test_input/scenarios/loh-del_dna-normal-missing-merge-column/rna_tumor.maf'
        r2d2 = R2D2(analysis_type=AnalysisTypes.all_inputs, dna_normal=dn, dna_tumor=dt, rna_normal=rn, rna_tumor=rt,
                    load_data=False)
        self.assertRaises(R2D2ParsingException, r2d2.load_data)

    def test_r2d2_analyze_all_inputs_germline_mosaic(self):
        dn = '../test_data/all_inputs/germline_mosaic/small-dna-normal.maf'
        dt = '../test_data/all_inputs/germline_mosaic/small-dna-tumor.maf'
        rn = '../test_data/all_inputs/germline_mosaic/small-rna-normal.maf'
        rt = '../test_data/all_inputs/germline_mosaic/small-rna-tumor.maf'
        r2d2 = R2D2(analysis_type=AnalysisTypes.all_inputs, dna_normal=dn, dna_tumor=dt, rna_normal=rn, rna_tumor=rt)
        r2d2.analyze()

    # Test that the correct analysis types are returned
    def test_get_analysis_type_4_files(self):
        dn = 'a'
        dt = 'b'
        rn = 'c'
        rt = 'd'

        self.assertEqual(get_analysis_type(dn, dt, rn, rt), [AnalysisTypes.all_inputs])

    def test_get_analysis_type_3_files_canonical(self):
        dn = 'a'
        dt = 'b'
        rn = None
        rt = 'd'

        self.assertEqual(get_analysis_type(dn, dt, rn, rt), [AnalysisTypes.no_rna_normal])

    def test_get_analysis_type_3_files_other(self):
        dn = 'a'
        dt = None
        rn = 'c'
        rt = 'd'

        self.assertEqual(get_analysis_type(dn, dt, rn, rt), [AnalysisTypes.normal_only])

        dn = None
        dt = 'b'
        rn = 'c'
        rt = 'd'

        self.assertEqual(get_analysis_type(dn, dt, rn, rt), [AnalysisTypes.tumor_only])

        dn = 'a'
        dt = 'b'
        rn = 'c'
        rt = None

        self.assertEqual(get_analysis_type(dn, dt, rn, rt), [AnalysisTypes.dna_only, AnalysisTypes.normal_only])

    def test_get_analysis_type_2_files(self):
        dn = 'a'
        dt = 'b'
        rn = None
        rt = None

        self.assertEqual(get_analysis_type(dn, dt, rn, rt), [AnalysisTypes.dna_only])

        dn = None
        dt = 'b'
        rn = None
        rt = 'd'

        self.assertEqual(get_analysis_type(dn, dt, rn, rt), [AnalysisTypes.tumor_only])

        dn = 'a'
        dt = None
        rn = 'c'
        rt = None

        self.assertEqual(get_analysis_type(dn, dt, rn, rt), [AnalysisTypes.normal_only])

    def test_get_analysis_type_1_file(self):
        dn = 'a'
        dt = None
        rn = None
        rt = None

        self.assertEqual(get_analysis_type(dn, dt, rn, rt), [])

    def test_get_analysis_type_0_files(self):
        dn = None
        dt = None
        rn = None
        rt = None

        self.assertEqual(get_analysis_type(dn, dt, rn, rt), [])