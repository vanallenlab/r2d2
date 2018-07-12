import unittest
import numpy as np
import os
import re
from collections import defaultdict
from scenarios import Condition, Scenario
from r2d2 import get_analysis_type, AnalysisTypes


class TestR2D2(unittest.TestCase):
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