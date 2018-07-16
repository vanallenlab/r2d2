import unittest
import numpy as np
import os
import re
from collections import defaultdict
from scenarios import Condition, ScenarioCalculator, Event
from r2d2 import get_analysis_type, AnalysisTypes

SCENARIOS_INI = './test_scenarios.ini'


class TestScenarioCalculator(unittest.TestCase):
    def setUp(self):
        self.sc = ScenarioCalculator(SCENARIOS_INI)

    def test_raises_no_scenario_exception(self):
        self.assertRaises(self.sc.NoScenarioException, self.sc.categorize, 'bad scenario', {})

    def test_all_inputs_germline_mosaic(self):
        self.assertEqual(self.sc.categorize(AnalysisTypes.all_inputs, {'dna_normal': 0.25,
                                                                       'dna_tumor': 0.20,
                                                                       'rna_normal': 0.25,
                                                                       'rna_tumor': 0.15}), Event.germline_mosaic)

    def test_all_inputs_tumor_in_normal(self):
        self.assertEqual(self.sc.categorize(AnalysisTypes.all_inputs, {'dna_normal': 0.25,
                                                                       'dna_tumor': 0.75,
                                                                       'rna_normal': 0.25,
                                                                       'rna_tumor': 0.75}), Event.tumor_in_normal)

    def test_all_inputs_vse_all_inputs(self):
        self.assertEqual(self.sc.categorize(AnalysisTypes.all_inputs, {'dna_normal': 0.55,
                                                                       'dna_tumor': 0.35,
                                                                       'rna_normal': 0.95,
                                                                       'rna_tumor': 0.90}), Event.vse_all_inputs)

    def test_all_inputs_t_vse(self):
        self.assertEqual(self.sc.categorize(AnalysisTypes.all_inputs, {'dna_normal': 0.50,
                                                                       'dna_tumor': 0.48,
                                                                       'rna_normal': 0.10,
                                                                       'rna_tumor': 0.98}), Event.t_vse)

    def test_all_inputs_vsl(self):
        self.assertEqual(self.sc.categorize(AnalysisTypes.all_inputs, {'dna_normal': 0.41,
                                                                       'dna_tumor': 0.69,
                                                                       'rna_normal': 0.05,
                                                                       'rna_tumor': 0.10}), Event.vsl)

    def test_all_inputs_t_vsl(self):
        self.assertEqual(self.sc.categorize(AnalysisTypes.all_inputs, {'dna_normal': 0.58,
                                                                       'dna_tumor': 0.68,
                                                                       'rna_normal': 0.8,
                                                                       'rna_tumor': 0.02}), Event.t_vsl)

    def test_all_inputs_loh_alt(self):
        self.assertEqual(self.sc.categorize(AnalysisTypes.all_inputs, {'dna_normal': 0.32,
                                                                       'dna_tumor': 0.96,
                                                                       'rna_normal': 0.44,
                                                                       'rna_tumor': 0.93}), Event.loh_alt)

    def test_all_inputs_loh_ref(self):
        self.assertEqual(self.sc.categorize(AnalysisTypes.all_inputs, {'dna_normal': 0.45,
                                                                       'dna_tumor': 0.08,
                                                                       'rna_normal': 0.62,
                                                                       'rna_tumor': 0.07}), Event.loh_ref)

    def test_all_inputs_germline(self):
        self.assertEqual(self.sc.categorize(AnalysisTypes.all_inputs, {'dna_normal': 0.01,
                                                                       'dna_tumor': 0.72,
                                                                       'rna_normal': 0.59,
                                                                       'rna_tumor': 0.03}), Event.germline)

    def test_all_inputs_rnaed(self):
        self.assertEqual(self.sc.categorize(AnalysisTypes.all_inputs, {'dna_normal': 0,
                                                                       'dna_tumor': 0,
                                                                       'rna_normal': 0.35,
                                                                       'rna_tumor': 0.45}), Event.rnaed)

    def test_all_inputs_t_rnaed(self):
        self.assertEqual(self.sc.categorize(AnalysisTypes.all_inputs, {'dna_normal': 0,
                                                                       'dna_tumor': 0,
                                                                       'rna_normal': 0,
                                                                       'rna_tumor': 0.55}), Event.t_rnaed)

    def test_all_inputs_somatic(self):
        self.assertEqual(self.sc.categorize(AnalysisTypes.all_inputs, {'dna_normal': 0,
                                                                       'dna_tumor': 0.20,
                                                                       'rna_normal': 0,
                                                                       'rna_tumor': 0.15}), Event.somatic)

    def test_all_inputs_unclassified(self):
        self.assertEqual(self.sc.categorize(AnalysisTypes.all_inputs, {'dna_normal': 0,
                                                                       'dna_tumor': 0.20,
                                                                       'rna_normal': 0.25,
                                                                       'rna_tumor': 0.15}), Event.unclassified)