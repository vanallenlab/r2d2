import unittest
from scenarios import ScenarioCalculator, Event
from r2d2 import AnalysisTypes

SCENARIOS_INI = './test_scenarios.ini.txt'


class TestScenarioCalculator(unittest.TestCase):
    def setUp(self):
        self.sc = ScenarioCalculator(SCENARIOS_INI)

    # Test not-categorizable scenarios
    def test_raises_no_scenario_exception(self):
        self.assertRaises(self.sc.NoScenarioException, self.sc.categorize, 'bad scenario', {})

    # Test that using incorrect VAF structure throws error
    def test_wrong_vaf_exception_1(self):
        self.assertRaises(self.sc.WrongVAFValuesException, self.sc.categorize, AnalysisTypes.tumor_only,
                          {'dna_normal': 0.25,
                           'dna_tumor': 0.20,
                           'rna_normal': 0.25,
                           'rna_tumor': 0.15})

    def test_wrong_vaf_exception_2(self):
        self.assertRaises(self.sc.WrongVAFValuesException, self.sc.categorize, AnalysisTypes.all_inputs,
                          {'dna_normal': 0.25,
                           'rna_tumor': 0.15})

    def test_wrong_vaf_exception_3(self):
        self.assertRaises(self.sc.WrongVAFValuesException, self.sc.categorize, AnalysisTypes.no_rna_normal,
                          {'dna_normal': 0.25,
                           'rna_tumor': 0.15})

    # Test all_inputs scenarios
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

    def test_all_inputs_vse(self):
        self.assertEqual(self.sc.categorize(AnalysisTypes.all_inputs, {'dna_normal': 0.55,
                                                                       'dna_tumor': 0.35,
                                                                       'rna_normal': 0.95,
                                                                       'rna_tumor': 0.90}), Event.vse)

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

    # Test no_rna_normal scenarios
    def test_no_rna_normal_germline_mosaic(self):
        self.assertEqual(self.sc.categorize(AnalysisTypes.no_rna_normal, {'dna_normal': .25,
                                                                          'dna_tumor': 0.20,
                                                                          'rna_tumor': 0.15}), Event.germline_mosaic)

    def test_no_rna_normal_tumor_in_normal(self):
        self.assertEqual(self.sc.categorize(AnalysisTypes.no_rna_normal, {'dna_normal': 0.25,
                                                                          'dna_tumor': 0.75,
                                                                          'rna_tumor': 0.75}), Event.tumor_in_normal)

    def test_no_rna_normal_vse(self):
        self.assertEqual(self.sc.categorize(AnalysisTypes.no_rna_normal, {'dna_normal': 0.55,
                                                                          'dna_tumor': 0.35,
                                                                          'rna_tumor': 0.90}), Event.vse)

    def test_no_rna_normal_vsl(self):
        self.assertEqual(self.sc.categorize(AnalysisTypes.no_rna_normal, {'dna_normal': 0.41,
                                                                          'dna_tumor': 0.69,
                                                                          'rna_tumor': 0.10}), Event.vsl)

    def test_no_rna_normal_loh_alt(self):
        self.assertEqual(self.sc.categorize(AnalysisTypes.no_rna_normal, {'dna_normal': 0.32,
                                                                          'dna_tumor': 0.96,
                                                                          'rna_tumor': 0.93}), Event.loh_alt)

    def test_no_rna_normal_loh_ref(self):
        self.assertEqual(self.sc.categorize(AnalysisTypes.no_rna_normal, {'dna_normal': 0.45,
                                                                          'dna_tumor': 0.08,
                                                                          'rna_tumor': 0.07}), Event.loh_ref)

    def test_no_rna_normal_germline(self):
        self.assertEqual(self.sc.categorize(AnalysisTypes.no_rna_normal, {'dna_normal': 0.01,
                                                                          'dna_tumor': 0.72,
                                                                          'rna_tumor': 0.03}), Event.germline)

    def test_no_rna_normal_rnaed(self):
        self.assertEqual(self.sc.categorize(AnalysisTypes.no_rna_normal, {'dna_normal': 0,
                                                                          'dna_tumor': 0,
                                                                          'rna_tumor': 0.45}), Event.rnaed)

    def test_no_rna_normal_somatic(self):
        self.assertEqual(self.sc.categorize(AnalysisTypes.no_rna_normal, {'dna_normal': 0,
                                                                          'dna_tumor': 0.20,
                                                                          'rna_tumor': 0.15}), Event.somatic)

    def test_no_rna_normal_unclassified(self):
        self.assertEqual(self.sc.categorize(AnalysisTypes.no_rna_normal, {'dna_normal': 0,
                                                                          'dna_tumor': 0,
                                                                          'rna_tumor': 0.15}), Event.unclassified)

    # Test dna only scenarios
    def test_dna_only_germline_mosaic(self):
        self.assertEqual(self.sc.categorize(AnalysisTypes.dna_only, {'dna_normal': .25,
                                                                     'dna_tumor': 0.20}), Event.germline_mosaic)

    def test_dna_only_tumor_in_normal(self):
        self.assertEqual(self.sc.categorize(AnalysisTypes.dna_only, {'dna_normal': 0.25,
                                                                     'dna_tumor': 0.75}), Event.tumor_in_normal)

    def test_dna_only_loh_alt(self):
        self.assertEqual(self.sc.categorize(AnalysisTypes.dna_only, {'dna_normal': 0.32,
                                                                     'dna_tumor': 0.96}), Event.loh_alt)

    def test_dna_only_loh_ref(self):
        self.assertEqual(self.sc.categorize(AnalysisTypes.dna_only, {'dna_normal': 0.45,
                                                                     'dna_tumor': 0.08}), Event.loh_ref)

    def test_dna_only_germline(self):
        self.assertEqual(self.sc.categorize(AnalysisTypes.dna_only, {'dna_normal': 0.4,
                                                                     'dna_tumor': 0.72}), Event.germline)

    def test_dna_only_somatic(self):
        self.assertEqual(self.sc.categorize(AnalysisTypes.dna_only, {'dna_normal': 0,
                                                                     'dna_tumor': 0.20}), Event.somatic)

    # Test normal only scenarios
    def test_normal_only_vse(self):
        self.assertEqual(self.sc.categorize(AnalysisTypes.normal_only, {'dna_normal': .35,
                                                                        'rna_normal': 0.98}), Event.vse)

    def test_normal_only_vsl(self):
        self.assertEqual(self.sc.categorize(AnalysisTypes.normal_only, {'dna_normal': .45,
                                                                        'rna_normal': 0.07}), Event.vsl)

    def test_normal_only_germline(self):
        self.assertEqual(self.sc.categorize(AnalysisTypes.normal_only, {'dna_normal': .35,
                                                                        'rna_normal': 0.88}), Event.germline)

    # Test tumor_only scenarios
    def test_tumor_only_vse(self):
        self.assertEqual(self.sc.categorize(AnalysisTypes.tumor_only, {'dna_tumor': .35,
                                                                       'rna_tumor': 0.98}), Event.vse)

    def test_tumor_only_vsl(self):
        self.assertEqual(self.sc.categorize(AnalysisTypes.tumor_only, {'dna_tumor': .45,
                                                                       'rna_tumor': 0.10}), Event.vsl)

    def test_tumor_only_somatic(self):
        self.assertEqual(self.sc.categorize(AnalysisTypes.tumor_only, {'dna_tumor': .25,
                                                                       'rna_tumor': 0.98}), Event.somatic)