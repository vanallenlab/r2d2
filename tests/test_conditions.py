import unittest
import numpy as np
import os
import re
from collections import defaultdict
from scenarios import Condition
from r2d2 import get_analysis_type, AnalysisTypes


class TestCondition(unittest.TestCase):
    def test_incorrect_thresholds_throw_error(self):
        def setup_condition():
            Condition('= 2.5 3.8')
        self.assertRaises(ValueError, setup_condition)

    def test_incorrect_equality_formatting_throws_error(self):
        def setup_condition():
            Condition('= .8')
        self.assertRaises(ValueError, setup_condition)

    def test_condition_less_than(self):
        condition = Condition('< 0.3')

        # Should return true for inputs less than 0.3
        self.assertEqual(condition.test(0), True)
        self.assertEqual(condition.test(-1), True)
        self.assertEqual(condition.test(.1), True)
        self.assertEqual(condition.test(.2999), True)

        # Should return false for inputs greater than or equal to 0.3
        self.assertEqual(condition.test(.3), False)
        self.assertEqual(condition.test(1), False)
        self.assertEqual(condition.test(10000), False)

    def test_condition_less_than_or_equal_to(self):
        condition = Condition('<= 0.3')

        # Should return true for inputs less than or equal to 0.3
        self.assertEqual(condition.test(0), True)
        self.assertEqual(condition.test(-1), True)
        self.assertEqual(condition.test(.1), True)
        self.assertEqual(condition.test(.2999), True)
        self.assertEqual(condition.test(.3), True)

        # Should return false for inputs greater than 0.3
        self.assertEqual(condition.test(1), False)
        self.assertEqual(condition.test(10000), False)

    def test_condition_greater_than(self):
        condition = Condition('> 0.3')

        # Should return false for inputs less than or equal to 0.3
        self.assertEqual(condition.test(0), False)
        self.assertEqual(condition.test(-1), False)
        self.assertEqual(condition.test(.1), False)
        self.assertEqual(condition.test(.2999), False)
        self.assertEqual(condition.test(.3), False)

        # Should return true for inputs greater than 0.3
        self.assertEqual(condition.test(1), True)
        self.assertEqual(condition.test(10000), True)

    def test_condition_greater_than_or_equal_to(self):
        condition = Condition('>= 0.3')

        # Should return false for inputs less than or equal to 0.3
        self.assertEqual(condition.test(0), False)
        self.assertEqual(condition.test(-1), False)
        self.assertEqual(condition.test(.1), False)
        self.assertEqual(condition.test(.2999), False)

        # Should return true for inputs greater than 0.3
        self.assertEqual(condition.test(.3), True)
        self.assertEqual(condition.test(1), True)
        self.assertEqual(condition.test(10000), True)

    def test_condition_between(self):
        condition = Condition('= 0.4 0.6')

        # Should return false for inputs less than 0.4 or greater than 0.6
        self.assertEqual(condition.test(.3), False)
        self.assertEqual(condition.test(.7), False)

        # Should return true for inputs in between 0.4 and 0.6
        self.assertEqual(condition.test(.5), True)

        # Should return true for 0.4 and 0.6
        self.assertEqual(condition.test(.4), True)
        self.assertEqual(condition.test(.6), True)

    def test_condition_not_equal_to(self):
        condition = Condition('<> 0')

        # Should return true for inputs greater than 0
        self.assertEqual(condition.test(1), True)
        self.assertEqual(condition.test(.5), True)

        # Should return false for inputs equal to 0
        self.assertEqual(condition.test(0), False)
        self.assertEqual(condition.test(0.0), False)

    def test_condition_equal_to(self):
        condition = Condition('= 0.0 0.0')

        # Should return true for inputs equal to 0
        self.assertEqual(condition.test(0), True)
        self.assertEqual(condition.test(0.0), True)

    def test_condition_equal_to_5(self):
        condition = Condition('= 0.50 0.50')

        # Should return true for inputs equal to 0
        self.assertEqual(condition.test(.5), True)
        self.assertEqual(condition.test(0.50), True)

