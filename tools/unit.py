#!/usr/bin/env python

import unittest

class TestQuarticSolver(unittest.TestCase):

    def setUp(self):
        from rmhd import _lib
        self.solve_quartic_equation = _lib.solve_quartic_equation
        self.new_QuarticEquation    = _lib.new_QuarticEquation

    def test1(self):
        from ctypes import c_double, c_int, byref
        self.new_QuarticEquation(384, 1200, -480, -1620, 756)

        roots = [c_double() for i in range(4)]
        numrt = [c_int()    for i in range(2)]

        self.solve_quartic_equation(*tuple([byref(a) for a in roots+numrt]))

        test_soln = sorted([a.value for a in roots])
        true_soln = [-3.0, -3.0/2.0, 1.0/2.0, 7.0/8.0]

        for test, true in zip(test_soln, true_soln):
            self.assertAlmostEqual(test, true)


if __name__ == "__main__":

    unittest.main()
