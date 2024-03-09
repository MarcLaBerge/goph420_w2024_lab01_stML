#Will contain unit tests for function performing integration using
#Newton-Cotes rules

from goph420_lab01.integration import integrate_newton

import unittest
import numpy as np

class TestNewtonCotes(unittest.TestCase):
    #Testing trap rule with linear functions
    def testTrapValue(self):
        self.x = [0,1,2,3]
        self.f = [0,1,2,3]
        self.newInt = integrate_newton(self.x, self.f, alg = "trap")
        #Check that they're equal
        expected = 4.5
        self.assertAlmostEqual(self.newInt, expected)

    #Test simp rule with odd number of data points
    def testSimpOdd(self):
        self.x = [0,1,2]
        self.f = -1 * self.x ** 2
        self.newInt = integrate_newton(self.x,self.f, alg = "simp")
        #Check that they're equal
        expected = -1 * (8 / 3)
        self.assertAlmostEqual(self.newInt, expected)

    #Test simp rule with even number of data points
    def testSimpEven(self):
        self.x = [0,1,2,3]
        self.f = -1 * self.x ** 2
        self.newInt = integrate_newton(self.x, self.f, alg = "simp")
        #Check that they're equal
        expected = - 9
        self.assertAlmostEqual(self.newInt, expected)





