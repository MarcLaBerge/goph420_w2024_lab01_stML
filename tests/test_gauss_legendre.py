#Will contain unit tests for functions performing integration using
#Gauss-Legendre guadrature rules

from goph420_lab01.integration import integrate_gauss

import unittest
import numpy as np
import scipy as sc

#Testing Gauss with data from up to the 9th order polynomial
class TestGaussOrders(unittest.TestCase):
    

    def testGaussFirstOrder(self):
        #Test random 1st order function
        def f(x):
            return x + 2
        self.integral = integrate_gauss(f, lims = [0,2], npts = 1)
        #Expected value (excel calculator)
        expected = 6
        self.assertAlmostEqual(self.integral, expected)

    def testGaussSecondOrder(self):
        #Test random 2nd order function
        def f(x):
            return x ** 2 + 2
        self.integral = integrate_gauss(f, lims = [0,2], npts = 1)
        #Expected value (excel calculator)
        expected = 6
        self.assertAlmostEqual(self.integral, expected)
    
    def testGaussThirdOrder(self):
        #Test random 3rd order function
        def f(x):
            return x ** 3 + 2
        self.integral = integrate_gauss(f, lims = [0,2], npts = 2)
        #Expected value (excel calculator)
        expected = 8
        self.assertAlmostEqual(self.integral, expected)

    def testGaussFourthOrder(self):
        #Test random 4th order function
        def f(x):
            return x ** 4 + 2
        self.integral = integrate_gauss(f, lims = [0,2], npts = 2)
        #Expected value (excel calculator)
        expected = 92 / 9
        self.assertAlmostEqual(self.integral, expected)

    def testGaussFifthOrder(self):
        #Test random 5th order function
        def f(x):
            return x ** 5 + 2
        self.integral = integrate_gauss(f, lims = [0,2], npts = 3)
        #Expected value (excel calculator)
        expected = 132/9
        self.assertAlmostEqual(self.integral, expected)

    def testGaussSixthOrder(self):
        #Test random 6th order function
        def f(x):
            return x ** 6 + 2
        self.integral = integrate_gauss(f, lims = [0,2], npts = 3)
        #Expected value (excel calculator)
        expected = 556/25
        self.assertAlmostEqual(self.integral, expected)

    def testGaussSeventhOrder(self):
        #Test random 7th order function
        def f(x):
            return x ** 7 + 2
        self.integral = integrate_gauss(f, lims = [0,2], npts = 4)
        #Expected value (excel calculator)
        expected = 36
        self.assertAlmostEqual(self.integral, expected)

    def testGaussEighthOrder(self):
        #Test random 8th order function
        def f(x):
            return x ** 8 +2
        self.integral = integrate_gauss(f, lims = [0,2], npts = 4)
        #Expected value (excel calculator)
        expected = 60.877279
        self.assertAlmostEqual(self.integral, expected)

    def testGaussNinthOrder(self):
        #Test random 9th order function
        def f(x):
            return x ** 9
        self.integral = integrate_gauss(f, lims = [0,2], npts = 5)
        #Expected value (excel calculator)
        expected = 102.4
        self.assertAlmostEqual(self.integral, expected)