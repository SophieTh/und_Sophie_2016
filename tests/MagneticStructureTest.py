# coding: utf-8
# /*##########################################################################
#
# Copyright (c) 2016 European Synchrotron Radiation Facility
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#
# ###########################################################################*/

__authors__ = ["S Thery, M Glass, M Sanchez del Rio - ESRF ISDD Advanced Analysis and Modelling"]
__license__ = "MIT"
__date__ = "31/08/2016"

import unittest
import numpy as np
from pySRU.MagneticStructureUndulatorPlane import MagneticStructureUndulatorPlane as Undulator


class MagneticStructureTest(unittest.TestCase):

    def test_copy_undulator(self):
        K = 1.87
        period_length = 0.035
        length = 0.035 * 14
        undulator = Undulator(K=K, period_length=period_length, length=length)
        undulator2 = undulator.copy()

        undulator2.K=1.5
        self.assertEqual(undulator.K,1.87)
        self.assertEqual(undulator2.K ,1.5)

        self.assertEqual(undulator.period_number(), 14)



    def test_create_magnetic_field(self):
        K = 1.87
        lambda_u = 0.035
        L = 0.035 * 12
        undulator = Undulator(K=K, period_length=lambda_u, length=L)
        Z=np.linspace(-1.,1.,101)
        Y=0.0
        X=0.0
        B=undulator.create_magnetic_field(harmonic_number=1)
        By=B.By(Z,Y,X)
        self.assertTrue(By.shape==((101,)))

        Y=Z
        By=B.By(Z,Y,X)
        self.assertEqual(By.shape , ((101,)))



#TODO faire Bending magnet
