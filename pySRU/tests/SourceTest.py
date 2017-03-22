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
from pySRU.ElectronBeam import ElectronBeam
from pySRU.SourceUndulatorPlane import SourceUndulatorPlane


class UndulatorParameterTest(unittest.TestCase):

    def test_copy(self):
        K = 1.87
        lambda_u = 0.035
        L = 0.035 * 14
        E=1.3
        I=1.0
        undulator = Undulator(K=K, period_length=lambda_u, length=L)
        beam=ElectronBeam(Electron_energy=E, I_current=I)
        source=SourceUndulatorPlane(undulator=undulator, electron_beam=beam)

        source2=source.copy()
        source2.electron_beam.I_current=0.3
        self.assertEqual(source.I_current(),1.0)
        self.assertEqual(source2.I_current() ,0.3)

        self.assertAlmostEqual(source.harmonic_frequency(1)/1e17, 2.5346701615509917,5)
        self.assertAlmostEqual(source.Lorentz_factor(), 2544.0367765521196,2 )
        self.assertAlmostEqual(source.electron_speed(), 0.99999992274559524,5)
        self.assertEqual(source.magnetic_structure.period_number(), 14)
        self.assertAlmostEqual(source.choose_distance_automatic(2),49.000000000000,10 )



    #TODO completer ?
