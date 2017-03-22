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
import scipy.constants as codata
from pySRU.MagneticStructureUndulatorPlane import MagneticStructureUndulatorPlane as Undulator
from pySRU.SourceUndulatorPlane import SourceUndulatorPlane
from pySRU.ElectronBeam import ElectronBeam
from pySRU.TrajectoryFactory import TrajectoryFactory,TRAJECTORY_METHOD_ANALYTIC,   TRAJECTORY_METHOD_ODE

class TrajectoryFactoryTest(unittest.TestCase):

    #TODO a completer
    def test_create_trajectory(self):

        # test le trajectoire ana lytic (des valeurs special
        # le Beta doit est constant pour certaine methode
        # les produi scalaire de l'acc avec la vitesse est ...
        # difference max entre deux trajectoire

        K = 1.87
        lambda_u = 0.035
        L = 0.035 * 14
        E=1.3
        I=1.0
        undulator = Undulator(K=K, period_length=lambda_u, length=L)
        beam=ElectronBeam(Electron_energy=E, I_current=I)
        source=SourceUndulatorPlane(undulator=undulator, electron_beam=beam)

        fact_test = TrajectoryFactory(Nb_pts=201,method=TRAJECTORY_METHOD_ODE)
        traj_test = fact_test.create_from_source(source=source)

        self.assertFalse(fact_test.initial_condition is None)

        self.assertTrue(all(fact_test.initial_condition==source.choose_initial_contidion_automatic()))
        self.assertTrue(fact_test.method==TRAJECTORY_METHOD_ODE)

        scalar_product=traj_test.v_x*traj_test.a_x + traj_test.a_z*traj_test.v_z

        self.assertAlmostEqual(np.abs(scalar_product).max(),0.0,3)







