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
from pySRU.ElectronBeam import ElectronBeam
from pySRU.SourceUndulatorPlane import SourceUndulatorPlane
from pySRU.TrajectoryFactory import TrajectoryFactory,TRAJECTORY_METHOD_ANALYTIC, TRAJECTORY_METHOD_ODE
from pySRU.RadiationFactory import RadiationFactory , RADIATION_METHOD_APPROX_FARFIELD,RADIATION_METHOD_NEAR_FIELD

class RadiationFactoryTest(unittest.TestCase):

    #TODO des print sont cache qql part
    def test_create_radiation_undulator(self):
        undulator_test = Undulator(K=1.87, period_length=0.035, length=0.035 * 14)
        electron_beam_test = ElectronBeam(Electron_energy=1.3, I_current=1.0)
        source_test=SourceUndulatorPlane(undulator=undulator_test, electron_beam=electron_beam_test)
        traj_fact=TrajectoryFactory(Nb_pts=1001, method=TRAJECTORY_METHOD_ANALYTIC)
        traj=traj_fact.create_from_source(source_test)

        rad_fact = RadiationFactory(photon_frequency=source_test.harmonic_frequency(1), method=RADIATION_METHOD_NEAR_FIELD)

        rad=rad_fact.create_for_one_relativistic_electron(trajectory=traj, source=source_test)
        self.assertFalse(rad.X is None)
        self.assertFalse(rad.Y is None)
        self.assertFalse(rad.distance is None)





        rad_fact.method=RADIATION_METHOD_APPROX_FARFIELD

        rad2=rad_fact.create_for_one_relativistic_electron(trajectory=traj, source=source_test)
        self.assertTrue(rad2.distance == None)


        rad2=rad_fact.create_for_one_relativistic_electron(trajectory=traj, source=source_test,distance=rad.distance)
        self.assertFalse(rad.distance == None)
        err=rad.difference_with(rad2)

        self.assertTrue(rad.XY_are_similar_to(rad2))
        self.assertTrue(rad.XY_are_similar_to(err))
        self.assertTrue(rad.distance == rad2.distance)
        self.assertTrue(err.distance == rad2.distance)
        self.assertGreaterEqual(err.intensity.min(),0.0)
        self.assertLessEqual(err.max(), rad.max()*1e-1) # srio changed 1e-3 by 1e-1


        traj_test2=TrajectoryFactory(Nb_pts=1001, method=TRAJECTORY_METHOD_ODE,
                                     initial_condition=traj_fact.initial_condition).create_from_source(source_test)

        rad3=rad_fact.create_for_one_relativistic_electron(trajectory=traj_test2,source=source_test,distance=rad.distance)
        err = rad2.difference_with(rad3)
        self.assertLessEqual(err.max(),rad2.max()*1e-3)
