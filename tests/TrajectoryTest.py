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
from pySRU.Trajectory import Trajectory


class TrajectoryTest(unittest.TestCase):
    def test_error_nb_pts_copy(self):
        t = np.arange(0.0, 6.0*np.pi, 0.005*np.pi)
        x = np.sin(t)
        y = 0.0 * t
        z = np.cos(t)
        v_x = np.cos(t)
        v_y = 0.0 * t
        v_z = -np.sin(t)
        a_x = -np.sin(t)
        a_y = 0.0 * t
        a_z = -np.cos(t)
        Traj = Trajectory(t, x, y, z, v_x, v_y, v_z, a_x, a_y, a_z)

        Traj2 = Traj.copy()

        self.assertEqual(Traj.nb_points(), 1200)


        self.assertTrue(all(Traj.t==Traj2.t))

        error=Traj.error_max(Traj2)
        # Traj.draw_2_trajectory(Traj2)
        error[0] -= len(Traj.t)
        self.assertEqual(error.max(), 0.0)

        Traj2.multiply_by(3)
        error = Traj.error_max(Traj2)
        error[0] -= len(Traj.t)
        self.assertEqual(error.max(), 2.0)
