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
from scipy import integrate
from pySRU.Radiation import Radiation


class RadiationTest(unittest.TestCase):

    def test_first(self):

        X=np.linspace(0.0,0.5,101)
        Y=np.linspace(0.0,0.5,101)
        distance=100
        X_grid,Y_grid=np.meshgrid(X,Y)
        intensity_gird = X_grid ** 2 * Y_grid ** 2
        intensity_list= X ** 2 * Y ** 2
        radGrid=Radiation(intensity=intensity_gird, X=X_grid, Y=Y_grid,distance=distance)
        radList=Radiation(intensity=intensity_list, X=X, Y=Y,distance=distance)

        self.assertEqual(radGrid.max(), radList.max())


        # test RadiationGrid

        # test copy
        radGrid2=radGrid.copy()
        self.assertTrue(radGrid.XY_are_similar_to(radGrid2))

        err=radGrid.error_max(radGrid2)
        self.assertEqual(err.max(),0.0)

        # test modification of the screen
        radGrid2.X += 1
        self.assertFalse(radGrid.XY_are_similar_to(radGrid2))

        radGrid3=radGrid.copy()
        radGrid3.change_Nb_pts(Nb_pts=120)
        self.assertFalse(radGrid.XY_are_similar_to(radGrid3))
        self.assertTrue(len(radGrid3.X.shape)==2)
        self.assertTrue(len(radGrid3.Y.shape) == 2)
        self.assertTrue(radGrid.X.max()==radGrid3.X.max())
        self.assertTrue(radGrid.X.min() == radGrid3.X.min())
        self.assertTrue(radGrid.Y.max() == radGrid3.Y.max())
        self.assertTrue(radGrid.Y.min() == radGrid3.Y.min())

        #error
        radGrid2.intensity= radGrid2.X ** 2 * radGrid2.Y ** 2

        err=radGrid.relativ_error(radGrid2)
        self.assertEqual(err,8.0)




        # test RadiationList

        # test copy
        radList2=radList.copy()
        self.assertTrue(radList.XY_are_similar_to(radList2))

        err = radList.error_max(radList2)
        self.assertEqual(err.max(), 0.0)

        # test modification of the screen
        radList2.X += 1
        self.assertFalse(radList.XY_are_similar_to(radList2))

        radList3=radList.copy()
        radList3.change_Nb_pts(Nb_pts=120)
        self.assertFalse(radList.XY_are_similar_to(radList3))
        self.assertTrue(len(radList3.X.shape)==1)
        self.assertTrue(len(radList3.Y.shape) == 1 )
        self.assertTrue(radList.X.max()==radList3.X.max())
        self.assertTrue(radList.X.min() == radList3.X.min())
        self.assertTrue(radList.Y.max() == radList3.Y.max())
        self.assertTrue(radList.Y.min() == radList3.Y.min())


        #test error
        radList2.intensity = radList2.X ** 2 * radList2.Y ** 2

        err=radGrid.relativ_error(radGrid2)
        self.assertEqual(err,8.0)



    def test_integration(self):
        X = np.linspace(0.0, 1., 101)
        Y = np.linspace(0.0, 1., 101)
        distance = 100
        X_grid, Y_grid = np.meshgrid(X, Y)
        intensity_gird = X_grid + Y_grid
        X *= 0.0
        intensity_list = X + Y
        radGrid = Radiation(intensity=intensity_gird, X=X_grid, Y=Y_grid, distance=distance)
        radList = Radiation(intensity=intensity_list, X=X, Y=Y, distance=distance)

        integartionGrid=radGrid.integration(use_flux_per_mrad2_or_mm2=0)
        integrationList=radList.integration(use_flux_per_mrad2_or_mm2=0)
        self.assertEqual(integrationList,0.5)
        self.assertEqual(integartionGrid,0.5+integrationList)

