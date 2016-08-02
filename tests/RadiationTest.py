import unittest
import numpy as np
from pySRU.RadiationGrid import RadiationGrid
from pySRU.RadiationList import RadiationList

class RadiationTest(unittest.TestCase):

    def test_first(self):

        X=np.linspace(0.0,0.5,101)
        Y=np.linspace(0.0,0.5,101)
        distance=100
        X_grid,Y_grid=np.meshgrid(X,Y)
        intensity_gird = X_grid ** 2 * Y_grid ** 2
        intensity_list= X ** 2 * Y ** 2
        radGrid=RadiationGrid(intensity=intensity_gird, X=X_grid, Y=Y_grid,distance=distance)
        radList=RadiationList(intensity=intensity_list, X=X, Y=Y,distance=distance)

        self.assertEqual(radGrid.max(), radList.max())


        # test RadiationGrid

        # test copy
        radGrid2=radGrid.copy()
        self.assertTrue(radGrid.XY_are_like_in(radGrid2))
        err=radGrid.error_max(radGrid2)
        self.assertEqual(err.max(),0.0)

        radGrid2.X += 1
        self.assertFalse(radGrid.XY_are_like_in(radGrid2))

        radGrid2.intensity= radGrid2.X ** 2 * radGrid2.Y ** 2

        err=radGrid.relativ_error(radGrid2)
        self.assertEqual(err,8.0)

        # test RadiationList

        # test copy
        radList2=radList.copy()
        self.assertTrue(radList.XY_are_like_in(radList2))
        err = radList.error_max(radList2)
        self.assertEqual(err.max(), 0.0)

        radList2.X += 1
        self.assertFalse(radList.XY_are_like_in(radList2))

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
        radGrid = RadiationGrid(intensity=intensity_gird, X=X_grid, Y=Y_grid, distance=distance)
        radList = RadiationList(intensity=intensity_list, X=X, Y=Y, distance=distance)

        integartionGrid=radGrid.integration()
        integrationList=radList.integration()
        self.assertEqual(integrationList,0.5)
        self.assertEqual(integartionGrid,0.5+integrationList)

