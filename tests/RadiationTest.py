import unittest
import numpy as np
from pySRU.RadiationGrid import RadiationGrid

class RadiationTest(unittest.TestCase):

    def test_first(self):

        X=np.linspace(0.0,0.5,101)
        Y=np.linspace(0.0,0.5,101)
        distance=100
        X_grid,Y_grid=np.meshgrid(X,Y)
        intensity = X_grid ** 2 * Y_grid ** 2
        rad=RadiationGrid(intensity=intensity, X=X_grid, Y=Y_grid,distance=distance)

        #TODO a verifier qd meme
        self.assertEqual(rad.max(), 0.0625)

        rad2=rad.copy()
        self.assertTrue(rad.XY_are_like_in(rad2))

        err=rad.error_max(rad2)

        self.assertEqual(err.max(),0.0)

        rad2.X += 1

        self.assertFalse(rad.XY_are_like_in(rad2))

        rad2.intensity= rad2.X ** 2 * rad2.Y ** 2

        err=rad.relativ_error(rad2)
        self.assertEqual(err,8.0)
