import unittest
import numpy as np
from pySRU.Radiation import Radiation

class RadiationTest(unittest.TestCase):

    def test_first(self):

        X=np.linspace(0.0,0.5,101)
        Y=np.linspace(0.0,0.5,101)
        distance=100
        X_grid,Y_grid=np.meshgrid(X,Y)
        intensity = X_grid ** 2 * Y_grid ** 2
        rad=Radiation(intensity=intensity, X=X, Y=Y,distance=distance)

        self.assertEqual(rad.max(), 0.0625)

        rad2=rad.copy()
        err=rad.error_max(rad2)

        self.assertEqual(err.max(),0.0)

        rad2.X += 1

        self.assertFalse(np.all(rad.X==rad2.X))

        rad2.intensity= rad2.X ** 2 * rad2.Y ** 2

        err=rad.relativ_error(rad2)
        self.assertEqual(err,9.0)
