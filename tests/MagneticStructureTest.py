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
