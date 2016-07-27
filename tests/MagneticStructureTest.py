import unittest
import numpy as np
from pySRU.MagneticStructureUndulatorPlane import MagneticStructureUndulatorPlane as Undulator


class UndulatorParameterTest(unittest.TestCase):

    def test_copy_parameter(self):
        K = 1.87
        lambda_u = 0.035
        L = 0.035 * 12
        undulator = Undulator(K=K, lambda_u=lambda_u, L=L)
        undulator2 = undulator.copy()

        undulator2.K=1.5
        self.assertEqual(undulator.K,1.87)
        self.assertEqual(undulator2.K ,1.5)

        # self.assertEqual(undulator.omega1(), 2.534659271012565e+17)
        # self.assertEqual(undulator.gamma(), 2544.031311154599 )
        # self.assertEqual(undulator.Beta(), 0.99999992274526328)
        self.assertEqual(undulator.Nb_period(), 12)
        self.assertEqual(undulator.D_min(2),21.000000000000004 )


    def test_create_magnetic_field(self):
        K = 1.87
        lambda_u = 0.035
        L = 0.035 * 12
        undulator = Undulator(K=K, lambda_u=lambda_u, L=L)
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
