import unittest
import numpy as np
from pySRU.ParameterPlaneUndulator import ParameterPlaneUndulator as Undulator


class UndulatorParameterTest(unittest.TestCase):

    def test_copy_parameter(self):
        K = 1.87
        E = 1.3e9
        lambda_u = 0.035
        L = 0.035 * 12
        I = 1.0
        undulator = Undulator(K=K, E=E, lambda_u=lambda_u, L=L, I=I)
        undulator2 = undulator.copy()

        undulator2.I=0.2
        self.assertEqual(undulator.I,1.0)
        self.assertEqual(undulator2.I ,0.2)

        self.assertEqual(undulator.omega1(), 2.534659271012565e+17)
        self.assertEqual(undulator.gamma(), 2544.031311154599 )
        self.assertEqual(undulator.Beta(), 0.99999992274526328)
        self.assertEqual(undulator.Nb_period(), 12)
        self.assertEqual(undulator.D_max_plane_undulator(2),21.000000000000004 )


    def test_create_magnetic_field(self):
        K = 1.87
        E = 1.3e9
        lambda_u = 0.035
        L = 0.035 * 12
        I = 1.0
        undulator = Undulator(K=K, E=E, lambda_u=lambda_u, L=L, I=I)
        Z=np.linspace(-1.,1.,101)
        Y=0.0
        X=0.0
        B=undulator.create_magnetic_field(Z=Z,Y=Y,X=X,harmonic_number=1)
        By=B.By(Z,Y,X)
        self.assertTrue(By.shape==((101,)))

        Y=Z
        By=B.By(Z,Y,X)
        self.assertEqual(By.shape , ((101,)))