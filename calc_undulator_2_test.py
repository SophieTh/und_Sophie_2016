import unittest
import numpy as np
import scipy.constants as codata
import scipy.integrate as integrate

from calc_undulator_2 import radiation, trajectory_undulator_reference, trajectory_undulator

class ClacUndulator2Test(unittest.TestCase):
    def test_radiation(self):
        K = 1.87
        E = 1.3e9
        lambda_u = 0.035
        Nb_period = 10
        Nb_pts = 20

        gamma = E / 0.511e6
        Beta = np.sqrt(1.0 - (1.0 / gamma ** 2))
        Beta_et = Beta * (1.0 - (K / (2.0 * gamma)) ** 2)

        T = trajectory_undulator_reference(K=K, gamma=gamma, lambda_u=lambda_u, Nb_period=Nb_period, Nb_point=Nb_pts,
                                           Beta_et=Beta_et)

        X = np.arange(0.0, 0.00101, 0.00001)
        Y = np.arange(0.0, 0.00101, 0.00001)
        Z2 = radiation(K=K, E=E, trajectory=T, X=X, Y=Y)

        self.assertAlmostEqual(Z2.max()/1e14, 9.1716, 3)
