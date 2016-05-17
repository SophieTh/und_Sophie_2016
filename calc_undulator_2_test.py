import unittest
import numpy as np
import scipy.constants as codata
import scipy.integrate as integrate

from calc_undulator_2 import radiation, undulator_trajectory

class ClacUndulator2Test(unittest.TestCase):
    
    def test_radiation(self):
        K = 1.87
        E = 1.3e9
        lambda_u = 0.035
        Nb_period = 12
        Nb_points=100

        T = undulator_trajectory(K, E, lambda_u, Nb_period,Nb_points, type_trajectory=0)

        X = np.arange(0.0, 0.00101, 0.00001)
        Y = np.arange(0.0, 0.00101, 0.00001)
        XY = np.meshgrid(X, Y)
        Z2 = radiation(K=K, E=E, trajectory=T, X=XY[0], Y=XY[1])

        self.assertAlmostEqual(Z2.max() / 1e15, 1.314, 3)

    def test_trajectory(self):
        K = 1.87
        E = 1.3e9
        lambda_u = 0.035
        Nb_period = 12
        Nb_pts = 100

        # recuperation des donnees de B en array en fonction de z
        reference = np.load("x_ray_booklet_field.npz")
        Z_By = np.zeros((2, len(reference['ct'])))
        Z_By[0] = reference['ct']
        Z_By[0] -= (Z_By[0][len(Z_By[0]) - 1]) / 2.0
        Z_By[1] = reference['B_y']

        T = undulator_trajectory(K, E, lambda_u, Nb_period, Nb_pts, Z_By, type_trajectory=2)

        X = np.arange(0.0, 0.00101, 0.00001)
        Y = np.arange(0.0, 0.00101, 0.00001)
        XY = np.meshgrid(X,Y)
        Z2 = radiation(K=K, E=E, trajectory=T, X=XY[0], Y=XY[1])

        self.assertAlmostEqual(Z2.max()/1e15, 1.22700, 3)
