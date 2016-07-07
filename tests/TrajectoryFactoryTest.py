import unittest
import numpy as np
import scipy.constants as codata
from pySRU.ParameterPlaneUndulator import ParameterPlaneUndulator as Undulator
from pySRU.TrajectoryFactory import TrajectoryFactory,TRAJECTORY_METHOD_ANALYTIC,TRAJECTORY_METHOD_INTEGRATION,\
                                                        TRAJECTORY_METHOD_ODE

class TrajectoryFactoryTest(unittest.TestCase):

    #TODO
    def test_create_trajectory(self):

        # test le trajectoire ana lytic (des valeurs special
        # le Beta doit est constant pour certaine methode
        # les produi scalaire de l'acc avec la vitesse est ...
        # difference max entre deux trajectoire

        und_test = Undulator(K=1.87, E=1.3e9, lambda_u=0.035, L=0.035 * 12, I=1.0)
        fact_test = TrajectoryFactory(Nb_pts=201,method=TRAJECTORY_METHOD_INTEGRATION)
        traj_test = fact_test.create_for_plane_undulator_ideal(undulator=und_test)
        self.assertFalse(fact_test.initial_condition==None)
        self.assertTrue(fact_test.method==TRAJECTORY_METHOD_ANALYTIC)
        traj_test=traj_test.convert()
        ps=traj_test.v_x*traj_test.a_x + traj_test.a_z*traj_test.v_z

        self.assertLessEqual(np.abs(ps).max(),3.55e-3)

        Beta_et = 1.0 - (1.0 / (2.0 * und_test.gamma() ** 2)) * (1.0 + (und_test.K ** 2) / 2.0)

        ku = 2.0 * np.pi / und_test.lambda_u

        xo = und_test.K / (und_test.gamma() * Beta_et * ku)
        zo = Beta_et * codata.c * und_test.L / (2.0 * codata.c * Beta_et)
        vzo = Beta_et * codata.c - ((und_test.K / und_test.gamma()) ** 2)

        initial_condition = np.array([0.0, 0.0, vzo, xo, 0.0, zo])
        fact_test2= TrajectoryFactory(Nb_pts=201, method=TRAJECTORY_METHOD_INTEGRATION,
                                      initial_condition=initial_condition)
        B=und_test.create_magnetic_field(Z=traj_test.t*codata.c,Y=0.0,X=0.0,harmonic_number=1)

        traj_test2 = fact_test2.create_for_parameter( parameter=und_test,B=B)






