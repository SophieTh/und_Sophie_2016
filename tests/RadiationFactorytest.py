import unittest
import numpy as np
import scipy.constants as codata
from pySRU.ParameterPlaneUndulator import ParameterPlaneUndulator as Undulator
from pySRU.TrajectoryFactory import TrajectoryFactory,TRAJECTORY_METHOD_ANALYTIC,TRAJECTORY_METHOD_INTEGRATION,\
                                                        TRAJECTORY_METHOD_ODE
from pySRU.RadiationFactory import RadiationFactory , RADIATION_METHOD_APPROX_FARFIELD,RADIATION_METHOD_NEAR_FIELD ,\
                                RADIATION_METHOD_FARFIELD

class RadiationFactoryTest(unittest.TestCase):


    def test_create_radiation(self):
        und_test = Undulator(K=1.87, E=1.3e9, lambda_u=0.035, L=0.035 * 12, I=1.0)
        traj_test = TrajectoryFactory(Nb_pts=1001, method=TRAJECTORY_METHOD_ANALYTIC).create_for_plane_undulator_ideal(
            undulator=und_test)
        traj_test=traj_test.convert()
        fact_rad = RadiationFactory(omega=und_test.omega1(), method=RADIATION_METHOD_APPROX_FARFIELD, Nb_pts=101)

        rad=fact_rad.create_for_single_electron(trajectory=traj_test, parameter=und_test)
        self.assertFalse(rad.X == None)
        self.assertFalse(rad.Y == None)
        self.assertFalse(rad.distance == None)


        fact_rad.method=RADIATION_METHOD_FARFIELD

        rad2=fact_rad.create_for_single_electron(trajectory=traj_test, parameter=und_test)
        err=rad.difference_with(rad2)

        self.assertTrue(np.all(rad.X==rad2.X))
        self.assertTrue(np.all(err.X == rad2.X))
        self.assertTrue(np.all(rad.distance == rad2.distance))
        self.assertTrue(np.all(err.distance == rad2.distance))
        self.assertGreaterEqual(err.intensity.min(),0.0)
        self.assertLessEqual(err.max(), 1.9e13)

        B=und_test.create_magnetic_field(Z=traj_test.z*codata.c,Y=0.0,X=0.0,harmonic_number=1)
        traj_test2=TrajectoryFactory(Nb_pts=1001, method=TRAJECTORY_METHOD_INTEGRATION).create_for_parameter(
            parameter=und_test,B=B )
        rad3=fact_rad.create_for_single_electron(trajectory=traj_test2, parameter=und_test)

        err = rad2.difference_with(rad3)
        self.assertLessEqual(err.max(), 1.1e11)
