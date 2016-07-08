import unittest
import numpy as np
import scipy.integrate as integrate
import scipy.constants as codata
import matplotlib.pyplot as plt
from pySRU.ParameterPlaneUndulator import ParameterPlaneUndulator as Undulator
from pySRU.ParameterBendingMagnet import BENDING_MAGNET as BM
from pySRU.TrajectoryFactory import TrajectoryFactory,TRAJECTORY_METHOD_ANALYTIC,TRAJECTORY_METHOD_INTEGRATION,\
                                                        TRAJECTORY_METHOD_ODE
from pySRU.RadiationFactory import RadiationFactory , RADIATION_METHOD_APPROX_FARFIELD,RADIATION_METHOD_NEAR_FIELD ,\
                                RADIATION_METHOD_FARFIELD
from pySRU.Simulation import Simulation,create_simulation

class MagneticFieldTest(unittest.TestCase):



    def create_magn_field_test(self, und_test, method_traj,formule):
        distance = 250
        theta_1_1 = und_test.theta(1, 1)
        Xmax = distance * theta_1_1
        Ymax = distance * theta_1_1

        traj_test = TrajectoryFactory(Nb_pts=2000, method=method_traj)
        rad_test = RadiationFactory(method=RADIATION_METHOD_APPROX_FARFIELD, omega=und_test.omega1(), Nb_pts=100,
                                    formula=formule)
        sim_test = create_simulation(parameter=und_test, trajectory_fact=traj_test, radiation_fact=rad_test,
                                     distance=distance, X_max=Xmax, Y_max=Ymax)
        if method_traj==TRAJECTORY_METHOD_ANALYTIC :
            Zo=und_test.Zo_analitic()
        else :
            Zo= und_test.Zo_symetry()
        self.assertEqual(sim_test.magnetic_field.z[0], Zo)
        self.assertEqual(sim_test.magnetic_field.z[-1], -Zo)
        self.assertEqual(len(sim_test.magnetic_field.z), sim_test.trajectory_fact.Nb_pts)
        Z_test = np.linspace(und_test.Zo_analitic(), -und_test.Zo_analitic(), sim_test.trajectory_fact.Nb_pts)
        diff_mf = np.abs(sim_test.magnetic_field.By(z=Z_test, y=0.0, x=0.0) -
                         und_test.Bo() * np.cos((2.0 * np.pi / und_test.lambda_u) * Z_test))
        self.assertTrue(all(diff_mf < np.abs(und_test.Bo()) * 1e-3))
        int1 = integrate.quad((lambda z: sim_test.magnetic_field.By(z=z, y=0.0, x=0.0)), sim_test.magnetic_field.z[0],
                              sim_test.magnetic_field.z[-1])
        self.assertAlmostEqual(int1[0], 0.0, 5)
        int2 = integrate.quad((lambda z1: integrate.quad((lambda z: sim_test.magnetic_field.By(z=z, y=0.0, x=0.0)),
                                                         sim_test.magnetic_field.z[0], z1)[0])
                              , sim_test.magnetic_field.z[0], sim_test.magnetic_field.z[-1])
        self.assertAlmostEqual(int2[0], 0.0, 5)

    # doen't work ...
    def create_magn_field_test_BM(self, BM, method_traj, formule=1):
        distance = 250.0
        Xmax = distance * 1e-3
        Ymax = distance * 1e-3

        traj_test = TrajectoryFactory(Nb_pts=2000, method=method_traj)
        rad_test = RadiationFactory(method=RADIATION_METHOD_APPROX_FARFIELD, omega=BM.omega1(), Nb_pts=100,
                                    formula=formule)
        sim_test = create_simulation(parameter=BM, trajectory_fact=traj_test, radiation_fact=rad_test,
                                     distance=distance, X_max=Xmax, Y_max=Ymax)
        if method_traj == TRAJECTORY_METHOD_ANALYTIC:
            Zo = BM.Zo_analitic()
        else:
            Zo = BM.Zo_symetry()
        self.assertEqual(sim_test.magnetic_field.z[0], Zo)
        self.assertEqual(sim_test.magnetic_field.z[-1], -Zo)
        self.assertEqual(len(sim_test.magnetic_field.z), sim_test.trajectory_fact.Nb_pts)
        Z_test = np.linspace(BM.Zo_analitic(), -BM.Zo_analitic(), sim_test.trajectory_fact.Nb_pts)
        diff_mf = np.abs(sim_test.magnetic_field.By(z=Z_test, y=0.0, x=0.0) -
                         BM.Bo * np.cos((2.0 * np.pi / BM.lambda_u) * Z_test))
        self.assertTrue(all(diff_mf < np.abs(BM.Bo) * 1e-3))
        int1 = integrate.quad((lambda z: sim_test.magnetic_field.By(z=z, y=0.0, x=0.0)),
                              sim_test.magnetic_field.z[0],
                              sim_test.magnetic_field.z[-1])
        self.assertAlmostEqual(int1[0], 0.0, 5)
        int2 = integrate.quad((lambda z1: integrate.quad((lambda z: sim_test.magnetic_field.By(z=z, y=0.0, x=0.0)),
                                                         sim_test.magnetic_field.z[0], z1)[0])
                              , sim_test.magnetic_field.z[0], sim_test.magnetic_field.z[-1])
        self.assertAlmostEqual(int2[0], 0.0, 5)

    def test_magn_field(self):
        und_test = Undulator(K=1.87, E=1.3e9, lambda_u=0.035, L=0.035 * 14, I=1.0)
        ESRF18 = Undulator(K=1.68, E=6.0e9, lambda_u=0.018, L=0.018 * 111, I=0.2)
        #ESRFBM = BM(E=6.0e9, Bo=0.8, div=5.0e-3, R=25.0, I=0.2)

        self.create_magn_field_test(und_test=und_test, method_traj=TRAJECTORY_METHOD_ODE,formule=1)
        print("und_test, ok")

        self.create_magn_field_test(und_test=ESRF18, method_traj=TRAJECTORY_METHOD_ODE,formule=1)
        print("esrf18, ok")

        # self.create_magn_field_test_BM(BM=ESRFBM, method_traj=TRAJECTORY_METHOD_ODE, formule=1)
        # print("esrf BM ok")