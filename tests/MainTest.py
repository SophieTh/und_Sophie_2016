import unittest
import numpy as np
import scipy.integrate as integrate
import scipy.constants as codata
import matplotlib.pyplot as plt
from pySRU.ParameterPlaneUndulator import ParameterPlaneUndulator as Undulator
from pySRU.TrajectoryFactory import TrajectoryFactory,TRAJECTORY_METHOD_ANALYTIC,TRAJECTORY_METHOD_INTEGRATION,\
                                                        TRAJECTORY_METHOD_ODE
from pySRU.RadiationFactory import RadiationFactory , RADIATION_METHOD_APPROX_FARFIELD,RADIATION_METHOD_NEAR_FIELD ,\
                                RADIATION_METHOD_FARFIELD,RADIATION_METHOD_APPROX
from pySRU.Simulation import Simulation,create_simulation

class MainTest(unittest.TestCase):

    def simul_undulator_near_to_farfield(self,und_test, method_traj,formule=1):

        distance = und_test.D_max_plane_undulator(2) * 2.0
        theta_max = und_test.theta(n=1, l=1)
        Xmax = distance * theta_max
        Ymax = distance * theta_max

        Nb_pts = int(2.0*und_test.lambda_u * 1e3 * und_test.Nb_period())


        traj_test = TrajectoryFactory(Nb_pts=Nb_pts, method=method_traj)
        rad_test = RadiationFactory(method=RADIATION_METHOD_NEAR_FIELD, omega=und_test.omega1(), Nb_pts=100,
                                    formula=formule)
        sim_test = create_simulation(parameter=und_test, trajectory_fact=traj_test, radiation_fact=rad_test,
                                     distance=distance, X_max=Xmax, Y_max=Ymax)

        ref=sim_test.copy()
        rad_max=(sim_test.radiation.max())

        sim_test.change_radiation_method(RADIATION_METHOD_APPROX)
        rad_err=ref.radiation.error_max(sim_test.radiation)
        self.assertLessEqual(np.abs(sim_test.radiation.max()-rad_max)/rad_max,5e-2)
        self.assertLessEqual(rad_err / rad_max, 5e-2)

        sim_test.change_radiation_method(RADIATION_METHOD_FARFIELD)
        rad_err=ref.radiation.error_max(sim_test.radiation)
        self.assertLessEqual(np.abs(sim_test.radiation.max()-rad_max)/rad_max,5e-2)
        self.assertLessEqual(rad_err  / rad_max, 5e-2)

        sim_test.change_radiation_method(RADIATION_METHOD_APPROX_FARFIELD)
        rad_err=ref.radiation.error_max(sim_test.radiation)
        self.assertLessEqual(np.abs(sim_test.radiation.max()-rad_max)/rad_max,5e-2)
        self.assertLessEqual(rad_err / rad_max, 5e-2)

    def simul_undulator_traj_method(self, und_test, method_rad, formule=1):
        distance = und_test.D_max_plane_undulator(2) * 2.0
        theta_max = und_test.theta(n=1, l=1)
        Xmax = distance * theta_max
        Ymax = distance * theta_max

        Nb_pts = int(2.0*und_test.lambda_u * 1e3 * und_test.Nb_period())

        traj_test = TrajectoryFactory(Nb_pts=Nb_pts, method=TRAJECTORY_METHOD_ANALYTIC)
        rad_test = RadiationFactory(method=method_rad, omega=und_test.omega1(), Nb_pts=100,
                                    formula=formule)
        sim_test = create_simulation(parameter=und_test, trajectory_fact=traj_test, radiation_fact=rad_test,
                                     distance=distance, X_max=Xmax, Y_max=Ymax)

        ref = sim_test.copy()
        rad_max = (sim_test.radiation.max())

        sim_test.change_trajectory_method(TRAJECTORY_METHOD_ODE)
        rad_err = ref.radiation.error_max(sim_test.radiation)
        traj_err = ref.trajectory.error_max(sim_test.trajectory)
        self.assertLessEqual(np.abs(sim_test.radiation.max() - rad_max) / rad_max, 1e-2)
        self.assertLessEqual(traj_err[1] / ref.trajectory.x.max(), 5e-3)
        self.assertLessEqual(traj_err[3] / ref.trajectory.z.max(), 5e-3)
        self.assertLessEqual(traj_err[4] / ref.trajectory.v_x.max(), 5e-3)
        self.assertLessEqual(traj_err[6] / ref.trajectory.v_z.max(), 5e-3)
        self.assertLessEqual(traj_err[7] / ref.trajectory.a_x.max(), 5e-3)
        self.assertLessEqual(traj_err[9] / ref.trajectory.a_z.max(), 5e-3)
        self.assertLessEqual(np.abs(sim_test.radiation.max()-rad_max)/rad_max,1e-2)
        self.assertLessEqual(rad_err / rad_max, 5e-2)

        sim_test.change_trajectory_method(TRAJECTORY_METHOD_INTEGRATION)
        rad_err = ref.radiation.error_max(sim_test.radiation)
        traj_err=ref.trajectory.error_max(sim_test.trajectory)
        self.assertLessEqual(np.abs(sim_test.radiation.max() - rad_max) / rad_max, 1e-2)
        self.assertLessEqual(traj_err[1] / (ref.trajectory.x.max()*und_test.L), 5e-2)
        self.assertLessEqual(traj_err[3] / ref.trajectory.z.max(), 5e-2)
        self.assertLessEqual(traj_err[4] / ref.trajectory.v_x.max(), 5e-2)
        self.assertLessEqual(traj_err[6] / ref.trajectory.v_z.max(), 5e-2)
        self.assertLessEqual(traj_err[7] / ref.trajectory.a_x.max(), 5e-2)
        self.assertLessEqual(traj_err[9] / ref.trajectory.a_z.max(), 5e-2)
        self.assertLessEqual(np.abs(sim_test.radiation.max() - rad_max) / rad_max, 1e-2)
        self.assertLessEqual(rad_err / rad_max, 1e-2)




    # simul_und_analitic_formule_near_farfield(und_test=und_test,formule=1)
    #TODO
    # def simul_BM_analitic_near_farfield(BM, formule):
    #     distance = 1000
    #     Xmax = distance * 1e-2
    #     Ymax = distance * 1e-2
    #
    #     traj_test = TrajectoryFactory(Nb_pts=2000, method=TRAJECTORY_METHOD_ANALYTIC)
    #     rad_test = RadiationFactory(method=RADIATION_METHOD_APPROX_FARFIELD, omega=BM.omega1(), Nb_pts=100,
    #                                 formula=formule)
    #     sim_test = create_simulation(parameter=BM, trajectory_fact=traj_test, radiation_fact=rad_test,
    #                                  distance=distance, X_max=Xmax, Y_max=Ymax)
    #     ref = sim_test.copy()
    #
    #     sim_test.change_trajectory_method(TRAJECTORY_METHOD_INTEGRATION)
    #     sim_test.trajectory.plot_2_trajectory(ref.trajectory)
    #     err = sim_test.trajectory.error_rel_max(ref.trajectory)
    #     err.plot()
    #     print(ref.radiation.max())
    #     err2 = ref.radiation.difference_with(sim_test.radiation)
    #     err2.plot()

    #TODO
    #
    # def simul_und_analitic_erreur_near_farfield(und_test, formule):
    #     distance = 100
    #     theta_1_1 = np.sqrt((1.0 + und_test.K ** 2 / 2.0)) * (1.0 / und_test.gamma())
    #     Xmax = distance * theta_1_1 / 4.0
    #     Ymax = distance * theta_1_1 / 4.0
    #
    #     traj_test = TrajectoryFactory(Nb_pts=2000, method=TRAJECTORY_METHOD_ANALYTIC)
    #     rad_test = RadiationFactory(method=RADIATION_METHOD_APPROX_FARFIELD, omega=und_test.omega1(), Nb_pts=100,
    #                                 formula=formule)
    #     sim_test = create_simulation(parameter=und_test, trajectory_fact=traj_test, radiation_fact=rad_test,
    #                                  distance=distance, X_max=Xmax, Y_max=Ymax)
    #     D = np.linspace(und_test.D_max_plane_undulator(2) * 0.5, und_test.D_max_plane_undulator(2) * 2.0, 10)
    #     err = sim_test.error_radiation_method_distance(RADIATION_METHOD_APPROX, D)
    #     plt.plot(D, err)
    #     plt.show()

    def test_main(self):
        und_test = Undulator(K=1.87, E=1.3e9, lambda_u=0.035, L=0.035 * 14, I=1.0)
        ESRF18 = Undulator(K=1.68, E=6.0e9, lambda_u=0.018, L=0.018 * 111, I=0.2)
        #ESRFBM = BM(E=6.0e9, Bo=0.8, div=5e-3, R=25.0, I=0.2)

        print(' undulator test ')
        self.simul_undulator_near_to_farfield(und_test=und_test,method_traj=TRAJECTORY_METHOD_ANALYTIC,formule=1)
        print("TRAJECTORY_METHOD_ANALYTIC ok")
        self.simul_undulator_near_to_farfield(und_test=und_test,method_traj=TRAJECTORY_METHOD_INTEGRATION, formule=1)
        print("TRAJECTORY_METHOD_INTEGRATION ok")
        self.simul_undulator_near_to_farfield(und_test=und_test,method_traj=TRAJECTORY_METHOD_ODE, formule=1)
        print("TRAJECTORY_METHOD_ODE ok")

        self.simul_undulator_traj_method(und_test=und_test, method_rad=RADIATION_METHOD_APPROX_FARFIELD)
        print('APPROX FARFIELD ok')

        self.simul_undulator_traj_method(und_test=und_test, method_rad=RADIATION_METHOD_FARFIELD)
        print('FARFIELD ok')

        self.simul_undulator_traj_method(und_test=und_test, method_rad=RADIATION_METHOD_APPROX)
        print('APPROX ok')

        self.simul_undulator_traj_method(und_test=und_test, method_rad=RADIATION_METHOD_NEAR_FIELD)
        print('NEAR FIELD ok')
        print(' ')
        print('undulator ESRF18')
        self.simul_undulator_near_to_farfield(und_test=ESRF18,method_traj=TRAJECTORY_METHOD_ANALYTIC, formule=1)
        print("TRAJECTORY_METHOD_ANALYTIC ok")


        self.simul_undulator_near_to_farfield(und_test=ESRF18,method_traj=TRAJECTORY_METHOD_INTEGRATION, formule=1)
        print("TRAJECTORY_METHOD_INTEGRATION ok")

        self.simul_undulator_near_to_farfield(und_test=ESRF18, method_traj=TRAJECTORY_METHOD_ODE, formule=1)
        print("TRAJECTORY_METHOD_ODE ok")

        self.simul_undulator_traj_method(und_test=ESRF18, method_rad=RADIATION_METHOD_APPROX_FARFIELD)
        print('APPROX FARFIELD ok')

        self.simul_undulator_traj_method(und_test=ESRF18, method_rad=RADIATION_METHOD_FARFIELD)
        print('FARFIELD ok')

        self.simul_undulator_traj_method(und_test=ESRF18, method_rad=RADIATION_METHOD_APPROX)
        print('APPROX ok')

        self.simul_undulator_traj_method(und_test=ESRF18, method_rad=RADIATION_METHOD_NEAR_FIELD)
        print('NEAR FIELD ok')