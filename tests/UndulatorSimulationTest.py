import unittest
import numpy as np
import scipy.constants as codata
from pySRU.UndulatorParameter import UndulatorParameters as Undulator
from pySRU.TrajectoryFactory import TrajectoryFactory,TRAJECTORY_METHOD_ANALYTIC,TRAJECTORY_METHOD_INTEGRATION,\
                                                        TRAJECTORY_METHOD_ODE
from pySRU.RadiationFactory import RadiationFactory , RADIATION_METHOD_APPROX_FARFIELD,RADIATION_METHOD_NEAR_FIELD ,\
                                RADIATION_METHOD_FARFIELD
from pySRU.UndulatorSimulation import UndulatorSimulation,create_simulation

class UndulatorSimulationTest(unittest.TestCase):


    def test_simulation(self):
        und_test = Undulator(K=1.87, E=1.3e9, lambda_u=0.035, L=0.035 * 12, I=1.0)
        initial_cond = np.array(
            [0.0, 0.0, (und_test.Beta() * codata.c), 0.0, 0.0, -und_test.L / 2.0 - 5.0 * und_test.lambda_u])
        traj_test = TrajectoryFactory(Nb_pts=100, method=TRAJECTORY_METHOD_ANALYTIC,initial_condition=initial_cond)
        rad_test = RadiationFactory(method=RADIATION_METHOD_APPROX_FARFIELD, omega=und_test.omega1(), Nb_pts=101)
        distance = 100
        Xmax = distance * 1e-3
        Ymax = distance * 1e-3

        sim_test = create_simulation(undulator=und_test, trajectory_fact=traj_test,radiation_fact=rad_test,
                                     distance=distance, X_max=Xmax, Y_max=Ymax)

        self.assertFalse(np.all(sim_test.trajectory_fact.initial_condition==initial_cond))
        ref=sim_test.copy()

        # test change
        sim_test.change_radiation_method(RADIATION_METHOD_FARFIELD)
        self.assertEqual(sim_test.radiation_fact.method, RADIATION_METHOD_FARFIELD)
        self.assertFalse(ref.radiation_fact.method==sim_test.radiation_fact.method)
        self.assertFalse(np.all(ref.radiation.intensity == sim_test.radiation.intensity))
        self.assertAlmostEqual(ref.radiation.intensity[0][0]/1e15, sim_test.radiation.intensity[0][0]/1e15, 3)

        sim_test.change_trajectory_method(TRAJECTORY_METHOD_INTEGRATION)
        self.assertEqual(sim_test.trajectory_fact.method, TRAJECTORY_METHOD_INTEGRATION)
        self.assertFalse(ref.trajectory_fact.method==sim_test.trajectory_fact.method)
        time_diff=np.abs(ref.trajectory.t - sim_test.trajectory.t)
        self.assertTrue(np.all(time_diff<=1e-19))
        self.assertFalse(np.all(ref.trajectory.x == sim_test.trajectory.x))
        self.assertFalse(np.all(ref.radiation.intensity == sim_test.radiation.intensity))
        self.assertAlmostEqual(ref.radiation.intensity[0][0]/1e15, sim_test.radiation.intensity[0][0]/1e15, 1)

        sim_test.change_Nb_pts_trajectory(ref.trajectory_fact.Nb_pts+1)
        self.assertEqual(sim_test.trajectory_fact.Nb_pts,ref.trajectory_fact.Nb_pts+1)
        self.assertEqual(sim_test.trajectory.nb_points(), ref.trajectory_fact.Nb_pts+1)
        self.assertFalse(ref.trajectory_fact.Nb_pts == sim_test.trajectory_fact.Nb_pts)
        self.assertAlmostEqual(ref.radiation.intensity[0][0]/1e15,sim_test.radiation.intensity[0][0]/1e15,1)

        sim_test.change_Nb_pts_radiation(100)
        self.assertEqual(sim_test.radiation_fact.Nb_pts,100)
        self.assertFalse(ref.radiation_fact.Nb_pts == sim_test.radiation_fact.Nb_pts)
        self.assertFalse(np.all(ref.radiation.X == sim_test.radiation.X))
        self.assertTrue(ref.radiation.X[0] == sim_test.radiation.X[0])
        self.assertTrue(ref.radiation.X[-1] == sim_test.radiation.X[-1])
        self.assertFalse(len(ref.radiation.X) == len(sim_test.radiation.X))

        sim_test.change_distance(50)
        self.assertEqual(sim_test.radiation.distance,50)
        self.assertFalse(ref.radiation.distance == sim_test.radiation.distance)

        sim_test.change_omega(und_test.omega1()*0.8)
        self.assertEqual(sim_test.radiation_fact.omega,und_test.omega1()*0.8)
        self.assertFalse(ref.radiation_fact.omega == sim_test.radiation_fact.omega)


        nb_pts=np.arange(500,2001,500,dtype='int')
        err=sim_test.error_radiation_method_nb_pts_traj(RADIATION_METHOD_APPROX_FARFIELD,nb_pts=nb_pts)
        self.assertAlmostEqual(err.max()/1e9,2.94,1)

        distance=np.arange(20,101,20,dtype='int')
        err = sim_test.error_radiation_method_distance(RADIATION_METHOD_APPROX_FARFIELD, D=distance)
        self.assertAlmostEqual(err.max()/1e9, 6.65, 1)

        err_t,err_r = sim_test.error_trajectory_method(TRAJECTORY_METHOD_ANALYTIC, nb_pts=nb_pts)
        traj=ref.trajectory.copy()
        traj.convert(err_t)
        self.assertAlmostEqual(err_r.max(), 2000.0, 2)
        self.assertAlmostEqual(traj.x.max()/1e-14, 1.36, 2)
        self.assertAlmostEqual(traj.z.max()/1e-10, 7.00, 2)
        self.assertAlmostEqual(traj.v_x.max()/1e-4, 7.34, 2)
        self.assertAlmostEqual(traj.v_z.max(), 1.000, 2)
        self.assertAlmostEqual(traj.a_x.max()/1e7, 3.96, 2)
        self.assertAlmostEqual(traj.a_z.max()/1e4, 1.45, 2)
