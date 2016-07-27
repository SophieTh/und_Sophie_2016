import unittest
import numpy as np
import scipy.constants as codata

from pySRU.MagneticStructureUndulatorPlane import MagneticStructureUndulatorPlane as Undulator
from pySRU.ElectronBeam import ElectronBeam
from pySRU.Source import Source
from pySRU.TrajectoryFactory import TrajectoryFactory,TRAJECTORY_METHOD_ANALYTIC,TRAJECTORY_METHOD_INTEGRATION,\
                                                        TRAJECTORY_METHOD_ODE
from pySRU.RadiationFactory import RadiationFactory , RADIATION_METHOD_APPROX_FARFIELD,RADIATION_METHOD_NEAR_FIELD ,\
                                RADIATION_METHOD_FARFIELD
from pySRU.Simulation import Simulation,create_simulation

class UndulatorSimulationTest(unittest.TestCase):


    def test_simulation(self):
        undulator_test = Undulator(K=1.87, lambda_u=0.035, L=0.035 * 14)
        electron_beam_test = ElectronBeam(E=1.3e9, I=1.0)

        sim_test = create_simulation(magnetic_structure=undulator_test,electron_beam=electron_beam_test,
                                     traj_method=TRAJECTORY_METHOD_ANALYTIC)
        source_test=sim_test.source

        self.assertFalse(all(sim_test.trajectory_fact.initial_condition==np.array([0.0,0.0,source_test.Beta()*codata.c,
                                                     0.0,0.0,source_test.Zo_symetry()])))
        ref=sim_test.copy()
        rad_max = sim_test.radiation.max()

        # test change
        sim_test.change_radiation_method(RADIATION_METHOD_FARFIELD)
        self.assertEqual(sim_test.radiation_fact.method, RADIATION_METHOD_FARFIELD)
        self.assertFalse(ref.radiation_fact.method==sim_test.radiation_fact.method)
        self.assertFalse(np.all(ref.radiation.intensity == sim_test.radiation.intensity))
        self.assertAlmostEqual(ref.radiation.intensity[0][0]/rad_max, sim_test.radiation.intensity[0][0]/rad_max, 3)

        sim_test.change_trajectory_method(TRAJECTORY_METHOD_INTEGRATION)
        self.assertEqual(sim_test.trajectory_fact.method, TRAJECTORY_METHOD_INTEGRATION)
        self.assertFalse(ref.trajectory_fact.method==sim_test.trajectory_fact.method)
        time_diff=np.abs(ref.trajectory.t - sim_test.trajectory.t)
        self.assertTrue(np.all(time_diff<=1e-19))
        self.assertFalse(np.all(ref.trajectory.x == sim_test.trajectory.x))
        self.assertFalse(np.all(ref.radiation.intensity == sim_test.radiation.intensity))
        rad_max = sim_test.radiation.max()
        self.assertAlmostEqual(ref.radiation.intensity[0][0]/rad_max, sim_test.radiation.intensity[0][0]/rad_max, 1)

        sim_test.change_Nb_pts_trajectory(ref.trajectory_fact.Nb_pts+1)
        self.assertEqual(sim_test.trajectory_fact.Nb_pts,ref.trajectory_fact.Nb_pts+1)
        self.assertEqual(sim_test.trajectory.nb_points(), ref.trajectory_fact.Nb_pts+1)
        self.assertFalse(ref.trajectory_fact.Nb_pts == sim_test.trajectory_fact.Nb_pts)
        self.assertAlmostEqual(ref.radiation.intensity[0][0]/rad_max,sim_test.radiation.intensity[0][0]/rad_max,1)

        sim_test.change_Nb_pts_radiation(100)
        self.assertEqual(sim_test.radiation_fact.Nb_pts,100)
        self.assertFalse(ref.radiation_fact.Nb_pts == sim_test.radiation_fact.Nb_pts)
        self.assertFalse(np.all(ref.radiation.X == sim_test.radiation.X))
        self.assertTrue(ref.radiation.X.min() == sim_test.radiation.X.min())
        self.assertTrue(ref.radiation.X.max() == sim_test.radiation.X.max())
        self.assertTrue(ref.radiation.Y.min() == sim_test.radiation.Y.min())
        #TODO pb here
        #self.assertTrue(ref.radiation.Y.max() == sim_test.radiation.Y.max())
        self.assertFalse(len(ref.radiation.X) == len(sim_test.radiation.X))

        sim_test.change_distance(50)
        self.assertEqual(sim_test.radiation.distance,50)
        self.assertFalse(ref.radiation.distance == sim_test.radiation.distance)

        sim_test.change_omega(source_test.omega1()*0.8)
        self.assertEqual(sim_test.radiation_fact.omega,source_test.omega1()*0.8)
        self.assertFalse(ref.radiation_fact.omega == sim_test.radiation_fact.omega)


        nb_pts=np.arange(500,2001,500,dtype='int')
        err=sim_test.error_radiation_method_nb_pts_traj(RADIATION_METHOD_APPROX_FARFIELD,nb_pts=nb_pts)
        self.assertLessEqual((err.max() / rad_max), 1e-2, 1)

        distance=np.arange(20,101,20,dtype='int')
        err = sim_test.error_radiation_method_distance(RADIATION_METHOD_APPROX_FARFIELD, D=distance)
        self.assertLessEqual(err.max()/rad_max, 1e-2, 1)


    #TODO


    # def test_error_ODE_analitic(self):
    #     und_test = Undulator(K=1.87, E=1.3e9, lambda_u=0.035, L=0.035 * 12, I=1.0)
    #     traj_test = TrajectoryFactory(Nb_pts=1000, method=TRAJECTORY_METHOD_ANALYTIC)
    #     rad_test = RadiationFactory(method=RADIATION_METHOD_APPROX_FARFIELD, omega=und_test.omega1(), Nb_pts=101)
    #     distance = 100
    #     Xmax = distance * 1e-3
    #     Ymax = distance * 1e-3
    #
    #     sim_test = create_simulation(parameter=und_test, trajectory_fact=traj_test, radiation_fact=rad_test,
    #                                  distance=distance, X_max=Xmax, Y_max=Ymax)
    #
    #     sim_test.change_trajectory_method(TRAJECTORY_METHOD_ANALYTIC)
    #     self.assertTrue(all(sim_test.trajectory_fact.initial_condition==
    #                      np.array([sim_test.trajectory.v_x[0],sim_test.trajectory.v_y[0],sim_test.trajectory.v_z[0],
    #                                sim_test.trajectory.x[0],sim_test.trajectory.y[0],sim_test.trajectory.z[0]])*codata.c))
    #     traj_ref = sim_test.trajectory.copy()
    #     rad_m=sim_test.radiation.max()
    #     nb_pts = np.arange(500, 2001, 500, dtype='int')
    #     rad_er,traj_er = sim_test.error_trajectory_method(TRAJECTORY_METHOD_ODE, nb_pts=nb_pts)
    #
    #     self.assertLessEqual(rad_er.max() / rad_m, 1e-2, 2)
    #     self.assertLessEqual(traj_er.x.max() / traj_ref.x.max(), 1e-2, 2)
    #     self.assertLessEqual(traj_er.z.max() / traj_ref.z.max(), 1e-2, 2)
    #     self.assertLessEqual(traj_er.v_x.max() / traj_ref.v_x.max(), 1e-2, 2)
    #     self.assertLessEqual(traj_er.v_z.max() / traj_ref.v_z.max(), 1e-2, 2)
    #     self.assertLessEqual(traj_er.a_x.max() / traj_ref.a_x.max(), 1e-2, 2)
    #     self.assertLessEqual(traj_er.a_z.max() / traj_ref.a_z.max(), 1e-2, 2)
