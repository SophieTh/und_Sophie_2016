import unittest
import numpy as np
import scipy.constants as codata

from pySRU.MagneticStructureUndulatorPlane import MagneticStructureUndulatorPlane as Undulator
from pySRU.ElectronBeam import ElectronBeam
from pySRU.Source import Source
from pySRU.TrajectoryFactory import TrajectoryFactory,TRAJECTORY_METHOD_ANALYTIC,    TRAJECTORY_METHOD_ODE
from pySRU.RadiationFactory import RadiationFactory , RADIATION_METHOD_APPROX_FARFIELD,RADIATION_METHOD_NEAR_FIELD
from pySRU.Simulation import Simulation,create_simulation

class UndulatorSimulationTest(unittest.TestCase):


    def test_simulation(self):
        electron_beam_test = ElectronBeam(Electron_energy=1.3, I_current=1.0)
        beam_ESRF = ElectronBeam(Electron_energy=6.0, I_current=0.2)
        undulator_test = Undulator(K=1.87, period_length=0.035, length=0.035 * 14)
        ESRF18 = Undulator(K=1.68, period_length=0.018, length=2.0)

        sim_test = create_simulation(magnetic_structure=undulator_test,electron_beam=electron_beam_test,
                                     traj_method=TRAJECTORY_METHOD_ANALYTIC,rad_method=RADIATION_METHOD_NEAR_FIELD)
        self.assertFalse(sim_test.radiation.distance == None)
        source_test=sim_test.source

        self.assertFalse(all(sim_test.trajectory_fact.initial_condition==
                             source_test.choose_initial_contidion_automatic()))
        ref=sim_test.copy()
        rad_max = sim_test.radiation.max()

        # test change
        sim_test.change_radiation_method(RADIATION_METHOD_APPROX_FARFIELD)
        self.assertEqual(sim_test.radiation_fact.method, RADIATION_METHOD_APPROX_FARFIELD)
        self.assertFalse(ref.radiation_fact.method==sim_test.radiation_fact.method)
        self.assertFalse(np.all(ref.radiation.intensity == sim_test.radiation.intensity))
        self.assertAlmostEqual(ref.radiation.intensity[0][0]/rad_max, sim_test.radiation.intensity[0][0]/rad_max, 3)

        sim_test.change_trajectory_method(TRAJECTORY_METHOD_ODE)
        self.assertEqual(sim_test.trajectory_fact.method, TRAJECTORY_METHOD_ODE)
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
        self.assertTrue(ref.radiation.Y.max() == sim_test.radiation.Y.max())
        self.assertFalse(len(ref.radiation.X) == len(sim_test.radiation.X))

        sim_test.change_distance(50)
        self.assertEqual(sim_test.radiation.distance,50)
        self.assertFalse(ref.radiation.distance == sim_test.radiation.distance)

        sim_test.change_omega(source_test.harmonic_frequency(1)*0.8)
        self.assertEqual(sim_test.radiation_fact.omega,source_test.harmonic_frequency(1)*0.8)
        self.assertFalse(ref.radiation_fact.omega == sim_test.radiation_fact.omega)


        # nb_pts=np.arange(500,2001,500,dtype='int')
        # err=sim_test.error_radiation_method_nb_pts_traj(RADIATION_METHOD_APPROX_FARFIELD,nb_pts=nb_pts)
        # self.assertLessEqual((err.max() / rad_max), 1e-2, 1)
        #
        # distance=np.arange(20,101,20,dtype='int')
        # err = sim_test.error_radiation_method_distance(RADIATION_METHOD_APPROX_FARFIELD, D=distance)
        # self.assertLessEqual(err.max()/rad_max, 1e-2, 1)

