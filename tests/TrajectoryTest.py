import unittest
import numpy as np
from pySRU.Trajectory import Trajectory


class TrajectoryTest(unittest.TestCase):
    def test_error_nb_pts_copy(self):
        t = np.arange(0.0, 6.0*np.pi, 0.005*np.pi)
        x = np.sin(t)
        y = 0.0 * t
        z = np.cos(t)
        v_x = np.cos(t)
        v_y = 0.0 * t
        v_z = -np.sin(t)
        a_x = -np.sin(t)
        a_y = 0.0 * t
        a_z = -np.cos(t)
        Traj = Trajectory(t, x, y, z, v_x, v_y, v_z, a_x, a_y, a_z)

        Traj2 = Traj.copy()

        self.assertEqual(Traj.nb_points(), 1200)


        self.assertTrue(all(Traj.t==Traj2.t))

        error=Traj.error_max(Traj2)
        # Traj.draw_2_trajectory(Traj2)
        error[0] -= len(Traj.t)
        self.assertEqual(error.max(), 0.0)

        Traj2.multiply_by(3)
        error = Traj.error_max(Traj2)
        error[0] -= len(Traj.t)
        self.assertEqual(error.max(), 2.0)