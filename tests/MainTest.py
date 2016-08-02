import unittest
import numpy as np
from pySRU.MagneticStructureUndulatorPlane import MagneticStructureUndulatorPlane as Undulator
from pySRU.ElectronBeam import ElectronBeam
from pySRU.TrajectoryFactory import TRAJECTORY_METHOD_ANALYTIC,TRAJECTORY_METHOD_INTEGRATION,\
                                                        TRAJECTORY_METHOD_ODE
from pySRU.RadiationFactory import RADIATION_METHOD_APPROX_FARFIELD,RADIATION_METHOD_NEAR_FIELD ,\
                                RADIATION_METHOD_FARFIELD,RADIATION_METHOD_APPROX
from pySRU.Simulation import create_simulation

class MainTest(unittest.TestCase):

    def simul_undulator_near_to_farfield(self,magnetic_struc,electron_beam, method_traj,formule=1):

        sim_test = create_simulation(magnetic_structure=magnetic_struc,electron_beam=electron_beam,
                                     traj_method=method_traj,formule=formule)

        theta_max=sim_test.source.angle_wave_number(harmonic_number=1,wave_number=1)*0.3
        xy_max=theta_max*sim_test.radiation.distance
        X=np.linspace(0.0,xy_max,sim_test.radiation_fact.Nb_pts)
        Y=np.linspace(0.0,xy_max,sim_test.radiation_fact.Nb_pts)
        sim_test.change_XY_radiation(X=X,Y=Y)


        ref=sim_test.copy()
        rad_max=(sim_test.radiation.max())

        sim_test.change_radiation_method(RADIATION_METHOD_FARFIELD)
        rad_err=ref.radiation.error_max(sim_test.radiation)
        self.assertLessEqual(np.abs(sim_test.radiation.max()-rad_max)/rad_max,1e-2)
        self.assertLessEqual(rad_err / rad_max, 1e-2)

        sim_test.change_radiation_method(RADIATION_METHOD_APPROX)
        rad_err=ref.radiation.error_max(sim_test.radiation)
        self.assertLessEqual(np.abs(sim_test.radiation.max()-rad_max)/rad_max,1e-2)
        self.assertLessEqual(rad_err  / rad_max, 1e-2)

        sim_test.change_radiation_method(RADIATION_METHOD_NEAR_FIELD)
        rad_err=ref.radiation.error_max(sim_test.radiation)
        self.assertLessEqual(np.abs(sim_test.radiation.max()-rad_max)/rad_max,1e-2)
        self.assertLessEqual(rad_err / rad_max, 1e-2)

    def simul_undulator_traj_method(self, magnetic_struc,electron_beam, method_rad, formule=1):

        sim_test = create_simulation(magnetic_structure=magnetic_struc, electron_beam=electron_beam,
                                     rad_method=method_rad, formule=formule)

        ref = sim_test.copy()
        rad_max = (sim_test.radiation.max())

        sim_test.change_trajectory_method(TRAJECTORY_METHOD_ODE)
        rad_err = ref.radiation.error_max(sim_test.radiation)
        traj_err = ref.trajectory.error_max(sim_test.trajectory)
        self.assertLessEqual(np.abs(sim_test.radiation.max() - rad_max) / rad_max, 1e-2,4)
        self.assertLessEqual(traj_err[1] / ref.trajectory.x.max(), 5e-3)
        self.assertLessEqual(traj_err[3] / ref.trajectory.z.max(), 5e-3)
        self.assertLessEqual(traj_err[4] / ref.trajectory.v_x.max(), 5e-3)
        self.assertLessEqual(traj_err[6] / ref.trajectory.v_z.max(), 5e-3)
        self.assertLessEqual(traj_err[7] / ref.trajectory.a_x.max(), 5e-3)
        self.assertLessEqual(traj_err[9] / ref.trajectory.a_z.max(), 5e-3)
        self.assertLessEqual(np.abs(sim_test.radiation.max()-rad_max)/rad_max,1e-2)
        self.assertLessEqual(rad_err / rad_max, 1e-2)

        sim_test.change_trajectory_method(TRAJECTORY_METHOD_INTEGRATION)
        rad_err = ref.radiation.error_max(sim_test.radiation)
        traj_err=ref.trajectory.error_max(sim_test.trajectory)
        self.assertLessEqual(np.abs(sim_test.radiation.max() - rad_max) / rad_max, 1e-2)
        self.assertLessEqual(traj_err[1] / (ref.trajectory.x.max()), 5e-2)
        self.assertLessEqual(traj_err[3] / ref.trajectory.z.max(), 5e-2)
        self.assertLessEqual(traj_err[4] / ref.trajectory.v_x.max(), 5e-2)
        self.assertLessEqual(traj_err[6] / ref.trajectory.v_z.max(), 5e-2)
        self.assertLessEqual(traj_err[7] / ref.trajectory.a_x.max(), 5e-2)
        self.assertLessEqual(traj_err[9] / ref.trajectory.a_z.max(), 5e-2)
        self.assertLessEqual(np.abs(sim_test.radiation.max() - rad_max) / rad_max, 1e-2)
        self.assertLessEqual(rad_err / rad_max, 1e-2)



    def test_main(self):
        beam_test = ElectronBeam(Electron_energy=1.3, I_current=1.0)
        beam_ESRF = ElectronBeam(Electron_energy=6.0, I_current=0.2)
        und_test = Undulator(K=1.87, period_length=0.035, length=0.035 * 14)
        ESRF18 = Undulator(K=1.68, period_length=0.018, length=2.0)

        ##ESRFBM = BM(E=6.0e9, Bo=0.8, div=5e-3, R=25.0, I=0.2)

        # print(' undulator test ')
        # self.simul_undulator_near_to_farfield(magnetic_struc=und_test,electron_beam=beam_test
        #                                       ,method_traj=TRAJECTORY_METHOD_ANALYTIC,formule=1)
        # print("TRAJECTORY_METHOD_ANALYTIC ok")
        # self.simul_undulator_near_to_farfield(magnetic_struc=und_test, electron_beam=beam_test
        #                                       ,method_traj=TRAJECTORY_METHOD_INTEGRATION, formule=1)
        # print("TRAJECTORY_METHOD_INTEGRATION ok")
        # self.simul_undulator_near_to_farfield(magnetic_struc=und_test, electron_beam=beam_test
        #                         ,method_traj=TRAJECTORY_METHOD_ODE, formule=1)
        # print("TRAJECTORY_METHOD_ODE ok")
        #
        # self.simul_undulator_traj_method(magnetic_struc=und_test, electron_beam=beam_test
        #                                  , method_rad=RADIATION_METHOD_APPROX_FARFIELD)
        # print('APPROX FARFIELD ok')
        #
        # self.simul_undulator_traj_method(magnetic_struc=und_test,electron_beam=beam_test
        #                                  , method_rad=RADIATION_METHOD_FARFIELD)
        # print('FARFIELD ok')
        #
        # self.simul_undulator_traj_method(magnetic_struc=und_test, electron_beam=beam_test
        #                                  , method_rad=RADIATION_METHOD_APPROX)
        # print('APPROX ok')
        #
        # self.simul_undulator_traj_method(magnetic_struc=und_test, electron_beam=beam_test
        #                                  , method_rad=RADIATION_METHOD_NEAR_FIELD)
        # print('NEAR FIELD ok')
        print(' ')
        print('undulator ESRF18')
        # self.simul_undulator_near_to_farfield(magnetic_struc=ESRF18, electron_beam=beam_ESRF
        #                                       ,method_traj=TRAJECTORY_METHOD_ANALYTIC, formule=1)
        # print("TRAJECTORY_METHOD_ANALYTIC ok")

        # TODO ne marche pas
        # self.simul_undulator_near_to_farfield(magnetic_struc=ESRF18, electron_beam=beam_ESRF
        #                                     ,method_traj=TRAJECTORY_METHOD_INTEGRATION, formule=1)
        # print("TRAJECTORY_METHOD_INTEGRATION ok")

        self.simul_undulator_near_to_farfield(magnetic_struc=ESRF18, electron_beam=beam_ESRF
                                              , method_traj=TRAJECTORY_METHOD_ODE, formule=1)
        print("TRAJECTORY_METHOD_ODE ok")

        self.simul_undulator_traj_method(magnetic_struc=ESRF18, electron_beam=beam_ESRF
                                              , method_rad=RADIATION_METHOD_APPROX_FARFIELD)
        print('APPROX FARFIELD ok')

        self.simul_undulator_traj_method(magnetic_struc=ESRF18, electron_beam=beam_ESRF
                                         , method_rad=RADIATION_METHOD_FARFIELD)
        print('FARFIELD ok')

        self.simul_undulator_traj_method(magnetic_struc=ESRF18, electron_beam=beam_ESRF
                                         , method_rad=RADIATION_METHOD_APPROX)
        print('APPROX ok')

        self.simul_undulator_traj_method(magnetic_struc=ESRF18, electron_beam=beam_ESRF
                                         , method_rad=RADIATION_METHOD_NEAR_FIELD)
        print('NEAR FIELD ok')