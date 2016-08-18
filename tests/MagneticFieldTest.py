import unittest
import numpy as np
import scipy.integrate as integrate
import scipy.constants as codata
import matplotlib.pyplot as plt
from pySRU.MagneticStructureUndulatorPlane import MagneticStructureUndulatorPlane as Undulator
from pySRU.ElectronBeam import ElectronBeam
from pySRU.SourceBendingmagnet import BENDING_MAGNET as BM
from pySRU.TrajectoryFactory import TrajectoryFactory,TRAJECTORY_METHOD_ANALYTIC, TRAJECTORY_METHOD_ODE
from pySRU.RadiationFactory import RadiationFactory , RADIATION_METHOD_APPROX_FARFIELD,RADIATION_METHOD_NEAR_FIELD
from pySRU.Simulation import Simulation,create_simulation

class MagneticFieldTest(unittest.TestCase):



    def create_magn_field_undulator_test(self, magnetic_structure, electron_beam,method_traj):


        sim_test = create_simulation(magnetic_structure=magnetic_structure,electron_beam=electron_beam,
                                     traj_method=method_traj,rad_method=RADIATION_METHOD_APPROX_FARFIELD)
        print('create')
        source_test=sim_test.source
        Zo=sim_test.trajectory_fact.initial_condition[5]
        Bo=source_test.magnetic_field_strength()
        lambda_u=source_test.magnetic_structure.period_length

        self.assertTrue(Zo<0.)
        Zo_analitic=source_test.magnetic_structure.length*0.5
        Z_test = np.linspace(-Zo_analitic, Zo_analitic, sim_test.trajectory_fact.Nb_pts)
        diff_mf = np.abs(source_test.magnetic_field.By(z=Z_test, y=0.0, x=0.0) -
                         Bo * np.cos((2.0 * np.pi / lambda_u) * Z_test))
        self.assertTrue(all(diff_mf < np.abs(Bo) * 1e-3))

        print('integration 1')
        int1 = integrate.quad((lambda z: source_test.magnetic_field.By(z=z, y=0.0, x=0.0)),Zo,-Zo
                              ,limit=int(source_test.choose_nb_pts_trajectory(2)))
        self.assertAlmostEqual(int1[0], 0.0, 2)
        print('integration 2')#TODO marche pas donne tjr zeros meme qd c'est faux
        int2 = integrate.quad((lambda z: z*source_test.magnetic_field.By(z=z, y=0.0, x=0.0)),Zo,-Zo
                              ,limit=int(source_test.choose_nb_pts_trajectory(2)))
        print(int2[0])
        self.assertAlmostEqual(int2[0], 0.0)


    # doen't work ...
    # def create_magn_field_test_BM(self, BM, method_traj, formule=1):
    #     distance = 250.0
    #     Xmax = distance * 1e-3
    #     Ymax = distance * 1e-3
    #
    #     traj_test = TrajectoryFactory(Nb_pts=2000, method=method_traj)
    #     rad_test = RadiationFactory(method=RADIATION_METHOD_APPROX_FARFIELD, omega=BM.omega1(), Nb_pts=100,
    #                                 formula=formule)
    #     sim_test = create_simulation(parameter=BM, trajectory_fact=traj_test, radiation_fact=rad_test,
    #                                  distance=distance, X_max=Xmax, Y_max=Ymax)
    #     if method_traj == TRAJECTORY_METHOD_ANALYTIC:
    #         Zo = BM.Zo_analitic()
    #     else:
    #         Zo = BM.Zo_symetry()
    #     self.assertEqual(sim_test.magnetic_field.z[0], Zo)
    #     self.assertEqual(sim_test.magnetic_field.z[-1], -Zo)
    #     self.assertEqual(len(sim_test.magnetic_field.z), sim_test.trajectory_fact.Nb_pts)
    #     Z_test = np.linspace(BM.Zo_analitic(), -BM.Zo_analitic(), sim_test.trajectory_fact.Nb_pts)
    #     diff_mf = np.abs(sim_test.magnetic_field.By(z=Z_test, y=0.0, x=0.0) -
    #                      BM.Bo * np.cos((2.0 * np.pi / BM.lambda_u) * Z_test))
    #     self.assertTrue(all(diff_mf < np.abs(BM.Bo) * 1e-3))
    #     int1 = integrate.quad((lambda z: sim_test.magnetic_field.By(z=z, y=0.0, x=0.0)),
    #                           sim_test.magnetic_field.z[0],
    #                           sim_test.magnetic_field.z[-1])
    #     self.assertAlmostEqual(int1[0], 0.0, 5)
    #     int2 = integrate.quad((lambda z1: integrate.quad((lambda z: sim_test.magnetic_field.By(z=z, y=0.0, x=0.0)),
    #                                                      sim_test.magnetic_field.z[0], z1)[0])
    #                           , sim_test.magnetic_field.z[0], sim_test.magnetic_field.z[-1])
    #     self.assertAlmostEqual(int2[0], 0.0, 5)

    def test_magn_field(self):
        beam_test = ElectronBeam(Electron_energy=1.3, I_current=1.0)
        beam_ESRF = ElectronBeam(Electron_energy=6.0, I_current=0.2)
        und_test = Undulator(K=1.87, period_length=0.035, length=0.035 * 14)
        ESRF18 = Undulator(K=1.68, period_length=0.018, length=2.)

        self.create_magn_field_undulator_test(magnetic_structure=und_test,electron_beam=beam_test,
                                     method_traj=TRAJECTORY_METHOD_ODE)
        print("und_test, ok")

        self.create_magn_field_undulator_test(magnetic_structure=ESRF18,electron_beam=beam_ESRF,
                                    method_traj=TRAJECTORY_METHOD_ODE)
        print("esrf18, ok")

        # self.create_magn_field_test_BM(BM=ESRFBM, method_traj=TRAJECTORY_METHOD_ODE, formule=1)
        # print("esrf BM ok")