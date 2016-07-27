import unittest
import numpy as np
import scipy.constants as codata
from pySRU.MagneticStructureUndulatorPlane import MagneticStructureUndulatorPlane as Undulator
from pySRU.ElectronBeam import ElectronBeam
from pySRU.Source import Source
from pySRU.TrajectoryFactory import TrajectoryFactory,TRAJECTORY_METHOD_ANALYTIC,TRAJECTORY_METHOD_INTEGRATION,\
                                                        TRAJECTORY_METHOD_ODE

class TrajectoryFactoryTest(unittest.TestCase):

    #TODO a completer
    def test_create_trajectory(self):

        # test le trajectoire ana lytic (des valeurs special
        # le Beta doit est constant pour certaine methode
        # les produi scalaire de l'acc avec la vitesse est ...
        # difference max entre deux trajectoire

        undulator_test = Undulator(K=1.87, lambda_u=0.035, L=0.035 * 14)
        electron_beam_test = ElectronBeam(E=1.3e9, I=1.0)
        source_test=Source(magnetic_structure=undulator_test,electron_beam=electron_beam_test)

        fact_test = TrajectoryFactory(Nb_pts=201,method=TRAJECTORY_METHOD_INTEGRATION)
        traj_test = fact_test.create_from_source(source=source_test)

        self.assertFalse(fact_test.initial_condition==None)

        self.assertTrue(all(fact_test.initial_condition==np.array([0.0,0.0,source_test.Beta()*codata.c,
                                                     0.0,0.0,source_test.Zo_symetry()])))
        self.assertTrue(fact_test.method==TRAJECTORY_METHOD_INTEGRATION)

        scalar_product=traj_test.v_x*traj_test.a_x + traj_test.a_z*traj_test.v_z

        self.assertAlmostEqual(np.abs(scalar_product).max(),0.0,3)







