import unittest
import numpy as np
from pySRU.MagneticStructureUndulatorPlane import MagneticStructureUndulatorPlane as Undulator
from pySRU.ElectronBeam import ElectronBeam
from pySRU.Source import Source


class UndulatorParameterTest(unittest.TestCase):

    def test_copy(self):
        K = 1.87
        lambda_u = 0.035
        L = 0.035 * 14
        E=1.3e9
        I=1.0
        undulator = Undulator(K=K, lambda_u=lambda_u, L=L)
        beam=ElectronBeam(E=E,I=I)
        source=Source(magnetic_structure=undulator,electron_beam=beam)

        source2=source.copy()
        source2.electron_beam.I=0.3
        self.assertEqual(source.I(),1.0)
        self.assertEqual(source2.I() ,0.3)

        self.assertEqual(source.omega1(), 2.5346701615509917e+17)
        self.assertEqual(source.gamma(), 2544.0367765521196 )
        self.assertEqual(source.Beta(), 0.99999992274559524)
        self.assertEqual(source.Nb_period(), 14)
        self.assertEqual(source.D_min(2),24.500000000000004 )



