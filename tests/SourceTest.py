import unittest
import numpy as np
from pySRU.MagneticStructureUndulatorPlane import MagneticStructureUndulatorPlane as Undulator
from pySRU.ElectronBeam import ElectronBeam
from pySRU.SourceUndulatorPlane import SourceUndulatorPlane


class UndulatorParameterTest(unittest.TestCase):

    def test_copy(self):
        K = 1.87
        lambda_u = 0.035
        L = 0.035 * 14
        E=1.3
        I=1.0
        undulator = Undulator(K=K, period_length=lambda_u, length=L)
        beam=ElectronBeam(Electron_energy=E, I_current=I)
        source=SourceUndulatorPlane(undulator=undulator, electron_beam=beam)

        source2=source.copy()
        source2.electron_beam.I_current=0.3
        self.assertEqual(source.I_current(),1.0)
        self.assertEqual(source2.I_current() ,0.3)

        self.assertAlmostEqual(source.harmonic_frequency(1)/1e17, 2.5346701615509917,10)
        self.assertAlmostEqual(source.Lorentz_factor(), 2544.0367765521196,10 )
        self.assertAlmostEqual(source.electron_speed(), 0.99999992274559524,15)
        self.assertEqual(source.magnetic_structure.period_number(), 14)
        self.assertAlmostEqual(source.choose_distance_automatic(2),49.000000000000,10 )



    #TODO completer ?