import unittest

from pySRU.tests.MagneticFieldTest      import MagneticFieldTest
from pySRU.tests.MagneticStructureTest  import MagneticStructureTest
from pySRU.tests.MainTest               import MainTest
from pySRU.tests.RadiationFactorytest   import RadiationFactoryTest
from pySRU.tests.RadiationTest          import RadiationTest
from pySRU.tests.SimulationTest         import UndulatorSimulationTest
from pySRU.tests.SourceTest             import UndulatorParameterTest
from pySRU.tests.TrajectoryFactoryTest  import TrajectoryFactoryTest
from pySRU.tests.TrajectoryTest         import TrajectoryTest

def suite():
    """
    Gathers all the tests in a test suite.
    """
    suites = (
        # tests by Mark Glass.
        # unittest.makeSuite(MagneticFieldTest, 'test'),
        # unittest.makeSuite(MagneticStructureTest, 'test'),
        # unittest.makeSuite(MainTest, 'test'),
        # unittest.makeSuite(RadiationFactoryTest, 'test'),
        # unittest.makeSuite(RadiationTest, 'test'),
        # unittest.makeSuite(UndulatorSimulationTest, 'test'),
        unittest.makeSuite(UndulatorParameterTest, 'test'),
        # unittest.makeSuite(TrajectoryFactoryTest, 'test'),
        # unittest.makeSuite(TrajectoryTest, 'test'),

    )
    return unittest.TestSuite(suites)


if __name__ == "__main__":

    unittest.main(defaultTest="suite")
