import numpy as np
import scipy.constants as codata
from abc import abstractmethod
from pySRU.MagneticField import MagneticField
from pySRU.MagneticStructureUndulatorPlane import MagneticStructureUndulatorPlane as Undulator
from pySRU.MagneticStructureBendingMagnet import MagneticStructure
from pySRU.ElectronBeam import ElectronBeam

PLANE_UNDULATOR=0
BENDING_MAGNET=1


# magnetic field est il interressant ici ???
class Source(object):
    def __init__(self,magnetic_structure,electron_beam,magnetic_field=None):
        self.magnetic_structure=magnetic_structure
        self.electron_beam=electron_beam
        if magnetic_field==None :
            self.magnetic_field=self.magnetic_structure.create_magnetic_field()
        else :
            self.magnetic_field=magnetic_field


    def copy(self):
        return Source(magnetic_structure=self.magnetic_structure.copy(),
                      electron_beam=self.electron_beam.copy(),
                      magnetic_field=self.magnetic_field.copy())


    # gueter
    def I(self):
        return self.electron_beam.I

    def E(self):
        return self.electron_beam.E

    def L(self):
        return self.magnetic_structure.get_L()

    def K(self):
        return self.magnetic_structure.get_K()

    def lambda_u(self):
        return self.magnetic_structure.get_lambda_u()

    def Bo(self):
        return self.magnetic_structure.get_Bo()

    def Nb_period(self):
        return self.magnetic_structure.Nb_period()

    def magnet_type(self):
        return self.magnetic_structure.magnet_type

    def D_min(self, alpha):
        return self.magnetic_structure.D_min(alpha)

    def Zmax_no_symetry(self):
        return self.magnetic_structure.Zmax_no_symetry()

    def Zo_symetry(self):
        return self.magnetic_structure.Zo_symetry()

    def Zo_analitic(self):
        return self.magnetic_structure.Zo_analitic()

    def gamma(self):
        return self.electron_beam.gamma()

    def Beta(self):
        gamma = self.electron_beam.gamma()
        Beta = np.sqrt(1.0 - 1.0 / gamma ** 2)
        return Beta

    # source parameter
    def print_parameters(self) :
        if self.magnet_type()==PLANE_UNDULATOR :
            print(' Plane Undulator')
        else :
            print('Bending Magnet')

        self.magnetic_structure.print_parameters()

        self.electron_beam.print_parameters()

    # the arbitrary maximum angle use for calculate radiation
    def theta_max(self):
        return 2.0*self.K()/self.gamma()
        # if self.magnetic_structure.magnet_type==PLANE_UNDULATOR :
        #     theta_max= self.theta(1,1)
        # else :
        #     theta_max= self.magnetic_structure.get_K()/self.electron_beam.gamma()
        # return theta_max

    def n_min(self, alpha):
        racine = (np.pi * self.Nb_period() * 10 ** (alpha)) / (45.)
        n_mini = (2.0 * np.pi * self.Nb_period()) * np.exp(0.25 * np.log(racine))
        return n_mini

    def omega1(self):
        gamma = self.electron_beam.gamma()
        first_harm = ((2.0 * gamma ** 2) / (1.0 + (self.K() ** 2) / 2.0)) * (
            (2.0 * np.pi * codata.c) / self.lambda_u())
        return first_harm

    def Beta_et(self):
        Beta_et = 1.0 - (1.0 / (2.0 * self.gamma() ** 2)) * (1.0 + (self.K()** 2) / 2.0)
        return Beta_et

    def harmonic_number(self,omega):
        return omega/self.omega1()

    # for the harmonic n
    # the angle of the wave number l
    # useful above all for undulator
    def theta(self, n, l):
        if n == 0:
            n=1
        return np.sqrt((l / n) * (1.0 + self.K() ** 2 / 2.0)) * (1.0 / self.gamma())

    def flux_on_axis_theoric(self,n):
        spectre=self.magnetic_structure.flux_on_axis(n=n,electron_beam=self.electron_beam)
        return spectre



def Exemple1_undulator():

    undulator_test = Undulator(K=1.87, lambda_u=0.035, L=0.035 * 14)
    electron_beam_test=ElectronBeam(E=1.3e9,I=1.0)
    magnetic_field_test = undulator_test.create_magnetic_field()
    source=Source(magnetic_structure=undulator_test,
                  electron_beam=electron_beam_test,
                  magnetic_field=magnetic_field_test)


    print( " speed average in Z direction (m/s)/c")
    print(source.Beta_et())
    print(" first frequency")
    print(source.omega1())
    print(" angle of the first wave for the first number (radian) ")
    print(source.theta(1,1))


if __name__ == "__main__" :


    Exemple1_undulator()