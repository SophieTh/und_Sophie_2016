import numpy as np
import scipy.constants as codata
from abc import abstractmethod
from pySRU.MagneticField import MagneticField
from pySRU.ElectronBeam import ElectronBeam
from pySRU.MagneticStructure import MagneticStructure,PLANE_UNDULATOR,BENDING_MAGNET




# magnetic field est il interressant ici ???
class Source(object):
    def __init__(self,electron_beam,magnetic_field,magnetic_structure):
        self.electron_beam=electron_beam
        self.magnetic_structure=magnetic_structure
        if magnetic_field==None :
            self.magnetic_field=self.magnetic_structure.create_magnetic_field()
        else :
            self.magnetic_field=magnetic_field


    @abstractmethod
    def copy(self):
        return

    # gueter
    def I_current(self):
        return self.electron_beam.I_current

    def Electron_energy(self):
        return self.electron_beam.Electron_energy

    # utiles
    def Lorentz_factor(self):
        return self.electron_beam.Lorentz_factor()

    #utile ?
    def electron_speed(self):#BETA
        gamma = self.electron_beam.Lorentz_factor()
        Beta = np.sqrt(1.0 - 1.0 / gamma ** 2)
        return Beta

    def magnet_type(self):
        return self.magnetic_structure.magnet_type


    @abstractmethod
    def magnetic_field_strength(self):
        return

    # constante de construction

    @abstractmethod
    def choose_distance_automatic(self, alpha):
        return

    @abstractmethod
    def choose_nb_pts_trajectory(self, alpha):
        return

    @abstractmethod
    def choose_initial_contidion_automatic(self):
        return

    @abstractmethod
    def choose_photon_frequency(self):
        return

    @abstractmethod
    def angle_deflection_max(self):
        return

    @abstractmethod
    def angle_deflection_central_cone(self):
        return

    @abstractmethod
    def analytical_times_vector(self,Nb_pts):
        return

    @abstractmethod
    def construct_times_vector(self,initial_contition,Nb_pts):
        return


    # source parameter


    def print_parameters(self) :
        self.electron_beam.print_parameters()
        self.magnetic_structure.print_parameters()

    # theory
    @abstractmethod
    def flux_on_axis_theoric(self,omega):
        return








if __name__ == "__main__" :
    from SourceUndulatorPlane import SourceUndulatorPlane
    from MagneticStructureUndulatorPlane import MagneticStructureUndulatorPlane
    electron_beam_test=ElectronBeam(Electron_energy=1.3e9, I_current=1.0)
    undulator_test=MagneticStructureUndulatorPlane(K=1.87,period_length=0.035,length=0.49)
    source_test = SourceUndulatorPlane(undulator=undulator_test, electron_beam=electron_beam_test,
                                       magnetic_field=None)
    print(type(source_test.magnetic_field))
    print(type(source_test.magnetic_field.Bx))
    print(type(source_test.magnetic_field.By))
    print(source_test.magnetic_field.Bx(0.0, 0.0, 0.0))

    #Exemple1_undulator()