import numpy as np
import scipy.constants as codata


# container class
# parameter of the beam
class ElectronBeam(object):
    def __init__(self, Electron_energy, I_current):
        self.Electron_energy=Electron_energy
        self.I_current=I_current

    def copy(self):
        return ElectronBeam(Electron_energy=self.Electron_energy, I_current=self.I_current)

    def Lorentz_factor(self):
        return (self.Electron_energy / 0.51099890221)*1e3


    def electron_speed(self):
        gamma = self.Lorentz_factor()
        Beta = np.sqrt(1.0 - 1.0 / gamma ** 2)
        return Beta

    def print_parameters(self):
        print(' Electron Beam')
        print('    Electron energy (GeV): %.3f' %self.Electron_energy)
        print('    Intensity (A): %.3f' % self.I_current)
        print('    Electron speed : %.10f' %self.electron_speed())
        print('    Lorentz factor : %.7f' % self.Lorentz_factor())

if __name__ == "__main__" :
    K = 1.87
    E = 1.3
    I = 1.0
    beam=ElectronBeam(Electron_energy=E, I_current=I)

    beam.print_parameters()

    print('Lorentz factor: %f'%beam.Lorentz_factor())