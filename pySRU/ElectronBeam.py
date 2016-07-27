import numpy as np
import scipy.constants as codata


# container class
# parameter of the beam
class ElectronBeam(object):
    def __init__(self, E,I ):
        self.E=E
        self.I=I

    def copy(self):
        return ElectronBeam(E=self.E,I=self.I)

    def gamma(self):
        return self.E / 0.51099890221e6


    # the speed of the electron
    def Beta(self):
        gamma = self.gamma()
        Beta = np.sqrt(1.0 - 1.0 / gamma ** 2)
        return Beta

    def print_parameters(self):
        print(' Electron Beam')
        print('    Electron energy (eV): %.f' %self.E)
        print('    Intensity (A): %.f' % self.I)
        print('    Electron speed : %f' %self.Beta())

if __name__ == "__main__" :
    K = 1.87
    E = 1.3e9
    I = 1.0
    beam=ElectronBeam(E=E,I=I)

    print(beam.gamma())

    print(beam.Beta())