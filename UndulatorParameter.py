import numpy as np
import scipy.constants as codata
from MagneticField import MagneticField


class UndulatorParameters(object):
    def __init__(self, K, E, lambda_u, L, I):
        self.K = K
        self.E = E
        self.lambda_u = lambda_u
        self.L = L
        self.I = I

    def copy(self):
        return UndulatorParameters(K=self.K,E=self.E,lambda_u=self.lambda_u,L=self.L,I=self.I)

    def creation_magnetic_field_plane_undulator(self,z):
        Bo = self.K / (93.4 * self.lambda_u)
        By = -Bo * np.sin((2.0 * np.pi / self.lambda_u) * z)
        # Hamming windowing
        windpar = 1.0 / (2.0 * np.floor(self.L/self.lambda_u))
        zmin = z.min()
        apo1 = zmin + windpar
        apo2 = z.max() - windpar
        wind = np.ones(len(z))
        for i in range(len(z)):
            if z[i] <= apo1:
                wind[i] *= 1.08 - (.54 + 0.46 * np.cos(np.pi * (z[i] - zmin) / windpar))
            if z[i] >= apo2:
                wind[i] *= 1.08 - (.54 - 0.46 * np.cos(np.pi * (z[i] - apo2) / windpar))
        By *= wind
        B=MagneticField(0.0,0.0,z,0.0,By,0.0)
        return B


    def omega1(self) :
        gamma=self.E/0.511e6
        return ((2.0 * gamma ** 2) / (1.0 + (self.K ** 2) / 2.0)) * ((2.0 * np.pi * codata.c) / self.lambda_u)


