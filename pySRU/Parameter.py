import numpy as np
import scipy.constants as codata
from abc import abstractmethod
from pySRU.MagneticField import MagneticField

PLANE_UNDULATOR=0
BENDING_MAGNET=1

class Parameter(object):
    def __init__(self, E, I, type_magnet):
        self.type_magnet=type_magnet
        self.E = E
        self.I = I


    @abstractmethod
    def copy(self):
        return

    @abstractmethod
    def fct_magnetic_field(self, z, y,x, harmonic_number, coordonnee='y'):
        return

    def magnetic_field(self, Z, Y,x, harmonic_number=1, coordonnee='y'):
        if (type(Z) == np.ndarray):
            B = np.zeros_like(Z)
        else:
            if (type(Y) == np.ndarray):
                B = np.zeros_like(Y)
            else:
                B = self.fct_magnetic_field(z=Z, y=Y,x=x,
                                            harmonic_number=harmonic_number, coordonnee=coordonnee)

        if (type(Z) == np.ndarray):
            # B=np.zeros_like(Z)
            if (type(Y) == np.ndarray):
                for i, Zi in enumerate(Z):
                    B[i] = self.fct_magnetic_field(z=Zi, y=Y[i],x=x,
                                                   harmonic_number=harmonic_number, coordonnee=coordonnee)
            else:
                for i, Zi in enumerate(Z):
                    B[i] = self.fct_magnetic_field(z=Zi, y=Y,x=x,
                                                   harmonic_number=harmonic_number, coordonnee=coordonnee)
        else:

            if (type(Y) == np.ndarray):

                for i, Yi in enumerate(Y):
                    # B = np.zeros_like(Y)
                    B[i] = self.fct_magnetic_field(z=Z, y=Yi,x=x, harmonic_number=harmonic_number, coordonnee=coordonnee)

        return B

    def create_magnetic_field(self,Z,Y,X,harmonic_number):
        By = (lambda z,y,x: self.magnetic_field(Z=z,Y=y,x=X,harmonic_number=harmonic_number,coordonnee='y'))
        Bz = (lambda z, y,x: self.magnetic_field(Z=z, Y=y,x=X, harmonic_number=harmonic_number, coordonnee='z'))
        Bx=(lambda z, y,x: self.magnetic_field(Z=z, Y=y, x=X,harmonic_number=harmonic_number, coordonnee='x'))
        #B = MagneticField(0.0, Y, Z, fct_null, By, fct_null)
        B = MagneticField(X, Y, Z, Bx, By, Bz)
        return B

    @abstractmethod
    def get_L(self):
        return

    @abstractmethod
    def get_K(self):
        return

    @abstractmethod
    def get_lambda_u(self):
        return

    @abstractmethod
    def get_Bo(self):
        return

    @abstractmethod
    def theta_max(self):
        return self.get_K()/self.gamma()

    def gamma(self) :
        return self.E/0.511e6

    def Beta(self) :
        gamma=self.gamma()
        Beta=np.sqrt(1.0-1.0/gamma**2)
        return Beta

    def omega1(self):
        gamma = self.gamma()
        first_harm = ((2.0 * gamma ** 2) / (1.0 + (self.get_K() ** 2) / 2.0)) * ((2.0 * np.pi * codata.c) / self.get_lambda_u())
        return first_harm

    def Beta_et(self):
        Beta_et = 1.0 - (1.0 / (2.0 * self.gamma() ** 2)) * (1.0 + (self.get_K()** 2) / 2.0)
        return Beta_et

    def Nb_period(self):
        return np.floor(self.get_L() / self.get_lambda_u())


    def D_max(self, alpha):
        lim = self.get_L() / 2.0
        return lim * 10 ** alpha

    @abstractmethod
    def Zmax_no_symetry(self):
        return

    @abstractmethod
    def Zo_symetry(self):
        return

    @abstractmethod
    def Zo_analitic(self):
        return

    # n is the harmonic number
    # l un the wave number
    def theta(self, n, l):
        if n == 0:
            raise Exception('n must be != 0')
        return np.sqrt((l / n) * (1.0 + self.get_K() ** 2 / 2.0)) * (1.0 / self.gamma())
