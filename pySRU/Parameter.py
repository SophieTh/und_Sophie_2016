import numpy as np
import scipy.constants as codata
from abc import abstractmethod
from MagneticField import MagneticField

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
        By = (lambda z,y: self.magnetic_field(Z=z,Y=y,x=X,harmonic_number=harmonic_number,coordonnee='y'))
        Bz = (lambda z, y: self.magnetic_field(Z=z, Y=y,x=X, harmonic_number=harmonic_number, coordonnee='z'))
        Bx=(lambda z, y: self.magnetic_field(Z=z, Y=y, x=X,harmonic_number=harmonic_number, coordonnee='x'))
        #B = MagneticField(0.0, Y, Z, fct_null, By, fct_null)
        B = MagneticField(X, Y, Z, Bx, By, Bz)
        return B


    def gamma(self) :
        return self.E/0.511e6


    def Beta(self) :
        gamma=self.gamma()
        Beta=np.sqrt(1.0-1.0/gamma**2)
        return Beta

    @abstractmethod
    def Beta_et(self):
        return

    @abstractmethod
    def L(self):
        return

    @abstractmethod
    def Zo_symetry(self):
        return

    @abstractmethod
    def Zo_analitic(self):
        return

    @abstractmethod
    def Zmax_no_symetry(self):
        return

    @abstractmethod
    def omega1(self):
        return