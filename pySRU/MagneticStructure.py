import numpy as np
import scipy.constants as codata
from abc import abstractmethod
from pySRU.MagneticField import MagneticField

PLANE_UNDULATOR=0
BENDING_MAGNET=1

class MagneticStructure(object):
    def __init__(self, magnet_type ):
        self.magnet_type=magnet_type

    @abstractmethod
    def copy(self):
        return


    @abstractmethod
    def print_parameters(self):
        pass


    # magnetic field
    @abstractmethod
    def fct_magnetic_field(self, z, y, x, harmonic_number, coordonnee='y'):
        return 0.0


    # TODO la faire plus jolie ?
    def calculate_magnetic_field(self, Z, Y, X, harmonic_number=1, coordonnee='y'):
        if (type(Z) == np.ndarray):
            B = np.zeros_like(Z)
        else:
            if (type(Y) == np.ndarray):
                B = np.zeros_like(Y)
            else:
                if (type(X) == np.ndarray):
                    B = np.zeros_like(X)
                else:
                    B = self.fct_magnetic_field(z=Z, y=Y, x=X,
                                                harmonic_number=harmonic_number, coordonnee=coordonnee)

        if (type(Z) == np.ndarray):
            if (type(Y) == np.ndarray):
                if (len(Y) != len(Z)):
                    raise Exception(' Y and Z must have the same lenght')
                if (type(X) == np.ndarray):
                    if (len(X) != len(Z)):
                        raise Exception(' X and Z must have the same lenght')
                    for i, Zi in enumerate(Z):
                        B[i] = self.fct_magnetic_field(z=Zi, y=Y[i], x=X[i],
                                                   harmonic_number=harmonic_number, coordonnee=coordonnee)
                else:
                    for i, Zi in enumerate(Z):
                        B[i] = self.fct_magnetic_field(z=Zi, y=Y[i], x=X,
                                                       harmonic_number=harmonic_number, coordonnee=coordonnee)
            else:
                if (type(X) == np.ndarray):
                    if (len(X) != len(Z)):
                        raise Exception(' X and Z must have the same lenght')
                    for i, Zi in enumerate(Z):
                        B[i] = self.fct_magnetic_field(z=Zi, y=Y, x=X[i],
                                                       harmonic_number=harmonic_number, coordonnee=coordonnee)
                else:
                    for i, Zi in enumerate(Z):
                        B[i] = self.fct_magnetic_field(z=Zi, y=Y, x=X,
                                                       harmonic_number=harmonic_number, coordonnee=coordonnee)
        else:
            if (type(Y) == np.ndarray):
                if (type(X) == np.ndarray):
                    if (len(X) != len(Y)):
                        raise Exception(' X and Z must have the same lenght')
                for i, Yi in enumerate(Y):
                    B[i] = self.fct_magnetic_field(z=Z, y=Yi, x=X[i], harmonic_number=harmonic_number,
                                                   coordonnee=coordonnee)
            else:
                if (type(X) == np.ndarray):
                    for i, Xi in enumerate(X):
                        B[i] = self.fct_magnetic_field(z=Z, y=Y, x=Xi, harmonic_number=harmonic_number,
                                                       coordonnee=coordonnee)
                else:
                    B = self.fct_magnetic_field(z=Z, y=Y, x=X, harmonic_number=harmonic_number, coordonnee=coordonnee)
        return B


    # a Magnetic structur create a magnetic field
    # the object MagneticField is create like :
    # Bx , By , Bz are function R**(3) -> R
    # this function depend of the magnet type (BendingMagnet or undulator ..)
    # X, Y can be array or number, they describe the area where we want to work
    # Z must be an array , it will be use later
    # they have not necessary the same len
    # peux etre a changer et ne metrre que les fonctions ... oui
    def create_magnetic_field(self, harmonic_number=1):
        By = (lambda z, y, x: self.calculate_magnetic_field(Z=z, Y=y, X=x, harmonic_number=harmonic_number, coordonnee='y'))
        Bz = (lambda z, y, x: self.calculate_magnetic_field(Z=z, Y=y, X=x, harmonic_number=harmonic_number, coordonnee='z'))
        Bx = (lambda z, y, x: self.calculate_magnetic_field(Z=z, Y=y, X=x, harmonic_number=harmonic_number, coordonnee='x'))
        B = MagneticField(Bx, By, Bz)
        return B



if __name__ == "__main__" :
    pass
