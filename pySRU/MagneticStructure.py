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
    def Zmax_no_symetry(self):
        return

    @abstractmethod
    def Zo_symetry(self):
        return

    @abstractmethod
    def Zo_analitic(self):
        return

    @abstractmethod
    def print_parameters(self):
        pass

    def Nb_period(self):
        return np.floor(self.get_L() / self.get_lambda_u())

    def D_min(self, alpha):
        lim = self.get_L() / 2.0
        return lim * 10 ** alpha

    @abstractmethod
    def fct_magnetic_field(self, z, y,x, harmonic_number, coordonnee='y'):
        return

    @abstractmethod
    def flux_on_axis(self,n,electron_beam):
        return


    # function wich receive X,Y,Z in array or float
    # and give a array
    # if two of them are array, they must have the same lenght
    #TODO la faire plus jolie
    def magnetic_field(self, Z, Y,X, harmonic_number=1, coordonnee='y'):
        if (type(Z) == np.ndarray):
            B = np.zeros_like(Z)
        else:
            if (type(Y) == np.ndarray):
                B = np.zeros_like(Y)
            else:
                if (type(X) == np.ndarray):
                    B = np.zeros_like(X)
                else :
                    B = self.fct_magnetic_field(z=Z, y=Y,x=X,
                                            harmonic_number=harmonic_number, coordonnee=coordonnee)

        if (type(Z) == np.ndarray):
            if (type(Y) == np.ndarray):
                if (len(Y) != len(Z)) :
                    raise Exception (' Y and Z must have the same lenght')
                if (type(X) == np.ndarray):
                    if (len(X) != len(Z)):
                        raise Exception(' X and Z must have the same lenght')
                for i, Zi in enumerate(Z):
                    B[i] = self.fct_magnetic_field(z=Zi, y=Y[i], x=X[i],
                                                   harmonic_number=harmonic_number, coordonnee=coordonnee)
                else :
                    for i, Zi in enumerate(Z):
                        B[i] = self.fct_magnetic_field(z=Zi, y=Y[i],x=X,
                                                   harmonic_number=harmonic_number, coordonnee=coordonnee)
            else:
                if (type(X) == np.ndarray):
                    if (len(X) != len(Z)):
                        raise Exception(' X and Z must have the same lenght')
                    for i, Zi in enumerate(Z):
                        B[i] = self.fct_magnetic_field(z=Zi, y=Y,x=X[i],
                                                   harmonic_number=harmonic_number, coordonnee=coordonnee)
                else :
                    for i, Zi in enumerate(Z):
                        B[i] = self.fct_magnetic_field(z=Zi, y=Y, x=X,
                                           harmonic_number=harmonic_number, coordonnee=coordonnee)
        else:
            if (type(Y) == np.ndarray):
                if (type(X) == np.ndarray):
                    if (len(X) != len(Y)):
                        raise Exception(' X and Z must have the same lenght')
                for i, Yi in enumerate(Y):
                    B[i] = self.fct_magnetic_field(z=Z, y=Yi,x=X[i], harmonic_number=harmonic_number, coordonnee=coordonnee)
            else :
                if (type(X) == np.ndarray):
                    for i, Xi in enumerate(X):
                        B[i] = self.fct_magnetic_field(z=Z, y=Y,x=Xi, harmonic_number=harmonic_number, coordonnee=coordonnee)
                else :
                    B=self.fct_magnetic_field(z=Z, y=Y,x=X, harmonic_number=harmonic_number, coordonnee=coordonnee)

        return B

    # a Magnetic structur create a magnetic field
    # the object MagneticField is create like :
    # Bx , By , Bz are function R**(3) -> R
    # this function depend of the magnet type (BendingMagnet or undulator ..)
    # X, Y can be array or number, they describe the area where we want to work
    # Z must be an array , it will be use later
    # they have not necessary the same len
    # peux etre a changer et ne metrre que les fonctions ... oui
    def create_magnetic_field(self,harmonic_number=1):
        By = (lambda z,y,x: self.magnetic_field(Z=z,Y=y,X=x,harmonic_number=harmonic_number,coordonnee='y'))
        Bz = (lambda z, y,x: self.magnetic_field(Z=z, Y=y,X=x, harmonic_number=harmonic_number, coordonnee='z'))
        Bx = (lambda z, y,x: self.magnetic_field(Z=z, Y=y,X=x,harmonic_number=harmonic_number, coordonnee='x'))
        B = MagneticField(Bx, By, Bz)
        return B


def Exemple1_undulator():
    from MagneticStructureUndulatorPlane import MagneticStructureUndulatorPlane as Undulator

    undulator_test = Undulator(K=1.87, lambda_u=0.035, L=0.035 * 14)

    print(' intensity maximum magnetic fiel (T) : ')
    print(undulator_test.get_Bo())
    print(' K : ')
    print(undulator_test.get_K())
    magnetic_field_test = undulator_test.create_magnetic_field()

    Z = np.linspace(undulator_test.Zo_symetry(), -undulator_test.Zo_symetry(), 401)

    magnetic_field_test.plot_z(X=0.0, Y=0.0, Z=Z)
    magnetic_field_test.plot_z(X=0.1, Y=0.0, Z=Z)
    magnetic_field_test.plot_z(X=0.0, Y=0.2, Z=Z)

    X = np.linspace(-0.01,0.01,401)
    Y = np.linspace(-0.01, 0.01, 401)
    magnetic_field_test.plot_z(X=X, Y=Y, Z=Z)
    # magnetic_field_test.plot_y(X=X, Y=Y, Z=Z)
    # magnetic_field_test.plot_x(X=X, Y=Y, Z=Z)

def Exemple2_BM():
    from MagneticStructureBendingMagnet import MagneticStructureBendingMagnet as BM

    undulator_test = BM(Bo=0.8,div=5e-3,R=25.0)

    print(' magnet lenght (m): ')
    print(undulator_test.get_L())
    print(' magnet lenght (m): ')
    print(undulator_test.get_L())

    Z = np.linspace(undulator_test.Zo_symetry(), -undulator_test.Zo_symetry(), 101)
    magnetic_field_test = undulator_test.create_magnetic_field()
    magnetic_field_test.plot_z(X=0.0, Y=0.0, Z=Z)
    magnetic_field_test.plot_z(X=0.1, Y=0.0, Z=Z)
    magnetic_field_test.plot_z(X=0.0, Y=0.2, Z=Z)

    X = np.linspace(-0.01, 0.01,401)
    Y = np.linspace(-0.01, 0.01, 401)
    magnetic_field_test.plot_z(X=X, Y=Y, Z=Z)
    magnetic_field_test.plot_z(X=X, Y=Y, Z=Z)
    # magnetic_field_test.plot_y(X=X, Y=Y, Z=Z)
    # magnetic_field_test.plot_x(X=X, Y=Y, Z=Z)



if __name__ == "__main__" :

    Exemple1_undulator()
    Exemple2_BM()
