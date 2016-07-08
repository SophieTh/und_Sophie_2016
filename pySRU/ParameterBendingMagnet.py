import numpy as np
import scipy.constants as codata
from MagneticField import MagneticField
from ParameterPlaneUndulator import ParameterPlaneUndulator
from pySRU.Parameter import Parameter , PLANE_UNDULATOR,BENDING_MAGNET



class ParameterBendingMagnet(Parameter):
    def __init__(self, E, Bo, R, I, div):
        super(self.__class__, self).__init__(E=E, I=I, type_magnet=BENDING_MAGNET)
        self.R=R
        self.Bo = Bo
        self.div = div

    def copy(self):
        return ParameterBendingMagnet(R=self.R, E=self.E,Bo=self.Bo,div=self.div, I=self.I)

    def fct_magnetic_field2(self, z, y, x, harmonic_number, coordonnee='y'):
        lambda_inv = self.get_L() / 4.0
        lambda_h =lambda_inv
        ku = 2.0 * np.pi / self.get_lambda_u()

        if coordonnee == 'x':
            dB = 0.0
        else:
            if coordonnee == 'y':
                # codata.m_e * codata.c / codata.e= 0.00170450894933
                Bo = self.Bo* np.cosh(ku * y)
                # print(Bo)
                f_base = np.cos
            else:  # coordonnee == 'z' :
                Bo = -self.Bo* np.sinh(ku * y)
                f_base = np.sin

            # we enelarge the real effect of the magnetic field by 4 lambda_u
            # on each extremity of the undulator
            lambda_inv=self.L()/4.0

            L_magn_field = self.L() / 2.0 + 4.0*lambda_inv

            # we know considere that if -(L/2 + lambda_u/4) < Z < (L/2 + lambda_u/4)
            # the magnetic field if a classic cosinus (or sinus)
            L_cosinus_part = self.L() / 2.0

            if ((z < -L_magn_field) or (z > L_magn_field)):
                dB = 0.0

            else:
                if (z < -L_cosinus_part or z > L_cosinus_part):

                    # in this case, we work with a gaussian,
                    # so we shift the coordinate frame for use a central gaussian
                    if z < -L_cosinus_part:
                        sign = 1
                    else:  # z> L_cosinus_part
                        sign = -1

                    shift_z = z + sign * L_cosinus_part

                    p = 2.0 * np.pi ** 2 / (3.0 * lambda_h ** 2)
                    dB = (4.0*np.pi)*((2.0 * np.pi * Bo / lambda_h) * shift_z) * (
                    1.0 - 4.0 * np.pi ** 2 * shift_z ** 2 / (9.0 * lambda_h ** 2)) * np.exp(-p * shift_z ** 2)

                    # test du signe
                    z_test = sign * (-L_cosinus_part + lambda_h / 4.0)
                    test_trig = f_base(ku * z_test)
                    if (sign * test_trig < 0):
                        dB = -dB

                else:
                    # in this case we work in the cosinus part
                    dB = Bo

        return dB
    def fct_magnetic_field(self, z, y, x, harmonic_number, coordonnee='y'):
        if coordonnee == 'y':
                # codata.m_e * codata.c / codata.e= 0.00170450894933
                Bo = -self.Bo
        else:  # coordonnee == 'z' :
                Bo =0.0

            # we enelarge the real effect of the magnetic field by 4 lambda_u
            # on each extremity of the undulator
        L_magn_field = self.get_L()*1.5

            # we know considere that if -(L/2 + lambda_u/4) < Z < (L/2 + lambda_u/4)
            # the magnetic field if a classic cosinus (or sinus)
        L_cosinus_part = -self.Zo_analitic()

        if ((z < -L_magn_field) or (z > L_magn_field)):
            dB = 0.0

        else:
            if (z < -L_cosinus_part or z > L_cosinus_part):

                    # in this case, we work with a gaussian,
                    # so we shift the coordinate frame for use a central gaussian
                if z>0.0 :
                    z_shift=z-L_cosinus_part
                else :
                    z_shift = z + L_cosinus_part

                sigma=(L_magn_field-L_cosinus_part)/5.0
                #dB= Bo*(1.0-z_shift*self.get_L()**2/2.0)*np.exp(-0.5*(z_shift/sigma)**2)
                dB = Bo*np.exp(-0.5 * (z_shift / sigma) ** 2)
            else:
                    # in this case we work in the cosinus part
                dB = Bo
        return dB


    def get_L(self):
        return self.R*self.div

    def get_Bo(self):
        return self.Bo

    def get_K(self):
        return (codata.e*self.Bo*self.get_lambda_u())/(2.0*np.pi*codata.m_e*codata.c)


    def get_lambda_u(self):
        return self.get_L()*4.0

    def theta_max(self):
        return 5.0*self.get_K()/self.gamma()


    def Zmax_no_symetry(self):
        return self.get_L() + 10.0*self.get_L()/4.0


    def Zo_symetry(self):
        return -(self.get_L() / 2.0 + 5.0*self.get_L()/4.0)


    def Zo_analitic(self):
        return -self.get_L() / 2.0