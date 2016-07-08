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

    def fct_magnetic_field(self, z, y, x, harmonic_number, coordonnee='y'):
        lambda_h = self.L()/ harmonic_number
        if coordonnee == 'y':
                # codata.m_e * codata.c / codata.e= 0.00170450894933
                Bo = -self.Bo
        else:  # coordonnee == 'z' :
                Bo =0.0

            # we enelarge the real effect of the magnetic field by 4 lambda_u
            # on each extremity of the undulator
        L_magn_field = self.L()

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

                sigma=(L_magn_field-L_cosinus_part)/4.0
                dB= Bo*(1.0-z_shift*self.L()**2/2.0)*np.exp(-0.5*(z_shift/sigma)**2)
            else:
                    # in this case we work in the cosinus part
                dB = Bo
        return dB

    # def fct_magnetic_field(self, z, y,x, harmonic_number, coordonnee='y'):
    #     if z <= self.L()/2.0 and z>= -self.L()/2.0 :
    #         if coordonnee=='y' :
    #             dB=-self.Bo
    #         else :
    #             dB=0.0
    #     else :
    #         dB=0.0
    #     return dB

    def omega1(self):
        gamma = self.gamma()
        #??????????????
        return 1.0

    def L(self):
        return self.R*self.div

    def K(self):
        return (codata.e*self.Bo*self.L())/(2.0*np.pi*codata.m_e*codata.c)

    def Beta_et(self) :
        Beta_et = 1.0 - (1.0 / (2.0 * self.gamma() ** 2)) * (1.0 + (self.K()** 2) / 2.0)
        return Beta_et

    def omega1(self) :
        gamma=self.gamma()
        omega1=((2.0 * gamma ** 2) / (1.0 + (self.K() ** 2) / 2.0)) * ((2.0 * np.pi * codata.c) / self.L())
        return omega1



    def Zmax_no_symetry(self):
        return self.L()*1.5


    def Zo_symetry(self):
        return -self.L()*1.25


    def Zo_analitic(self):
        return -self.L() / 2.0