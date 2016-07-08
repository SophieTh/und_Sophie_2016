import numpy as np
import scipy.constants as codata
from MagneticField import MagneticField
from pySRU.Parameter import Parameter , PLANE_UNDULATOR,BENDING_MAGNET



class ParameterPlaneUndulator(Parameter):
    def __init__(self, K, E, lambda_u, L, I):
        super(self.__class__, self).__init__(E=E,I=I,type_magnet=PLANE_UNDULATOR)
        self.K = K
        self.lambda_u = lambda_u
        self.L = L

    def copy(self):
        return ParameterPlaneUndulator(K=self.K,E=self.E,lambda_u=self.lambda_u,L=self.L,I=self.I)

    def fct_magnetic_field(self, z, y,x, harmonic_number, coordonnee='y') :
        lambda_h= self.lambda_u / harmonic_number
        ku=2.0*np.pi/self.lambda_u

        if coordonnee=='x' :
            dB=0.0
        else :
            if coordonnee=='y' :
                #codata.m_e * codata.c / codata.e= 0.00170450894933
                Bo = (self.K * ku * 0.00170450894933)*np.cosh(ku*y)
                #print(Bo)
                f_base=np.cos
            else:#  coordonnee == 'z' :
                Bo=-(self.K * ku * 0.00170450894933)*np.sinh(ku*y)
                f_base=np.sin

            #we enelarge the real effect of the magnetic field by 4 lambda_u
            # on each extremity of the undulator
            L_magn_field=self.L/2.0+4.0*self.lambda_u

            # we know considere that if -(L/2 + lambda_u/4) < Z < (L/2 + lambda_u/4)
            # the magnetic field if a classic cosinus (or sinus)
            if coordonnee == 'y':
                a1=self.L/self.lambda_u-self.Nb_period()
                a2=(0.25-a1/2.0)
                L_cosinus_part=self.L/2.0 + self.lambda_u*a2

            else :
                L_cosinus_part = self.L / 2.0

            if ((z < -L_magn_field) or (z > L_magn_field)) :
                dB=0.0

            else :
                if (z < -L_cosinus_part or z> L_cosinus_part) :

                # in this case, we work with a gaussian,
                # so we shift the coordinate frame for use a central gaussian
                    if z < -L_cosinus_part :
                        sign=1
                    else : # z> L_cosinus_part
                        sign = -1

                    shift_z = z + sign * L_cosinus_part

                    p=2.0*np.pi**2/(3.0*lambda_h**2)
                    dB=((2.0*np.pi*Bo/lambda_h)*shift_z)*(1.0-4.0*np.pi**2*shift_z**2/(9.0*lambda_h**2))*np.exp(-p*shift_z**2)

                    # test du signe
                    z_test=sign*(-L_cosinus_part + lambda_h/4.0)
                    test_trig = f_base(ku*z_test)
                    if (sign * test_trig < 0) :
                        dB = -dB

                else :
                    # in this case we work in the cosinus part
                    dB=Bo*f_base(ku*z)

        return dB


    def get_Bo(self):
        return (2.0*np.pi*codata.m_e * codata.c*self.K)/( self.lambda_u*codata.e)

    def get_L(self):
        return  self.L

    def get_K(self):
        return self.K

    def theta_max(self):
        theta11=self.theta(1,1)
        return theta11

    def get_lambda_u(self):
        return self.lambda_u

    def Zmax_no_symetry(self):
        return self.get_L() + 10.0 * self.get_lambda_u()

    def Zo_symetry(self):
        return -(self.get_L() / 2.0 + 5.0 * self.get_lambda_u())

    def Zo_analitic(self):
        return -self.get_L() / 2.0



if __name__ == "__main__" :
    K = 1.87
    E = 1.3e9
    lambda_u = 0.035
    L = 0.035 * 12
    I = 1.0
    undulator=ParameterPlaneUndulator(K=K,E=E,lambda_u=lambda_u,L=L,I=I)
    print("The frequency for first harmotic is %f"%undulator.omega1())

    Bo=undulator.K / (93.4 * undulator.lambda_u)
    print(codata.e*Bo/(undulator.gamma()*codata.m))