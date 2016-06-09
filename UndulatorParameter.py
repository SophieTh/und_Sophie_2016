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

    ## transcription SRW's code
    # we considere tht in the middle of the undulator z =0
    #and y=0 !!!!!!!!
    def SRW_fct_magnetic_field_plane_undulator(self,z,harmonic_nb) :
        Bo = self.K / (93.4 * self.lambda_u)

        lambda_h=self.lambda_u/harmonic_nb
        ku=2.0*np.pi/self.lambda_u

        #we enelarge the real effect of the magnetic field by 4 lambda_u
        # on each extremity of the undulator
        L_magn_field=self.L/2.0+4.0*self.lambda_u

        # we know considere that if -(L/2 + lambda_u/4) < Z < (L/2 + lambda_u/4)
        # the magnetic field if a classic cosinus
        L_cosinus_part=self.L/2.0 + self.lambda_u/4.0

        #for i_z, z in enumerate(Z) :
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
                test_trig = np.cos(ku*z_test)
                if (sign * test_trig < 0) :
                    dB = -dB

            else :
                #print('ok')
                # in this case we work in the cosinus part
                dB=Bo*np.cos(ku*z)

        return dB


    def SRW_magnetic_field(self,Z,harmonic_number) :
        if (type(Z)==np.ndarray) :
            By=np.zeros_like(Z)
            for i, Zi in enumerate (Z)  :
                By[i]=self.SRW_fct_magnetic_field_plane_undulator(Zi,harmonic_number)
        else :
            By=self.SRW_fct_magnetic_field_plane_undulator(Z,harmonic_number)
        return By

    def create_magnetic_field_plane_undulator(self,Z,harmonic_number):
        By = (lambda z: self.SRW_magnetic_field(Z=z, harmonic_number=harmonic_number))
        fct_null= (lambda x : 0.0)
        B = MagneticField(0.0, 0.0, Z, fct_null, By, fct_null)
        return B

    # old
    def creation_magnetic_field_plane_undulator(self, z):
        Bo = self.K / (93.4 * self.lambda_u)
        By=np.cos((2.0*np.pi/self.lambda_u)*z)

        # Hamming windowing
        windpar = 1.0 / (2.0 * np.floor(self.L / self.lambda_u))
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

        return B

    def gamma(self) :
        return self.E/0.511e6

    def omega1(self) :
        gamma=self.gamma()
        first_harm=((2.0 * gamma ** 2) / (1.0 + (self.K ** 2) / 2.0)) * ((2.0 * np.pi * codata.c) / self.lambda_u)
        return first_harm

    def Beta(self) :
        gamma=self.E/0.511e6
        Beta=np.sqrt(1.0-1.0/gamma**2)
        return Beta

    def D_max_plane_undulator(self,alpha):
        lim= self.L/2.0
        return lim *10**alpha


if __name__ == "__main__" :
    K = 1.87
    E = 1.3e9
    lambda_u = 0.035
    L = 0.035 * 12
    I = 1.0
    undulator=UndulatorParameters(K=K,E=E,lambda_u=lambda_u,L=L,I=I)
    print("The frequency for first harmotic is %f"%undulator.omega1())

    a=np.array([0.,1.,52.])
    print(type(a))
    print(type(a)== np.ndarray)
