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
    def SRW_fct_magnetic_field_plane_undulator(self, z, y, harmonic_number, coordonnee='y') :
        lambda_h= self.lambda_u / harmonic_number
        ku=2.0*np.pi/self.lambda_u

        if coordonnee=='y' :
            #print('coordonnee y')
            Bo = self.K / (93.4 * self.lambda_u)*np.cosh(ku*y)
            #print(Bo)
            f_base=np.cos
        else : # coordonne = z
            Bo=-self.K / (93.4 * self.lambda_u)*np.sinh(ku*y)
            f_base=np.sin

        #we enelarge the real effect of the magnetic field by 4 lambda_u
        # on each extremity of the undulator
        L_magn_field=self.L/2.0+4.0*self.lambda_u

        # we know considere that if -(L/2 + lambda_u/4) < Z < (L/2 + lambda_u/4)
        # the magnetic field if a classic cosinus (or sinus)
        if coordonnee == 'y':
            L_cosinus_part=self.L/2.0 + self.lambda_u/4.0

        else :
            L_cosinus_part = self.L / 2.0

        #for i_z, z in enumerate(Z) :
        if ((z < -L_magn_field) or (z > L_magn_field)) :
            dB=0.0

        else :
            if (z < -L_cosinus_part or z> L_cosinus_part) :
                #print('hors cos part')
            # in this case, we work with a gaussian,
            # so we shift the coordinate frame for use a central gaussian
                if z < -L_cosinus_part :
                    sign=1
                else : # z> L_cosinus_part
                    sign = -1

                shift_z = z + sign * L_cosinus_part

                p=2.0*np.pi**2/(3.0*lambda_h**2)
                dB=((2.0*np.pi*Bo/lambda_h)*shift_z)*(1.0-4.0*np.pi**2*shift_z**2/(9.0*lambda_h**2))*np.exp(-p*shift_z**2)
                #print(dB)
                # test du signe
                z_test=sign*(-L_cosinus_part + lambda_h/4.0)
                test_trig = f_base(ku*z_test)
                if (sign * test_trig < 0) :
                    dB = -dB
                #print(dB)
            else :
                # print(' ')
                # print('ok')
                # print(' ')
                # in this case we work in the cosinus part
                dB=Bo*f_base(ku*z)
                # print(ku*z)
                # print(f_base(ku * z))
                # print(Bo)
                # print(dB)


        return dB

#### essai
    def SRW_fct_magnetic_field_plane_undulator2(self, z, y, harmonic_number, coordonnee='y'):
        lambda_h = self.lambda_u / harmonic_number
        ku = 2.0 * np.pi / self.lambda_u

        if coordonnee == 'y':
            # print('coordonnee y')
            Bo = self.K / (93.4 * self.lambda_u) * np.cosh(ku * y)
            # print(Bo)
            f_base = np.cos
        else:  # coordonne = z
            Bo = -self.K / (93.4 * self.lambda_u) * np.sinh(ku * y)
            f_base = np.sin

        # we enelarge the real effect of the magnetic field by 4 lambda_u
        # on each extremity of the undulator
        L_magn_field = self.L / 2.0 + 4.0 * self.lambda_u

        # we know considere that if -(L/2 + lambda_u/4) < Z < (L/2 + lambda_u/4)
        # the magnetic field if a classic cosinus (or sinus)
        if coordonnee == 'y':
            L_cosinus_part = self.L / 2.0 + self.lambda_u / 4.0

        else:
            L_cosinus_part = self.L / 2.0

        # for i_z, z in enumerate(Z) :
        if ((z < -L_magn_field) or (z > L_magn_field)):
            dB = 0.0

        else:
            if (z < -L_cosinus_part or z > L_cosinus_part):
                # print('hors cos part')
                # in this case, we work with a gaussian,
                # so we shift the coordinate frame for use a central gaussian
                if z < -L_cosinus_part:
                    sign = 1
                else:  # z> L_cosinus_part
                    sign = -1
                shift_z = z + sign * L_cosinus_part
                sigma2=(self.lambda_u)**2
                dB = (Bo*f_base(ku*z)/(np.sqrt(2.0*np.pi*sigma2)) * np.exp(-0.5*shift_z**2/sigma2))
                if (Bo*f_base(ku*z)*dB < 0):
                    dB = -dB
                # sigma2=((5.0*self.lambda_u)/16.0)**2
                # dB = ( Bo*shift_z*f_base(ku*z)/(np.sqrt(2.0*np.pi*sigma2)) )* np.exp(-0.5*shift_z ** 2/sigma2)
                # # print(dB)
                # # test du signe)
                # z_test = sign * (-L_cosinus_part + lambda_h / 4.0)
                # test_trig = f_base(ku * z_test)
                # if (sign * test_trig < 0):
                #     dB = -dB
            else:
                dB = Bo * f_base(ku * z)
                # print(ku*z)
                # print(f_base(ku * z))
                # print(Bo)
                # print(dB)

        return dB


    def SRW_magnetic_field(self,Z,Y,harmonic_number,coordonnee) :
        if (type(Z)==np.ndarray) :
            B=np.zeros_like(Z)
        else :
            if (type(Y) == np.ndarray ):
                B=np.zeros_like(Y)
            else :
                B = self.SRW_fct_magnetic_field_plane_undulator(Z, Y,
                                            harmonic_number = harmonic_number, coordonnee = coordonnee)

        if (type(Z)==np.ndarray) :
            #print('entree ici')
            #B=np.zeros_like(Z)
            if (type(Y)==np.ndarray) :
                for i, Zi in enumerate(Z):
                            B[i] = self.SRW_fct_magnetic_field_plane_undulator(Zi,Y[i],
                                            harmonic_number=harmonic_number,coordonnee=coordonnee)
            else :
                for i, Zi in enumerate (Z)  :
                    B[i]=self.SRW_fct_magnetic_field_plane_undulator(Zi,Y,
                                                        harmonic_number=harmonic_number,coordonnee=coordonnee)
        else :
            # if (type(Y) != np.ndarray) :
            #         B=self.SRW_fct_magnetic_field_plane_undulator(Z,Y,
            #                                             harmonic_number=harmonic_number, coordonnee=coordonnee)
            if (type(Y) == np.ndarray ) :
                #print('Y vector')
                for i, Yi in enumerate (Y)  :
                    #B = np.zeros_like(Y)
                    B[i]=self.SRW_fct_magnetic_field_plane_undulator(Z,Yi,
                                                        harmonic_number=harmonic_number,coordonnee=coordonnee)
                   # print(B[i])
        #print('B dans srw')
        #print(B)
        return B

    def create_magnetic_field_plane_undulator(self,Z,Y,harmonic_number):
        By = (lambda z,y: self.SRW_magnetic_field(Z=z,Y=y,harmonic_number=harmonic_number,coordonnee='y'))
        Bz = (lambda z, y: self.SRW_magnetic_field(Z=z, Y=y, harmonic_number=harmonic_number, coordonnee='z'))
        fct_null= (lambda z,y : 0.0)
        #B = MagneticField(0.0, Y, Z, fct_null, By, fct_null)
        B = MagneticField(0.0, Y, Z, fct_null, By, Bz)
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
