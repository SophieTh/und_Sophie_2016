import numpy as np
import scipy.constants as codata
from scipy.special import jn,yn,jv,yv
from pySRU.MagneticField import MagneticField
from pySRU.MagneticStructure import MagneticStructure , PLANE_UNDULATOR,BENDING_MAGNET



class MagneticStructureUndulatorPlane(MagneticStructure):
    def __init__(self, K, period_length, length):
        super(self.__class__, self).__init__(magnet_type=PLANE_UNDULATOR)
        self.K = K
        self.period_length = period_length
        self.length=length

    def copy(self):
        return MagneticStructureUndulatorPlane(K=self.K, period_length=self.period_length, length=self.length)

    def fct_magnetic_field(self, z, y, x, harmonic_number, coordonnee='y'):
        lambda_h = self.period_length / harmonic_number
        ku = 2.0 * np.pi / self.period_length

        if coordonnee == 'x':
            dB = 0.0
        else:
            if coordonnee == 'y':
                # codata.m_e * codata.c / codata.e= 0.00170450894933
                Bo = (self.K * ku * 0.00170450894933) * np.cosh(ku * y)
                # print(Bo)
                f_base = np.cos
            else:  # coordonnee == 'z' :
                Bo = -(self.K * ku * 0.00170450894933) * np.sinh(ku * y)
                f_base = np.sin

            # we enelarge the real effect of the magnetic field by 4 lambda_u
            # on each extremity of the undulator
            L_magn_field = self.length / 2.0 + 4.0 * self.period_length

            # we know considere that if -(L/2 + lambda_u/4) < Z < (L/2 + lambda_u/4)
            # the magnetic field if a classic cosinus (or sinus)
            if coordonnee == 'y':
                a1 = self.length / self.period_length \
                     - np.floor(self.period_number())
                a2 = (0.25 - a1 / 2.0)
                L_cosinus_part = self.length / 2.0 \
                                 + self.period_length * a2

            else:
                L_cosinus_part = self.length / 2.0

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
                    dB = ((2.0 * np.pi * Bo / lambda_h) * shift_z) * (
                    1.0 - 4.0 * np.pi ** 2 * shift_z ** 2 / (9.0 * lambda_h ** 2)
                    ) * np.exp(-p * shift_z ** 2)

                    # test du signe
                    z_test = sign * (-L_cosinus_part + lambda_h / 4.0)
                    test_trig = f_base(ku * z_test)
                    if (sign * test_trig < 0):
                        dB = -dB

                else:
                    # in this case we work in the cosinus part
                    dB = Bo * f_base(ku * z)

        return dB


    def print_parameters(self):
        print(' Magnetic Structure (Plane Undulator)')
        print('    K : %.2f'%(self.K))
        print('    period lenght : %.3f' % (self.period_length))
        print('    lenght : %.2f'% (self.length))
        print('    number of periods : %d (L/lambdau=%.2f)' % (int(self.period_number()),self.period_number()))
        print('    magntic field intensity: %.2f'% (self.magnetic_field_strength()))

    def period_number(self):
        return self.length / self.period_length

    def magnetic_field_strength(self):
        Bo = (2.0 * np.pi * codata.m_e * codata.c * self.K) / (
            self.period_length * codata.e)
        return Bo

if __name__ == "__main__" :
    und_test = MagneticStructureUndulatorPlane(K=1.87, period_length=0.035, length=0.035 * 14)
    ESRF18 = MagneticStructureUndulatorPlane(K=1.68, period_length=0.018, length=2.0)

    print('undulator test :')
    und_test.print_parameters()
    B=und_test.create_magnetic_field()
    Z=np.linspace(-und_test.length/2.0,und_test.length/2.0,10000)
    B.plot_z(Z=Z,Y=0.0,X=0.0)

    print('undulator ESRF18 :')
    ESRF18.print_parameters()
    B=ESRF18.create_magnetic_field()
    Z=np.linspace(-ESRF18.length/2.0,ESRF18.length/2.0,10000)
    B.plot_z(Z=Z,Y=0.0,X=0.0)