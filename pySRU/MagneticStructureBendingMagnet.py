# coding: utf-8
# /*##########################################################################
#
# Copyright (c) 2016 European Synchrotron Radiation Facility
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#
# ###########################################################################*/

__authors__ = ["S Thery, M Glass, M Sanchez del Rio - ESRF ISDD Advanced Analysis and Modelling"]
__license__ = "MIT"
__date__ = "31/08/2016"

import numpy as np
import scipy.constants as codata
import scipy.special as special
from pySRU.MagneticStructure import MagneticStructure, PLANE_UNDULATOR,BENDING_MAGNET


#TODO changer et mettre la divergeance en parametre ?
class MagneticStructureBendingMagnet(MagneticStructure):
    def __init__(self, Bo, horizontale_divergeance,electron_energy=1e9):
        super(self.__class__, self).__init__(magnet_type=BENDING_MAGNET)
        self.horizontal_div=horizontale_divergeance
        self.radius_curv=electron_energy*1e9/(codata.c*Bo)
        self.length=np.sin(horizontale_divergeance)*self.radius_curv
        self.Bo = Bo


    def copy(self):
        return MagneticStructureBendingMagnet(length=self.length, Bo=self.Bo)

    def fct_magnetic_field(self, z, y, x, harmonic_number, coordonnee='y'):
        Zo = self.length * 0.6
        if coordonnee == 'y' and z >= -Zo and z <= Zo:
            # codata.m_e * codata.c / codata.e= 0.00170450894933
            Bo = self.Bo

        else:  # coordonnee == 'z' :
            Bo = 0.0
        return Bo


    def fct_magnetic_field2(self, z, y, x, harmonic_number, coordonnee='y'):
        Zo = self.length * 0.5
        if coordonnee == 'y':
            # codata.m_e * codata.c / codata.e= 0.00170450894933
            Bo = self.Bo
        else:  # coordonnee == 'z' :
            Bo = 0.0

            # we enelarge the real effect of the magnetic field by 4 lambda_u
            # on each extremity of the undulator
        L_magn_field = self.length * 1.5

        # we know considere that if -(L/2 + lambda_u/4) < Z < (L/2 + lambda_u/4)
        # the magnetic field if a classic cosinus (or sinus)
        L_cosinus_part = -Zo

        if ((z < -L_magn_field) or (z > L_magn_field)):
            dB = 0.0

        else:
            if (z < -L_cosinus_part or z > L_cosinus_part):

                # in this case, we work with a gaussian,
                # so we shift the coordinate frame for use a central gaussian
                if z > 0.0:
                    z_shift = z - L_cosinus_part
                else:
                    z_shift = z + L_cosinus_part

                sigma = (L_magn_field - L_cosinus_part) / 5.0
                # dB= Bo*(1.0-z_shift*self.get_L()**2/2.0)*np.exp(-0.5*(z_shift/sigma)**2)
                dB = Bo * np.exp(-0.5 * (z_shift / sigma) ** 2)
            else:
                # in this case we work in the cosinus part
                dB = Bo
        return dB



    def print_parameters(self):
        print(" Magnetic Structure (Bending Magnet):")
        print('    length : %.5f (m)'%self.length)
        print('    magnetic_field_strength : %.5f'%self.Bo)


if __name__ == "__main__" :

    BM = MagneticStructureBendingMagnet(Bo=0.85, length=1e-3 * 20.0)

    print('bending magnet test :')
    BM.print_parameters()
    B = BM.create_magnetic_field()
    Z=np.linspace(-BM.length/2.0,BM.length/2.0,10000)
    B.plot_z(Z=Z,Y=0.0,X=0.0)

