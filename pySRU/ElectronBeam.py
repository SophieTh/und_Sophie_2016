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


# container class
# parameter of the beam
class ElectronBeam(object):
    def __init__(self, Electron_energy, I_current):
        self.Electron_energy=Electron_energy
        self.I_current=I_current

    def copy(self):
        return ElectronBeam(Electron_energy=self.Electron_energy, I_current=self.I_current)

    def Lorentz_factor(self):
        return (self.Electron_energy / 0.51099890221)*1e3


    def electron_speed(self):
        gamma = self.Lorentz_factor()
        Beta = np.sqrt(1.0 - 1.0 / gamma ** 2)
        return Beta

    def print_parameters(self):
        print(' Electron Beam')
        print('    Electron energy (GeV): %.3f' %self.Electron_energy)
        print('    Intensity (A): %.3f' % self.I_current)
        print('    Electron speed : %.10f' %self.electron_speed())
        print('    Lorentz factor : %.7f' % self.Lorentz_factor())

if __name__ == "__main__" :
    K = 1.87
    E = 1.3
    I = 1.0
    beam=ElectronBeam(Electron_energy=E, I_current=I)

    beam.print_parameters()

    print('Lorentz factor: %f'%beam.Lorentz_factor())
