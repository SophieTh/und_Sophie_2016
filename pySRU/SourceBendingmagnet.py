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
from pySRU.Source import Source,PLANE_UNDULATOR,BENDING_MAGNET
from pySRU.ElectronBeam import ElectronBeam



#TODO
class SourceBendingMagnet(Source):
    def __init__(self, electron_beam, magnetic_structure, magnetic_field=None):
        super(self.__class__, self).__init__(electron_beam=electron_beam,
                                             magnetic_field=magnetic_field,
                                             magnetic_structure=magnetic_structure)


    def magnetic_field_strength(self):
        return self.magnetic_structure.Bo


    def copy(self):
        return SourceBendingMagnet(electron_beam=self.electron_beam.copy(),
                                    magnetic_field=self.magnetic_field.copy(),
                                    magnetic_structure=self.magnetic_structure)



    def arc_length(self):
        return self.magnetic_structure.horizontal_div* self.magnetic_structure.radius_curv

    def horizontal_divergence(self):
        return self.magnetic_structure.horizontal_div


    def radius_curvature(self,):
        return self.magnetic_structure.radius_curv



    #TODO faire le tri
    def critical_energy (self):
        Ec= 0.665 * (self.Electron_energy() ) ** 2 * self.magnetic_structure.Bo
        return Ec*1e-6

    def critical_frequency2(self):
        critical_energy=self.critical_energy()
        return critical_energy/codata.hbar

    def critical_frequency(self):
        gamma=self.Lorentz_factor()
        res=3.*(gamma**3)*codata.c/(2.*self.radius_curvature())
        return res

    def H2(self,y):
        K23=special.kv((2./3.),0.5*y)
        return (y*K23)**2


    def theoretical_flux_on_axis(self,n):
        #n=frequency/self.critical_frequency()
        #print('y= %.3f'%y)
        result= 1.327e13*((self.Electron_energy())**2
                          )*self.I_current()*self.H2(n)
        return result

    def radiation_theoric(self,omega,observation_angle):
        gamma=self.Lorentz_factor()
        X=gamma*observation_angle
        y=omega/self.critical_frequency()
        xi=y*0.5*np.sqrt((1.+X**2)**3)
        cst=(3.*codata.alpha*(gamma**2)*1e-3*1e-6*self.I_current()*y**2)/(codata.e*4.*np.pi**2)
        rad=((1.+X**2)**2)*((special.kv((2./3.),xi))**2+((X**2)/(1.+X**2))*(special.kv((2./3.),xi))**2)
        return rad*cst






    def choose_distance_automatic(self, alpha=2,photon_frequency=None):
        return self.magnetic_structure.length*10**(alpha)


    def choose_nb_pts_trajectory(self, alpha=0.01,photon_frequency=None):
        return self.magnetic_structure.length*4.*10**(alpha)


    def choose_initial_contidion_automatic(self):
        Zo=-self.magnetic_structure.length*1.5
        ic=np.array([0.0,0.0,self.electron_speed()*codata.c,0.0,0.0,Zo])
        return ic

    def choose_photon_frequency(self):
        return self.critical_frequency()*0.1


    def choose_angle_deflection_max(self):
        theta=5.0/self.Lorentz_factor()
        return theta


    def angle_deflection_central_cone(self):
        return 1.0/self.Lorentz_factor()


    def analytical_times_vector(self, Nb_pts):
        to=-self.magnetic_structure.length*0.5/(self.electron_speed()*codata.c)
        t1=to+self.arc_length()/(self.electron_speed()*codata.c)
        time=np.linspace(to,t1,Nb_pts)
        return time


    def construct_times_vector(self, initial_contition, Nb_pts):
    #TODO a changer ne marche pas ???
        electron_speed=(self.electron_speed() * codata.c)
        to = -self.magnetic_structure.length * 0.5 / electron_speed
        t1 = to + self.arc_length() / electron_speed
        time_start=initial_contition[5]/ electron_speed
        if time_start <to :
            time=np.linspace(time_start,-time_start,Nb_pts)
        else :
            delta_t=t1-time_start
            if delta_t >0.0 :
                time=np.linspace(time_start,time_start+2.*delta_t,Nb_pts)
            else :#TODO a tester
                time=np.linspace(time_start,2.*time_start,Nb_pts)
        return time

    def rtol_for_ODE_method(self):
        return 1e-11

    def atol_for_ODE_method(self):
        gamma=self.Lorentz_factor()
        atol_vx= self.electron_speed()*codata.c*1e-11
        atol_vz = self.electron_speed()*codata.c*1e-11
        atol_x = (2.*self.radius_curvature())/(3.*gamma**2)*1e-2
        atol_z = (2.*self.radius_curvature())/(3.*gamma**2)*1e-2
        return np.array([atol_vx,1e-10,atol_vz,atol_x,1e-10,atol_z])

# source parameter

    def print_parameters(self):
        super(self.__class__, self).print_parameters()
        print('Bending Magnet')
        print('    length : %.5f (m)'%self.magnetic_structure.length)
        print('    magnetic_field_strength : %.5f (T)'%self.magnetic_structure.Bo)
        print('    horizontal divergeance : %.5f (rad ?)' % self.horizontal_divergence())
        print('    radius curvature : %.3f (m)' % self.radius_curvature())
        print('    critical frequency : %f *1e20 (unite ?)' %(self.critical_frequency()/1e20))





if __name__ == "__main__" :
    pass




