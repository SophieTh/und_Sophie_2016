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
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import scipy.integrate as integrate
import scipy.interpolate as interpolate
import time
from pySRU.MagneticStructureUndulatorPlane import MagneticStructureUndulatorPlane as Undulator
from pySRU.MagneticStructureBendingMagnet import MagneticStructureBendingMagnet  as BM
from pySRU.ElectronBeam import ElectronBeam
from pySRU.Source import Source
from pySRU.Simulation import Simulation ,create_simulation
from pySRU.TrajectoryFactory import TRAJECTORY_METHOD_ANALYTIC,TRAJECTORY_METHOD_ODE
from pySRU.RadiationFactory import RADIATION_METHOD_NEAR_FIELD, RADIATION_METHOD_APPROX_FARFIELD,RADIATION_METHOD_APPROX


eV_to_J=1.602176487e-19
######################################################

#BM_test=MagneticStructureBendingMagnet(Bo=0.8, div=5e-3, R=25.0)

E_1=7876.0

beam_test=ElectronBeam(Electron_energy=1.3, I_current=1.0)
beam_ESRF=ElectronBeam(Electron_energy=6.0, I_current=0.2)
und_test=Undulator(  K = 1.87,  period_length= 0.035, length=0.035 * 14)
ESRF18=Undulator( K = 1.68, period_length = 0.018, length=2.0)
ESRFBM=BM(Bo=0.8,horizontale_divergeance=0.005,electron_energy=6.0)



vx= 2e-4
vz= np.sqrt(beam_test.electron_speed()**2-vx**2)*codata.c

# initial_cond=np.array([ vx*codata.c,  0.00000000e+00 ,vz , 0.0 , 0.0 ,-0.42,])
#X=np.linspace(-0.02,0.02,150)
#Y=np.linspace(-0.02,0.02,150)
sim_test = create_simulation(magnetic_structure=und_test, electron_beam=beam_test, traj_method=TRAJECTORY_METHOD_ANALYTIC,
                            rad_method=RADIATION_METHOD_APPROX_FARFIELD, Nb_pts_trajectory=1000,distance=100)
sim_test.calculate_on_central_cone()
sim_test.change_energy_eV(E=2000)
#sim_test.change_harmonic_number(5.5)
#print(sim_test.source.choose_nb_pts_trajectory())
#print(sim_test.source.choose_distance_automatic())
#sim_test.print_parameters()


#sim_test.trajectory.plot_3D(title="Analytical electron trajectory in bending magnet")
#sim_test.trajectory.plot(title="Analytical electron trajectory in undulator")
#sim_test.calculate_on_central_cone()
#sim_test.calculate_spectrum_on_axis()
# omega1=sim_test.source.harmonic_frequency(1)
# omega_array=np.arange(.9*omega1,3.06*omega1,0.01*omega1)
# sim_test.calculate_spectrum_central_cone(abscissas_array=omega_array)
# sim_test.change_harmonic_number(3.)
#print (sim_test.radiation.max())
# sim_test.radiation.plot(title="radiation intensity from ODE trajectory solution ")
sim_test.radiation.plot(title="radiation intensity , approx ff 1 formula")



def erreur_near_ff(und,beam):
    sim_test_analy = create_simulation(magnetic_structure=und, electron_beam=beam,
                                       traj_method=TRAJECTORY_METHOD_ANALYTIC,Nb_pts_radiation=31,
                                       rad_method=RADIATION_METHOD_NEAR_FIELD,
                                       Nb_pts_trajectory=5000,distance=100)
    #sim_test_analy.change_harmonic_number(5.5)
    #sim_test_analy.calculate_on_central_cone()
    rad_max=sim_test_analy.radiation.max()
    print(rad_max)
    rad_max_theo=sim_test_analy.source.theoretical_flux_on_axis(frequency=sim_test_analy.radiation_fact.photon_frequency)
    print(rad_max_theo)
    d = np.linspace(20,60, 21)
    error_rad=sim_test_analy.error_radiation_method_distance(method=RADIATION_METHOD_APPROX,D=d)
    plt.plot(d,error_rad)
    plt.xlabel("distance")
    plt.ylabel("error")
    plt.title('absolut error between two radiation method')
    plt.show()
    plt.plot(d,error_rad/rad_max)
    plt.xlabel("distance")
    plt.ylabel("error")
    plt.title('relative error between two radiation method')
    plt.show()

#erreur_near_ff(ESRF18,beam_ESRF)
#erreur_near_ff(und_test,beam_test)

def erreur_nb_pts_traj(und,beam):
    sim_test_analy = create_simulation(magnetic_structure=und, electron_beam=beam,
                                       traj_method=TRAJECTORY_METHOD_ANALYTIC,Nb_pts_radiation=50,
                                       rad_method=RADIATION_METHOD_NEAR_FIELD,
                                       Nb_pts_trajectory=20000,distance=20)
    rad_ref=sim_test_analy.radiation.intensity
    rad_max=sim_test_analy.radiation.max()
    print(rad_max)
    rad_max_theo=sim_test_analy.source.theoretical_flux_on_axis(frequency=sim_test_analy.radiation_fact.photon_frequency)
    # print(rad_max_theo)
    # print((rad_max-rad_max_theo))
    np_pts= np.arange(300,500,10)
    error_rad=sim_test_analy.error_radiation_nb_pts_traj(good_value=rad_ref,nb_pts=np_pts)
    plt.plot(np_pts,error_rad)
    plt.xlabel("number of point")
    plt.ylabel("error")
    plt.title('absolut error, near field formula')
    plt.show()
    plt.plot(np_pts,error_rad/rad_max)
    plt.xlabel("number of point")
    plt.ylabel("error")
    plt.title('relative error, near field formula')
    plt.show()

#erreur_nb_pts_traj(und_test,beam_test)

