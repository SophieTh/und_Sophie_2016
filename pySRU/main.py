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
from pySRU.RadiationFactory import RADIATION_METHOD_NEAR_FIELD, RADIATION_METHOD_APPROX_FARFIELD


eV_to_J=1.602176487e-19
######################################################

#BM_test=MagneticStructureBendingMagnet(Bo=0.8, div=5e-3, R=25.0)

E_1=7876.0

beam_test=ElectronBeam(Electron_energy=1.3, I_current=1.0)
beam_ESRF=ElectronBeam(Electron_energy=6.0, I_current=0.2)
und_test=Undulator(  K = 1.87,  period_length= 0.035, length=0.035 * 14)
ESRF18=Undulator( K = 1.68, period_length = 0.018, length=2.0)
ESRFBM=BM(Bo=0.8,L=0.1249994791673177)

#omega=E_1*eV_to_J/codata.hbar

# # recuperation des donnees de B en array en fonction de z
# reference=np.load("../tests/x_ray_booklet_field.npz")
# Z=reference['ct']
# Z -= (Z[-1])/2.0
# By=reference['B_y']


def compare_2_traj_method(magnetic_structure, electron_beam, rad_method, traj_method_ref, traj_method_test, distance=None):
#

    #TODO trouver un moyen de le faire directement sur le central cone...
    sim_test = create_simulation(magnetic_structure=magnetic_structure, electron_beam=electron_beam,
                             rad_method=rad_method, traj_method=traj_method_ref,distance=distance)
    sim_test.calculate_on_central_cone()


    ref = sim_test.copy()

    rad_max = (sim_test.radiation.max())
    print('rad max reference')
    print(rad_max)

    sim_test.change_trajectory_method(traj_method_test)

    traj_err=ref.trajectory.relativ_difference_with(sim_test.trajectory)

    #traj_err.plot()

    #sim_test.radiation.plot()
    #sim_test.radiation.plot()

    rad_err=ref.radiation.relativ_difference_with(sim_test.radiation)
    print('rad_max simulation test')
    print(rad_err.max())
    #rad_err.plot()

    spectre_ref,omega_array=sim_test.spectre()
    spectre_test=sim_test.spectre(omega_array=omega_array)[1]
    plt.plot(omega_array,spectre_ref)
    plt.plot(omega_array,spectre_test)
    plt.show()


def compare_interpolation(magnetic_structure, electron_beam,traj_method=TRAJECTORY_METHOD_ODE, distance=None):

    sim = create_simulation(magnetic_structure=magnetic_structure, electron_beam=electron_beam,
                            traj_method=traj_method,distance=100)

    z_ref = sim.trajectory.z * codata.c

    z0 = sim.trajectory.z[0] * codata.c
    z1 = sim.trajectory.z[-1] * codata.c
    z_test = np.linspace(z0, z1, 100)
    By_array = sim.source.magnetic_field.By(x=0.0, y=0.0, z=z_test)
    plt.plot(z_test, By_array)
    plt.show()

    By_linear = interpolate.interp1d(z_test, By_array, kind='linear')
    By_spline = interpolate.interp1d(z_test, By_array, kind='cubic')

    plt.plot(z_ref, sim.source.magnetic_field.By(x=0.0, y=0.0, z=z_ref))
    plt.plot(z_ref, By_linear(z_ref))
    plt.plot(z_ref, By_spline(z_ref))
    plt.show()


#

#
sim_test = create_simulation(magnetic_structure=ESRF18, electron_beam=beam_ESRF, traj_method=TRAJECTORY_METHOD_ANALYTIC,
                             rad_method=RADIATION_METHOD_APPROX_FARFIELD, Nb_pts_trajectory=10000,distance=100)


##test methode ODE et ANALitic sans les bord
# Nb_pts=np.linspace(1000,20000,200)
# error=sim_test.error_trajectory_method(TRAJECTORY_METHOD_ODE,nb_pts=Nb_pts)
# plt.plot(Nb_pts,error)
# plt.show()



#
#
# compare_2_traj_method(magnetic_structure=ESRF18,electron_beam=beam_ESRF,rad_method=RADIATION_METHOD_APPROX_FARFIELD,
#                       traj_method_ref=TRAJECTORY_METHOD_ANALYTIC,traj_method_test=TRAJECTORY_METHOD_ODE,distance=100)
# #
#
#
#
#
# # sim_test = create_simulation(magnetic_structure=und_test, electron_beam=beam_test,
# #                              traj_method=TRAJECTORY_METHOD_INTEGRATION,
# #                              Nb_pts_trajectory=10000)
# #
# # #sim_test.print_parameters()
# #
# # sim_test.change_omega(sim_test.source.critical_frequency()*0.5)
# #
# # print('omega')
# # print(sim_test.radiation_fact.omega)
# # print('omega critic 2')
# # print(sim_test.source.critical_frequency2())
# # print('arc length')
# # print(sim_test.source.arc_length()/codata.c)
# # print('physical length')
# # print((sim_test.trajectory.z[-1]-sim_test.trajectory.z[0]))
# # #
# #
# #
# # #sim_test.trajectory.plot()
# # print(type(sim_test.radiation))
# # sim_test.trajectory.plot_3D()
# # sim_test.radiation.plot()
# # # #
# # print('intensity[0][0]/1e13')
# # print(sim_test.radiation.intensity[0,0]/1e13)
# # print('radiation.max()/1e13')
# # print(sim_test.radiation.max()/1e13)
# # omega=sim_test.radiation_fact.omega
# # print('flux_on_axis_theoric(omega=omega)/1e13')
# # print(sim_test.source.flux_on_axis_theoric(omega=omega)/1e13)
# #
# # rad_theoric=sim_test.create_theoric_radiation()
# # print(type(rad_theoric))
# # print(sim_test.radiation.XY_are_like_in(rad_theoric))
# # print('rad_theoric.intensity[0,0]')
# # print(rad_theoric.intensity[0,0]/1e13)
# # rad_theoric.plot()
# #
# #
#
#
#
#
#
#
#
#
# #print(sim_test.source.magnetic_field_strength())
# z=sim_test.trajectory.z*codata.c
# x=sim_test.trajectory.x*codata.c
#
#
# # observation_angle=np.linspace(0.0,sim_test.source.angle_wave_number(1,2),401)
# # sim_test.calculate_for_observation_angles(observation_angle=observation_angle)
#
# #sim_test.change_distance(D=None)
# # nb_period=np.linspace(10,50,5)
# # print(nb_period)
# # error=sim_test.error_radiation_method_nb_period(method=RADIATION_METHOD_APPROX,nb_period=nb_period)
# #
# # plt.plot(nb_period,error)
# # plt.show()
#
#
# # sim_test.plot_magnetic_field_along_Z()
# # sim_test.trajectory.plot()
# # sim_test.trajectory.plot_3D()
# # sim_test.calculate_until_wave_number(harmonic_number=1,wave_number=2)
# sim_test.radiation.plot()
# # sim_test.radiation.plot_wave(Nb_pts=sim_test.radiation_fact.Nb_pts)
#
# # X_max=sim_test.radiation.X.max()
# # Y_max=sim_test.radiation.Y.max()
#
# sim_test.plot_spectre_on_axis()
# sim_test.plot_spectre_central_cone()

