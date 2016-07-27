import numpy as np
import scipy.constants as codata
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import scipy.integrate as integrate
import time
from pySRU.Trajectory import Trajectory
from pySRU.Radiation import Radiation
from pySRU.MagneticField import MagneticField
from pySRU.MagneticStructureUndulatorPlane import MagneticStructureUndulatorPlane  as Undulator
from pySRU.MagneticStructureBendingMagnet import MagneticStructureBendingMagnet  as BM
from pySRU.ElectronBeam import ElectronBeam
from pySRU.Source import Source
from pySRU.Simulation import Simulation ,create_simulation
from pySRU.TrajectoryFactory import TrajectoryFactory, TRAJECTORY_METHOD_ANALYTIC,TRAJECTORY_METHOD_ODE,\
                                        TRAJECTORY_METHOD_INTEGRATION
from pySRU.RadiationFactory import RadiationFactory ,RADIATION_METHOD_NEAR_FIELD, RADIATION_METHOD_APPROX,\
                                RADIATION_METHOD_FARFIELD, RADIATION_METHOD_APPROX_FARFIELD


eV_to_J=1.602176487e-19
######################################################



E_1=7876.0

omega=E_1*eV_to_J/codata.hbar
print("omega")
print(omega)


# # recuperation des donnees de B en array en fonction de z
# reference=np.load("../tests/x_ray_booklet_field.npz")
# Z=reference['ct']
# Z -= (Z[-1])/2.0
# By=reference['B_y']


#
# # # minimum  :
#
beam_test=ElectronBeam(E=1.3e9,I=1.0)
beam_ESRF=ElectronBeam(E=6.0e9,I=0.2)
und_test=Undulator(K = 1.87,lambda_u = 0.035,L=0.035*14)
ESRF18=Undulator(K = 1.68, lambda_u = 0.018, L=2.0)
ESRFBM=BM(Bo=0.8,div=5e-3,R=25.0)




def simul_und_analitic_formule_near_farfield(und_test,formule) :
    distance = 250
    theta_1_1 = np.sqrt((1.0 + und_test.K ** 2 / 2.0)) * (1.0 / und_test.gamma())
    Xmax = distance * theta_1_1
    Ymax = distance * theta_1_1

    traj_test=TrajectoryFactory(Nb_pts=2000,method=TRAJECTORY_METHOD_ODE)
    rad_test=RadiationFactory(method=RADIATION_METHOD_APPROX_FARFIELD,omega=und_test.omega1(),Nb_pts=100,formula=formule)
    sim_test = create_simulation(parameter=und_test, trajectory_fact=traj_test, radiation_fact=rad_test,
                                           distance=distance,X_max=Xmax,Y_max=Ymax)
    sim_test.magnetic_field.plot_z()
    print(sim_test.magnetic_field.z[0])
    print(sim_test.parameter.Bo())
    int1=integrate.quad((lambda z :sim_test.magnetic_field.By(z=z,y=0.0)),sim_test.magnetic_field.z[0],sim_test.magnetic_field.z[-1])
    print(int1)
    int2=integrate.quad((lambda z1 : integrate.quad((lambda z :sim_test.magnetic_field.By(z=z,y=0.0)),sim_test.magnetic_field.z[0],z1)[0] )
                   ,sim_test.magnetic_field.z[0],sim_test.magnetic_field.z[-1])
    print(int2)
    sim_test.trajectory.plot_3D()
    print(sim_test.radiation.max())
    sim_test.radiation.plot()

    sim_test.change_radiation_method(RADIATION_METHOD_FARFIELD)
    print(sim_test.radiation.max())
    sim_test.radiation.plot()
#simul_und_analitic_formule_near_farfield(und_test=und_test,formule=1)

def simul_BM_analitic_near_farfield(BM,formule) :
    distance = 1000
    Xmax = distance * 1e-2
    Ymax = distance * 1e-2

    traj_test=TrajectoryFactory(Nb_pts=2000,method=TRAJECTORY_METHOD_ANALYTIC)
    rad_test=RadiationFactory(method=RADIATION_METHOD_APPROX_FARFIELD,omega=BM.omega1(),Nb_pts=100,formula=formule)
    sim_test = create_simulation(parameter=BM, trajectory_fact=traj_test, radiation_fact=rad_test,
                                           distance=distance,X_max=Xmax,Y_max=Ymax)
    ref=sim_test.copy()

    sim_test.change_trajectory_method(TRAJECTORY_METHOD_INTEGRATION)
    sim_test.trajectory.plot_2_trajectory(ref.trajectory)
    err=sim_test.trajectory.error_rel_max(ref.trajectory)
    err.plot()
    print(ref.radiation.max())
    err2=ref.radiation.difference_with(sim_test.radiation)
    err2.plot()

def simul_und_analitic_erreur_near_farfield(und_test,formule) :
    distance = 100
    theta_1_1 = np.sqrt((1.0 + und_test.K ** 2 / 2.0)) * (1.0 / und_test.gamma())
    Xmax = distance * theta_1_1/4.0
    Ymax = distance * theta_1_1/4.0

    traj_test=TrajectoryFactory(Nb_pts=2000,method=TRAJECTORY_METHOD_ANALYTIC)
    rad_test=RadiationFactory(method=RADIATION_METHOD_APPROX_FARFIELD,omega=und_test.omega1(),Nb_pts=100,formula=formule)
    sim_test = create_simulation(parameter=und_test, trajectory_fact=traj_test, radiation_fact=rad_test,
                                           distance=distance,X_max=Xmax,Y_max=Ymax)
    D=np.linspace(und_test.D_max_plane_undulator(2)*0.5,und_test.D_max_plane_undulator(2)*2.0,10)
    err=sim_test.error_radiation_method_distance(RADIATION_METHOD_APPROX,D)
    plt.plot(D,err)
    plt.show()
#simul_und_analitic_erreur_near_farfield(und_test,1)

def n_min(und_test,alpha) :
    racine=(np.pi*und_test.Nb_period()*10**(alpha))/(45.)
    n_mini=(2.0*np.pi*und_test.Nb_period())*np.exp(0.25*np.log(racine))
    return n_mini

def n_min2(und_test,alpha) :
    racine=(10**(alpha)*und_test.L())/(6*und_test.Beta_et()*codata.c)
    n_mini=(2.0*np.pi*und_test.Nb_period())*np.exp(np.log(racine)/3.)
    return n_mini

def n_min3(und_test,alpha) :
    racine=((0.5-10**(-alpha))*180.*2.0*np.pi*und_test.Nb_period())
    n_mini=np.exp(0.2*np.log(1./racine))*(2.0*np.pi*und_test.Nb_period())
    return n_mini


#test1=simulation(ESRFBM,traj_method=TRAJECTORY_METHOD_ANALYTIC,omega=ESRFBM.omega1())

# spectre1,omega_array=test1.spectre()
# plt.plot(omega_array,spectre1)
# plt.show()

# spectre2,omega_array=test1.spectre2()
# plt.plot(omega_array,spectre2)
# plt.show()
#
# spectre3,omega_array=test1.spectre3(omega_array)
# plt.plot(omega_array,spectre3)
# plt.show()

test2=create_simulation(electron_beam=beam_test,magnetic_structure=und_test,
                        traj_method=TRAJECTORY_METHOD_ANALYTIC,rad_method=RADIATION_METHOD_APPROX_FARFIELD)
#test2.change_trajectory_method(TRAJECTORY_METHOD_ODE)

c1=test2.radiation.max()
print('rad max 2')
print(c1)
print('la theorie')
c2=source_test.flux_on_axis_theoric(n=1)
print(c2)
print('dif')
print((c1-c2)/c2)
# #
# test3=test2.copy()
# test3.change_trajectory_method(TRAJECTORY_METHOD_INTEGRATION)
# test2.trajectory.plot_2_trajectory(test3.trajectory)
# error=test2.trajectory.error_rel_max(test3.trajectory)
# print('error.v_x[-1]')
# print(error.v_x[-1])
# error.plot()
# print('radiation max test 3 ')
# print(test3.radiation.max())
# diff=test2.radiation.difference_with(test3.radiation)
# diff.plot()


# print('n min')
#
#
#
#
# alpha=4
# print(10**alpha)
# print("n mini 1")
# print(n_min(und_test=source_test,alpha=alpha))
# print(n_min(und_test=source_ESRF18,alpha=alpha))
# print("n mini2")
# print(n_min2(und_test=source_test,alpha=alpha))
# print(n_min2(und_test=source_ESRF18,alpha=alpha))
# print("n max")
# print(n_min3(und_test=source_test,alpha=alpha))
# print(n_min3(und_test=source_ESRF18,alpha=alpha))


# f=np.array([0.0,1./3.,2./3.,1.0])
# x=np.array([0.0,1./9.,4./9.,1.0])
# int=integrate.simps(f,x)
# print(int)
#
# rtue=und_test.Beta_et()*codata.c/(und_test.lambda_u*500)
# print(rtue)
# print(np.cos(2.0*np.pi*rtue))


#print(source_test.Beta_et()*codata.c*2.*np.pi/source_test.lambda_u())

