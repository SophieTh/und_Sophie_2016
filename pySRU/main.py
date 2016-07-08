import numpy as np
import scipy.constants as codata
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import scipy.integrate as integrate
import time
from Trajectory import Trajectory
from Radiation import Radiation
from MagneticField import MagneticField
from operator import *
from ParameterPlaneUndulator import ParameterPlaneUndulator  as Undulator
from ParameterBendingMagnet import ParameterBendingMagnet  as BM
from Simulation import Simulation ,create_simulation
from TrajectoryFactory import TrajectoryFactory, TRAJECTORY_METHOD_ANALYTIC,TRAJECTORY_METHOD_ODE,\
                                        TRAJECTORY_METHOD_INTEGRATION
from RadiationFactory import RadiationFactory ,RADIATION_METHOD_NEAR_FIELD, RADIATION_METHOD_APPROX,\
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
und_test=Undulator(K = 1.87,E = 1.3e9,lambda_u = 0.035,L=0.035*14,I=1.0)
ESRF18=Undulator(K = 1.68, E = 6.0e9,lambda_u = 0.018, L=2.0, I=0.2)
ESRFBM=BM(E=6.0e9,Bo=0.8,div=5e-3,R=25.0,I=0.2)


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

def simulation_undulator_method(und_test,traj_method,rad_method,distance=None,omega=None) :
    print('begin simulation')
    if distance==None :
        distance=und_test.D_max_plane_undulator(2)*2.0
    print('distance')
    print(distance)

    if omega==None :
        omega=und_test.omega1()
        print('omega')
    print(omega)

    num_har=np.floor(omega/und_test.omega1())
    if num_har==0 :
        num_har==1
    print('harmonic numb')
    print(num_har)
    theta_max=und_test.theta(n=num_har,l=1)
    Xmax = distance * theta_max
    Ymax = distance * theta_max

    Nb_pts=int(2.0*und_test.lambda_u*1e3*und_test.Nb_period())
    print('Nb point trajectory')
    print(Nb_pts)
    traj_test=TrajectoryFactory(Nb_pts=Nb_pts,method=traj_method)
    rad_test=RadiationFactory(method=rad_method,omega=omega,Nb_pts=100,formula=1)

    print('begin calcul')
    start_time=time.time()
    sim_test = create_simulation(parameter=und_test, trajectory_fact=traj_test, radiation_fact=rad_test,
                                           distance=distance,X_max=Xmax,Y_max=Ymax)
    delta_time=time.time()-start_time
    print('calcul time')
    print(delta_time)
    sim_test.magnetic_field.plot_z()
    sim_test.trajectory.plot_3D()
    print(sim_test.radiation.max())
    sim_test.radiation.plot()
    return sim_test

def simulation_undulator(und_test,distance=None,initial_condition=None,omega=None) :
    print('begin simulation')
    if distance==None :
        distance=und_test.D_max_plane_undulator(2)*2.0
    print('distance')
    print(distance)

    if omega==None :
        omega=und_test.omega1()
        print('omega')
    print(omega)

    num_har=np.floor(omega/und_test.omega1())
    if num_har==0 :
        num_har==1
    print('harmonic numb')
    print(num_har)
    theta_max=und_test.theta(n=num_har,l=1)
    Xmax = distance * theta_max
    Ymax = distance * theta_max

    Nb_pts=int(2.0*und_test.lambda_u*10**3*und_test.Nb_period())
    print('Nb point trajectory')
    print(Nb_pts)
    traj_test=TrajectoryFactory(Nb_pts=Nb_pts,method=TRAJECTORY_METHOD_ODE,initial_condition=initial_condition)
    rad_test=RadiationFactory(method=RADIATION_METHOD_APPROX_FARFIELD,omega=omega,Nb_pts=101,formula=1)

    print('begin calcul')
    start_time=time.time()
    sim_test = create_simulation(parameter=und_test, trajectory_fact=traj_test, radiation_fact=rad_test,
                                           distance=distance,X_max=Xmax,Y_max=Ymax)
    delta_time=time.time()-start_time
    print('calcul time')
    print(delta_time)
    sim_test.magnetic_field.plot_z()
    sim_test.trajectory.plot_3D()
    print(sim_test.radiation.max())
    sim_test.radiation.plot()
    return sim_test


#simulation_undulator(und_test=und_test)


test1=simulation_undulator_method(ESRF18,traj_method=TRAJECTORY_METHOD_ODE,
                                  rad_method=RADIATION_METHOD_APPROX_FARFIELD)

# test2=test1.copy()
# test2.change_trajectory_method(TRAJECTORY_METHOD_INTEGRATION)
#
# test2.radiation.plot()
#
#
# err=test1.radiation.difference_with(test2.radiation)
# err.plot()
# err.intensity *= (1.0/test1.radiation.max())
# err.plot()