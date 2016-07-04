import numpy as np
import scipy.constants as codata
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
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
from RadiationFactory import RadiationFactory ,RADIATION_METHOD_NEAR_FIELD, \
                                RADIATION_METHOD_FARFIELD, RADIATION_METHOD_APPROX_FARFIELD


eV_to_J=1.602176487e-19
######################################################



E_1=7876.0

omega=E_1*eV_to_J/codata.hbar
print("omega")
print(omega)


# recuperation des donnees de B en array en fonction de z
reference=np.load("../tests/x_ray_booklet_field.npz")
Z=reference['ct']
Z -= (Z[-1])/2.0
By=reference['B_y']
#SRW_magentic_field=MagneticField(x=0,y=0,z=Z,Bx=0,By=By,Bz=0)

# #Bo=1.87/(93.4*0.035)*np.ones_like(Z)
# cst_magentic_field=MagneticField(x=0,y=0,z=Z,Bx=,By=Bo,Bz=0)
#



#
# # # minimum  :
#
und_test=Undulator(K = 1.87,E = 1.3e9,lambda_u = 0.035,L=0.035*14,I=1.0)
ESRF18=Undulator(K = 1.68, E = 6.0e9,lambda_u = 0.018, L=0.018*111, I=0.2)
ESRFBM=BM(E=6.0e9,Bo=0.8,div=5e-3,R=25.0,I=0.2)




# bending magnet

# print(ESRFBM.L())
# traj_fact=TrajectoryFactory(Nb_pts=2000,method=TRAJECTORY_METHOD_ANALYTIC)
# rad_fact=RadiationFactory(method=RADIATION_METHOD_APPROX_FARFIELD,omega=ESRF18.omega1(),Nb_pts=100,formula=2)
# sim_BM=create_simulation(parameter=ESRFBM,trajectory_fact=traj_fact,radiation_fact=rad_fact,
#                          distance=distance,X_max=Xmax,Y_max=Ymax)
#
# sim_BM.trajectory.plot()
#
# sim_BM.trajectory.plot_3D()
# print(sim_BM.radiation.max())
# sim_BM.radiation.plot()
#
# # sim_BM.change_trajectory_method(TRAJECTORY_METHOD_ODE)
# # sim_BM.trajectory.plot_3D()
# # print(sim_BM.radiation.max())
# # sim_BM.radiation.plot()
#
# sim_BM.change_radiation_method(RADIATION_METHOD_FARFIELD)
# print(sim_BM.radiation.max())
# sim_BM.radiation.plot()


# undulator

def simul_und_analitic_formule_near_farfield(und_test,formule) :
    distance = 100
    theta_1_1 = np.sqrt((1.0 + und_test.K ** 2 / 2.0)) * (1.0 / und_test.gamma())
    Xmax = distance * theta_1_1
    Ymax = distance * theta_1_1

    traj_test=TrajectoryFactory(Nb_pts=2000,method=TRAJECTORY_METHOD_ANALYTIC)
    rad_test=RadiationFactory(method=RADIATION_METHOD_APPROX_FARFIELD,omega=und_test.omega1(),Nb_pts=100,formula=formule)
    sim_test = create_simulation(parameter=und_test, trajectory_fact=traj_test, radiation_fact=rad_test,
                                           distance=distance,X_max=Xmax,Y_max=Ymax)
    sim_test.magnetic_field.plot_z()
    sim_test.trajectory.plot_3D()
    print(sim_test.radiation.max())
    sim_test.radiation.plot()

    sim_test.change_radiation_method(RADIATION_METHOD_FARFIELD)
    print(sim_test.radiation.max())
    sim_test.radiation.plot()

#simul_und_analitic_formule_near_farfield(und_test=ESRF18,formule=1)

def simul_BM_analitic_near_farfield(BM,formule) :
    distance = 100
    Xmax = distance * 1e-2
    Ymax = distance * 1e-2

    traj_test=TrajectoryFactory(Nb_pts=2000,method=TRAJECTORY_METHOD_ANALYTIC)
    rad_test=RadiationFactory(method=RADIATION_METHOD_APPROX_FARFIELD,omega=BM.omega1(),Nb_pts=100,formula=formule)
    sim_test = create_simulation(parameter=BM, trajectory_fact=traj_test, radiation_fact=rad_test,
                                           distance=distance,X_max=Xmax,Y_max=Ymax)
    sim_test.magnetic_field.plot_z()
    sim_test.trajectory.plot_3D()
    print(sim_test.radiation.max())
    sim_test.radiation.plot()

    sim_test.change_radiation_method(RADIATION_METHOD_FARFIELD)
    print(sim_test.radiation.max())
    sim_test.radiation.plot()

simul_BM_analitic_near_farfield(BM=ESRFBM,formule=1)
#
# traj_test2=TrajectoryFactory(Nb_pts=2000,method=TRAJECTORY_METHOD_ODE,
#                              initial_condition=sim_test.trajectory_fact.initial_condition)
#
# sim_test2= create_simulation(parameter=ESRF18, trajectory_fact=traj_test2, radiation_fact=rad_test,
#                                        distance=distance,X_max=Xmax,Y_max=Ymax)
#
# sim_test2.magnetic_field.plot_z()
# sim_test2.trajectory.plot_3D()
# print(sim_test2.radiation.max())
# sim_test2.radiation.plot()
#
#
# sim_test2.change_trajectory_method(TRAJECTORY_METHOD_INTEGRATION)
# #sim_test2.magnetic_field.plot_z()
# sim_test2.trajectory.plot_3D()
# print(sim_test2.radiation.max())
# sim_test2.radiation.plot()
#
# initial_cond=np.array([0.0,0.0,ESRF18.Beta()*codata.c,0.0,0.0,ESRF18.Zo_symetry()])
# sim_test2.change_initial_condition(initial_cond=initial_cond)
# sim_test2.magnetic_field.plot_z()
# sim_test2.trajectory.plot_3D()
# print(sim_test2.radiation.max())
# sim_test2.radiation.plot()
#
# sim_test2.change_trajectory_method(TRAJECTORY_METHOD_ODE)
# sim_test2.trajectory.plot_3D()
# print(sim_test2.radiation.max())
# sim_test2.radiation.plot()
#
#
# err=sim_test.radiation.difference_with(sim_test2.radiation)
# err.plot()


#near field

# start_time=time.time()
# sim_test2.change_radiation_method(RADIATION_METHOD_FARFIELD)
# delta_t=time.time()-start_time
# print(delta_t)
# print(sim_test2.radiation.max())
# sim_test2.radiation.plot()
#
# start_time=time.time()
# sim_test2.change_radiation_method(RADIATION_METHOD_NEAR_FIELD)
# delta_t=time.time()-start_time
# print(delta_t)
# print(sim_test2.radiation.max())
# sim_test2.radiation.plot()



# sim_test.change_trajectory_method(TRAJECTORY_METHOD_INTEGRATION)
# nb_pts = np.arange(500, 2001, 500, dtype='int')
# err = sim_test.error_radiation_method_nb_pts_traj(RADIATION_METHOD_APPROX_FARFIELD, nb_pts=nb_pts)
# plt.plot(nb_pts,err)
# plt.show()



# initial_condition=np.array([0.0,0.0,ESRF18.Beta()*codata.c,0.0,0.0,-ESRF18.L/2.0-5*ESRF18.lambda_u])
# traj=TrajectoryFactory(Nb_pts=20001,method=TRAJECTORY_METHOD_ODE, initial_condition=initial_condition)
# rad=RadiationFactory(method=RADIATION_METHOD_APPROX_FARFIELD,omega=omega,Nb_pts=100)
#
# Xmax = distance * 1e-4
# Ymax = distance * 1e-4
# sim_ESRF_1=create_simulation(undulator=ESRF18, trajectory_fact=traj, radiation_fact=rad, distance=distance,
#                              X_max=Xmax, Y_max=Ymax)
# print(sim_ESRF_1.radiation.max())
# sim_ESRF_1.radiation.plot()
#
# sim_ESRF_2=sim_ESRF_1.copy()
# sim_ESRF_2.change_trajectory_method(TRAJECTORY_METHOD_ANALYTIC)
# print(sim_ESRF_2.radiation.max())
# sim_ESRF_2.radiation.plot()
#
# sim_cst=create_simulation_cst_magnecticfield(undulator=ESRF18, trajectory_fact=traj_test, radiation_fact=rad_test, distance=distance,
#                              X_max=Xmax, Y_max=Ymax)
#
# sim_cst.trajectory.plot_3D()
# sim_cst.radiation.plot()
# D=np.arange(1000,2000,200)
# print(len(D))
# err=sim_cst.error_radiation_method_distance(RADIATION_METHOD_FARFIELD,D=D)
# plt.plot(D,err)
# plt.show()
