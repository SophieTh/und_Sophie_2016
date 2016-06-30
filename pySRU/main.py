import numpy as np
import scipy.constants as codata
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import time
from Trajectory import Trajectory
from Radiation import Radiation
from MagneticField import MagneticField
from UndulatorParameter import UndulatorParameters as Undulator
from UndulatorSimulation import UndulatorSimulation ,create_simulation
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
und_test=Undulator(K = 1.87,E = 1.3e9,lambda_u = 0.035,L=0.035*12,I=1.0)
ESRF18=Undulator(K = 1.68, E = 6.0e9,lambda_u = 0.018, L=0.018*111, I=0.2)
######### comparaison erreur trajectoire

distance=100
Xmax = distance * 1e-3
Ymax = distance * 1e-3


traj_test=TrajectoryFactory(Nb_pts=1001,method=TRAJECTORY_METHOD_ANALYTIC)
rad_test=RadiationFactory(method=RADIATION_METHOD_APPROX_FARFIELD,omega=und_test.omega1(),Nb_pts=101)
sim_test = create_simulation(undulator=und_test, trajectory_fact=traj_test, radiation_fact=rad_test, distance=distance,
                             X_max=Xmax, Y_max=Ymax)



initial_condition=np.array([0.0,0.0,ESRF18.Beta()*codata.c,0.0,0.0,-ESRF18.L/2.0-5*ESRF18.lambda_u])
traj=TrajectoryFactory(Nb_pts=20001,method=TRAJECTORY_METHOD_ODE, initial_condition=initial_condition)
rad=RadiationFactory(method=RADIATION_METHOD_APPROX_FARFIELD,omega=omega,Nb_pts=100)

Xmax = distance * 1e-4
Ymax = distance * 1e-4
sim_ESRF_1=create_simulation(undulator=ESRF18, trajectory_fact=traj, radiation_fact=rad, distance=distance,
                             X_max=Xmax, Y_max=Ymax)
print(sim_ESRF_1.radiation.max())
sim_ESRF_1.radiation.plot()

sim_ESRF_2=sim_ESRF_1.copy()
sim_ESRF_2.change_trajectory_method(TRAJECTORY_METHOD_ANALYTIC)
print(sim_ESRF_2.radiation.max())
sim_ESRF_2.radiation.plot()


