import numpy as np
import scipy.constants as codata
from Trajectory import Trajectory
from Radiation import Radiation
from MagneticField import MagneticField
from UndulatorParameter import UndulatorParameters as Undulator
from UndulatorSimulation import UndulatorSimulation
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
reference=np.load("x_ray_booklet_field.npz")
Z=reference['ct']
Z -= (Z[len(Z)-1])/2.0
By=reference['B_y']
SRW_magentic_field=MagneticField(x=0,y=0,z=Z,Bx=0,By=By,Bz=0)

#




# minimum  :
und_test=Undulator(K = 1.87,E = 1.3e9,lambda_u = 0.035,L=0.035*12,I=1.0)
traj_test=TrajectoryFactory(Nb_pts=201,method=TRAJECTORY_METHOD_ANALYTIC)

sim_test=UndulatorSimulation(undulator=und_test,trajectory_fact=traj_test)

# sim_test.trajectory.draw()
# sim_test.radiation.draw()

# with a given magnetic field
sim_SRW=UndulatorSimulation(undulator=und_test,trajectory_fact=traj_test,magnetic_field=SRW_magentic_field)
sim_SRW.change_trajectory_method(TRAJECTORY_METHOD_ODE)

# sim_SRW.trajectory.draw()
sim_SRW.radiation.draw()


# complete :
ESRF18=Undulator(K = 1.68, E = 6.0e9,lambda_u = 0.018, L=2.0, I=0.2)
# [Vx,Vy,Vz,x,y,z]
initial_condition=np.array([0.0,0.0,0.9999999997*codata.c,0.0,0.0,-1.0])
traj=TrajectoryFactory(Nb_pts=1001,method=TRAJECTORY_METHOD_ANALYTIC, initial_condition=initial_condition)
rad=RadiationFactory(method=RADIATION_METHOD_APPROX_FARFIELD,omega=omega)
X=np.linspace(0.0,0.005,101)
Y=np.linspace(0.0,0.005,101)

sim_ESRF=UndulatorSimulation(undulator=ESRF18,trajectory_fact=traj,radiation_fact=rad,
                             distance=30.0,X=X,Y=Y)

#sim_ESRF.trajectory.draw()
sim_ESRF.radiation.draw()



# c'est decale ????

