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
reference=np.load("x_ray_booklet_field.npz")
Z=reference['ct']
Z -= (Z[len(Z)-1])/2.0
By=reference['B_y']
SRW_magentic_field=MagneticField(x=0,y=0,z=Z,Bx=0,By=By,Bz=0)

Bo=1.87/(93.4*0.035)*np.ones_like(Z)
cst_magentic_field=MagneticField(x=0,y=0,z=Z,Bx=0,By=Bo,Bz=0)
#




# # minimum  :

und_test=Undulator(K = 1.87,E = 1.3e9,lambda_u = 0.035,L=0.035*12,I=1.0)
initial_cond=np.array([1e3,-1e5,np.sqrt((und_test.Beta()*codata.c)**2-(1e3)**2 -(-1e4)**2),1e-1,1e-3,-und_test.L/2.0-5.0*und_test.lambda_u])
traj_test=TrajectoryFactory(Nb_pts=1001,method=TRAJECTORY_METHOD_ODE,initial_condition=initial_cond)
rad_test=RadiationFactory(method=RADIATION_METHOD_FARFIELD,omega=und_test.omega1())

sim_test=create_simulation(undulator=und_test,trajectory_fact=traj_test,radiation_fact=rad_test)
# #
# #
#
#
#sim_test.trajectory.draw()

mag=sim_test.magnetic_filed
print(len(mag.z))
L_magn_field = und_test.L / 2.0 + 4.0 * und_test.lambda_u
print('L_magn_field')
print(L_magn_field)
L_cosinus_part = und_test.L / 2.0 + und_test.lambda_u / 4.0
print('L_cosinus_part')
print(L_cosinus_part)

#Z=mag.z
#Y=np.linspace(-2.0*1e-7,2.0*1e-7,len(Z))

Z=sim_test.trajectory.z
Y=sim_test.trajectory.y
#print(sim_test.trajectory.y*codata.c)
X=sim_test.trajectory.x

print(Z)
By=mag.By(Z,1e-3)
plt.plot(Z,By)
plt.show()



fig = plt.figure()
Bz=mag.Bz(Z,1e-3)
# Bo=-und_test.K / (93.4 *und_test.lambda_u)*np.sinh((2.0*np.pi/und_test.lambda_u)*1e-3)
# Bz2=Bo*np.sin((2.0*np.pi/und_test.lambda_u)*Z)
# plt.plot(Z,Bz2)
plt.plot(Z,Bz)
plt.show()
#
#
print('fini pour z')


fig = plt.figure()
By=mag.By(0.1,Y)
plt.plot(Y,By)
plt.show()

fig = plt.figure()
Bz=mag.Bz(0.1,Y)
plt.plot(Y,Bz)
plt.show()



sim_test.trajectory.plot_3D()

# print(Y)
# print(type(Z))
# print(type(Y))
# print(type(mag.By))
# print("ok")

# print("ok")
# Z,Y = np.meshgrid(Z,Y)
# fig = plt.figure()
# ax = Axes3D(fig)
# ax.plot_surface (Z,Y, By_array, rstride=1, cstride=1)
# ax.set_xlabel("Z")
# ax.set_ylabel('Y')
# ax.set_zlabel("By")
# plt.show()




#print(By_array)
# fig = plt.figure()
# ax = fig.gca(projection='3d')
# ax.plot(Y, Z,By, label='parametric curve')
# ax.legend()

plt.show()

#sim_test.trajectory.draw()
# print(sim_test.radiation.max())
#sim_test.radiation.draw()

# with a given magnetic field

# sim_SRW=UndulatorSimulation(undulator=und_test,trajectory_fact=traj_test,magnetic_field=SRW_magentic_field)
# sim_SRW.change_trajectory_method(TRAJECTORY_METHOD_ODE)
# #
# # sim_SRW.trajectory.draw()
# sim_SRW.radiation.draw()

# sim_cst=UndulatorSimulation(undulator=und_test,trajectory_fact=traj_test,magnetic_field=cst_magentic_field)
# sim_cst.change_trajectory_method(TRAJECTORY_METHOD_ODE)
#
# sim_cst.trajectory.draw()
# sim_cst.radiation.draw()


# print("calcul du spectre")
# omega1 = sim_test.undulator.omega1()
# omega_array = np.arange(omega1 * (1.0- 1e-2), omega1 * (1.0 + 1e-2), omega1*1e-4)
# sim_test.spectre_max(omega_array=omega_array)
#c'est decale ????


## complete :
#
# ESRF18=Undulator(K = 1.68, E = 6.0e9,lambda_u = 0.018, L=2.0, I=0.2)
# # [Vx,Vy,Vz,x,y,z]
# initial_condition=np.array([0.0,0.0,ESRF18.Beta()*codata.c,0.0,0.0,-1.0])
# traj=TrajectoryFactory(Nb_pts=501,method=TRAJECTORY_METHOD_ANALYTIC, initial_condition=initial_condition)
# rad=RadiationFactory(method=RADIATION_METHOD_APPROX_FARFIELD,omega=omega)
#
# #
#
#
# sim_ESRF=create_simulation(undulator=ESRF18,trajectory_fact=traj,radiation_fact=rad)
# print(sim_ESRF.radiation.X)
# # print(ESRF18.D_max_plane_undulator(alpha=2))
# # print(type(sim_ESRF.radiation_fact))
# # D=np.linspace(90.0,120.0,100)
# # error_method=sim_ESRF.error_radiation_method(method=RADIATION_METHOD_NEAR_FIELD,D=D)
# # print(sim_ESRF.radiation.max())
# # plt.plot(D,error_method)
# # plt.show()
# #
# # #sim_ESRF.trajectory.draw()
# sim_ESRF.radiation.draw()



