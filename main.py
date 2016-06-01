import numpy as np
import matplotlib.pyplot as plt
from pylab import *
from mpl_toolkits.mplot3d import Axes3D
import scipy.constants as codata
import Trajectory
import Radiation
from TrajectoryFactory import TrajectoryFactory
from RadiationFactory import UndulatorRadiationFactory

from calc_undulator_2 import radiation_single_electron, draw_trajectory,undulator_trajectory

######################################################


K = 1.87
E = 1.3e9
lambda_u = 0.035
Nb_period = 12
Nb_pts = 20*12+1

gamma = E /0.511e6
print(gamma**2)
Beta = np.sqrt(1.0 - (1.0 / gamma ** 2))
print(1e0-Beta)
ku = 2.0 * np.pi / lambda_u
gamma = E / 0.511e6
Beta = np.sqrt(1.0 - (1.0 / gamma ** 2))
Beta_et = 1.0 - (1.0 / (2.0 * gamma ** 2)) * (1.0 + (K ** 2) / 2.0)
omega1 = ((2.0 * gamma ** 2) / (1.0 + (K ** 2) / 2.0)) * ((2.0 * np.pi * codata.c) / lambda_u)
L=Nb_period*lambda_u

# recuperation des donnees de B en array en fonction de z
reference=np.load("x_ray_booklet_field.npz")
print(reference.keys())
Z=reference['ct']
Z -= (Z[len(Z)-1])/2.0
By2=reference['B_y']

Z_By=np.zeros((2,len(Z)))
Z_By[0]=Z
Z_By[1]=By2

cst=codata.e * 1e-10 / (8.0 * np.pi ** 2 * codata.epsilon_0 * codata.c * codata.h)
print("cst")
print(cst)
# print(Beta_et*codata.c*ku*K/gamma)
# print(codata.c/gamma**2)
# print(lambda_u*codata.e*1e-10*L/(2.0*np.pi**2*codata.epsilon_0*codata.h*4.0*codata.c**2*(1.0-Beta)**2 *gamma))
# print (codata.epsilon_0)
# print (codata.h)
# print (codata.c**2)
# print (codata.e)
#
cst2=(2*ku*K*codata.c/(gamma**3*(1.0-Beta)**4))
print("cst2")
print(cst2)

cst3=cst2*cst*(2.0*30.0-L)**2/30.0
print("cst3")
print(cst3)
T=TrajectoryFactory(K,E,lambda_u,Nb_period,Nb_pts,0).create_for_plane_undulator()
print(type(T))
#T.draw()

Rad= UndulatorRadiationFactory(T,1).create_for_single_electron(distance=30.0)
print(type(Rad))
#Rad.draw()

Rad2= UndulatorRadiationFactory(T,3).create_for_single_electron(distance=30.0)
print(type(Rad2))
#Rad2.draw()

D=np.linspace(L*0.5+0.1,20.0,40)
error=Rad.compare_with_ditance(Rad2,D)
print(type(error))
print(error.shape)
print(error)
plt.plot(D,error)
plt.show()


