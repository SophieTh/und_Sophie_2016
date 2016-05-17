import numpy as np
import matplotlib.pyplot as plt
from pylab import *
from mpl_toolkits.mplot3d import Axes3D
import scipy.constants as codata
import scipy.integrate as integrate
from scipy.interpolate import interp1d


from calc_undulator_2 import radiation, draw_trajectory,undulator_trajectory
from fct_manuel import window_ftr


print("c: %g m/s**2"%codata.c)
print("e: %g C"%codata.e)
print("h: %g J s"%codata.h)
print("hbar: %g J s"%codata.hbar)
print("epsilon_0: %g F/m"%codata.epsilon_0)
print("m_e: %g kg"%codata.m_e)

######################################################


K = 1.87
E = 1.3e9
lambda_u = 0.035
Nb_period = 12
Nb_pts = 100


gamma = E /0.511e6
print('gamma =')
print(gamma)
Beta = np.sqrt(1.0 - (1.0 / gamma ** 2))
print('Beta = ')
print(Beta)
#Beta_et = 1.0-(1.0/(2.0*gamma**2))*(1.0+(K**2)/2.0)
Beta_et=Beta*(1.0-(K/(2.0*gamma))**2)
print('Beta* = ')
print(Beta_et)
Bo = K / (93.4 * lambda_u)
print('Bo = ')
print(Bo)
N = Nb_period * Nb_pts + 1
print(N)
omega1 = ((2.0 * gamma ** 2) / (1.0 + (K ** 2) / 2.0)) * ((2.0 * np.pi * codata.c) / lambda_u)
print('1st harmonic')
print(omega1)

##






# recuperation des donnees de B en array en fonction de z
reference=np.load("x_ray_booklet_field.npz")
reference.keys()
Z=reference['ct']
Z -= (Z[len(Z)-1])/2.0
By2=reference['B_y']
nb_enlarg=100

Z_By=np.zeros((2,len(Z)))
Z_By[0]=Z
Z_By[1]=By2

T=undulator_trajectory(K,E,lambda_u,Nb_period,Nb_pts,Z_By,type_trajectory=2)
#
#
print('ok')

draw_trajectory(T)
# #
####3D plot
fig = figure()
ax = Axes3D(fig)
X=np.arange(0.0, 0.00101, 0.00001)
Y=np.arange(0.0, 0.00101, 0.00001)
X,Y = np.meshgrid(X,Y)
Z2 = radiation(K=K,E=E,trajectory=T,X=X,Y=Y)
ax.plot_surface(X,Y, Z2, rstride=1, cstride=1)
ax.set_xlabel("X")
ax.set_ylabel('Y')
ax.set_zlabel("flux")
show()
