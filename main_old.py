import numpy as np
import matplotlib.pyplot as plt
from pylab import *
from mpl_toolkits.mplot3d import Axes3D
import scipy.constants as codata
import scipy.integrate as integrate
from scipy.interpolate import interp1d


from calc_undulator_2 import radiation_single_electron, draw_trajectory,undulator_trajectory,radiation_single_electron2



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
Nb_pts = 20


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

print(lambda_u/(Beta_et*codata.c))


# recuperation des donnees de B en array en fonction de z
reference=np.load("x_ray_booklet_field.npz")
print(reference.keys())
Z=reference['ct']
Z -= (Z[len(Z)-1])/2.0
By2=reference['B_y']

Z_By=np.zeros((2,len(Z)))
Z_By[0]=Z
Z_By[1]=By2

plt.plot(Z,By2)
plt.title(" SRW ")
plt.show()

#attention changer les entree
T=undulator_trajectory(K,E,lambda_u,Nb_period,Nb_pts,type_trajectory=2,Vx=0.0,Xo=0.0,Vz=Beta*codata.c)
#
#draw_trajectory(T)
print('ok')
#print(T[6]-Beta_et)


# fig = figure()
# ax = Axes3D(fig)
# X=np.arange(0.0, 0.00301, 0.00003)
# Y=np.arange(0.0, 0.00301, 0.00003)
# print(X.shape)
# print(len(X.shape))
# X,Y = np.meshgrid(X,Y)
# Z = radiation_single_electron(K=K,E=E,trajectory=T,X=X,Y=Y,D=3.0)
# print('plot')
# print(X.shape)
# print(Y.shape)
# ax.plot_surface(X,Y, Z, rstride=1, cstride=1)
# ax.set_xlabel("X")
# ax.set_ylabel('Y')
# ax.set_zlabel("flux")
# show()


