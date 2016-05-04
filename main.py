import numpy as np
import matplotlib.pyplot as plt
from pylab import *
from mpl_toolkits.mplot3d import Axes3D
import scipy.constants as codata
import scipy.integrate as integrate


from calc_undulator_2 import radiation, trajectory_undulator_reference, trajectory_undulator

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
Nb_period = 10
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
omega_u = Beta_et *codata.c * 2.0 * np.pi / lambda_u
a=(lambda_u / (codata.c * Beta_et)) * (Nb_period / 2)
t = np.linspace(-a, a, N)
By=Bo*np.sin(omega_u*t)
Vo=( codata.e*Beta_et*Bo/(gamma*codata.m_e*omega_u))*np.cos((-1.0)*Nb_period*np.pi)
Zo = (-lambda_u*Nb_period)/(2.0*codata.c )
Xo=0.0



T = trajectory_undulator_reference(K=K,gamma=gamma, lambda_u=lambda_u, Nb_period=Nb_period, Nb_point=Nb_pts, Beta_et=Beta_et)
T2 = trajectory_undulator(By=By,t=t,gamma=gamma, Beta=Beta, Beta_et=Beta_et,Vo=Vo , Xo=Xo,Zo=Zo)
print("number of points : %d" % (T.shape[1]))
print("number of points : %d" % (T2.shape[1]))

# # #
# # # #


# plt.plot(T[0],T[1])
# #plt.plot(T[0],T2[1])
# plt.title(" X = f(t) ")
# plt.xlabel('t')
# plt.ylabel('x')
# plt.show()
#
# test1=Beta_et*T[0]
# plt.plot(T[0],T[3])
# #plt.plot(T[0],T2[3])
# #plt.plot(T[0],test1)
# plt.title(" Z = f(t) ")
# plt.xlabel('t')
# plt.ylabel('Z')
# plt.show()
#
# plt.plot(T[0],T[4])
# #plt.plot(T[0],T2[4])
# plt.title(" Vx = f(t) ")
# plt.xlabel('t')
# plt.ylabel('Vx')
# plt.show()
#
# test2=Beta_et*np.ones(N)
# plt.plot(T[0],T[6])
# #plt.plot(T[0],T2[6])
# plt.plot(T[0],test2)from calc_undulator_2 import radiation

# plt.title(" Vz = f(t) ")
# plt.xlabel('t')
# plt.ylabel('Vz')
# plt.show()
#
#
#
# plt.plot(T[0], T[7])
# #plt.plot(T[0], T2[7])
# plt.title(" Ax = f(t) ")
# plt.xlabel('t')
# plt.ylabel('Ax')
# plt.show()
#
# plt.plot(T[0], T[9])
# #plt.plot(T[0], T2[9])
# plt.title(" Az = f(t) ")
# plt.xlabel('t')
# plt.ylabel('Az')
# plt.show()

#
# 3D plot
fig = figure()
ax = Axes3D(fig)
X=np.arange(0.0, 0.00101, 0.00001)
Y=np.arange(0.0, 0.00101, 0.00001)
Z2 = radiation(K=K,E=E,trajectory=T,X=X,Y=Y)
X,Y = np.meshgrid(X,Y)
print('debut du plot')
ax.plot_surface(X,Y, Z2, rstride=1, cstride=1)
ax.set_xlabel("X")
ax.set_ylabel('Y')
ax.set_zlabel("flux")
show()
