import numpy as np
import matplotlib.pyplot as plt
from pylab import *
from mpl_toolkits.mplot3d import Axes3D
import scipy.constants as codata
import scipy.integrate as integrate
from scipy.interpolate import interp1d


from calc_undulator_2 import draw_trajectory,undulator_trajectory, draw_2_trajectory,radiation_single_electron_2D


def erreur_trajectory( ref=np.zeros((10,101)),test=np.zeros((10,101))) :
    # print(' x max=')
    # print(abs(ref[1]).max())
    # print(' z max=')
    # print(abs(ref[3]).max())
    # print(' Vx max=')
    # print(abs(ref[4]).max())
    # print(' Vz max=')
    # print(abs(ref[6]).max())
    # print(' Ax max=')
    # print(abs(ref[7]).max())
    # print(' Az max=')
    # print(abs(ref[9]).max())

    print("erreur temps")
    egalite_temp = (T0[0] - T1[0])
    print(egalite_temp.max())
    erreur01 = np.zeros((10, N))
    erreur01[0] = ref[0]
    k = 1
    while k < 10:
        erreur01[k]+= abs(test[k]- ref[k])
        k += 1

    # print(' erreur x max=')
    # print(abs(erreur01[1]).max())
    # print(' erreur z max=')
    # print(abs(erreur01[3]).max())
    # print('erreur Vx max=')
    # print(abs(erreur01[4]).max())
    # print(' erreur Vz max=')
    # print(abs(erreur01[6]).max())
    # print('erreur Ax max=')
    # print(abs(erreur01[7]).max())
    # print('erreur Az max=')
    # print(abs(erreur01[9]).max())
    print("proportion erreur.max / ref.max : ")
    print(' erreur x max=')
    print((abs(erreur01[1])).max() / (abs(T0[1])).max())
    print(' erreur z max=')
    print((abs(erreur01[3])).max() / (abs(T0[3])).max())
    print('erreur Vx max=')
    print((abs(erreur01[4])).max() / (abs(T0[4])).max())
    print(' erreur Vz max=')
    print((abs(erreur01[4])).max() / (abs(T0[4])).max())
    print('erreur Ax max=')
    print((abs(erreur01[6])).max() / (abs(T0[6])).max())
    print('erreur Az max=')
    print((abs(erreur01[7])).max() / (abs(T0[7])).max())

    return erreur01

def erreur_of_radiation_for_trajectory(X=np.arange(0.0, 0.0301, 0.0003),Y=np.zeros(101),
                                       ref=np.zeros((10,101)),test=np.zeros((10,101))) :
    Z0, maxref = radiation_single_electron_2D(K=K, E=E, trajectory=ref, X=X, Y=Y, D=30.0)
    Z1, maxtest = radiation_single_electron_2D(K=K, E=E, trajectory=test, X=X, Y=Y, D=30.0)

    u=np.linspace(0.0,len(X),len(X))/len(X)
    print(u)
    plt.plot(u, Z0)
    plt.plot(u, Z1)
    plt.title(" ref and test = f(u) ")
    plt.show()

    plt.plot(u, (abs(Z1 - Z0) / maxref))
    plt.title(" relativ error of Z1 for Z0 ")
    plt.ylabel('|Z1-Z0|/max0')
    plt.xlabel('u')
    plt.show()

    return abs(Z1-Z0)



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
Beta_et = 1.0-(1.0/(2.0*gamma**2))*(1.0+(K**2)/2.0)
#Beta_et=Beta*(1.0-(K/(2.0*gamma))**2)
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

# a magnetic field
Z=np.linspace(-(lambda_u) * (Nb_period / 2), (lambda_u ) * (Nb_period / 2), N)
By = -Bo * np.sin((2.0 * np.pi / lambda_u) *Z)

Z_By=np.zeros((2,len(Z)))
Z_By[0]=Z
Z_By[1]=By

T0=undulator_trajectory(K,E,lambda_u,Nb_period,Nb_pts,Z_By,type_trajectory=0)

T1=undulator_trajectory(K,E,lambda_u,Nb_period,Nb_pts,Z_By,type_trajectory=1,
                        Vx=T0[4][0]*codata.c,Vz=T0[6][0]*codata.c,Xo=T0[1][0]*codata.c)
T2=undulator_trajectory(K,E,lambda_u,Nb_period,Nb_pts,Z_By,type_trajectory=2,
                        Vx=T0[4][0] * codata.c, Vz=T0[6][0] * codata.c, Xo=T0[1][0] * codata.c)


#draw_2_trajectory(T0,T1)
erreur01 =erreur_trajectory(ref=T0,test=T1)
# draw_trajectory(erreur01)


#draw_2_trajectory(T0,T2)
erreur02 =erreur_trajectory(ref=T0,test=T2)
#draw_trajectory(erreur02)


#draw_2_trajectory(T0,T2)
erreur12 =erreur_trajectory(ref=T1,test=T2)
#draw_trajectory(erreur02)


X=np.arange(0.0, 0.0301, 0.0003)
Y=np.zeros_like(X)

erreur_rad_01=erreur_of_radiation_for_trajectory(X=X,Y=Y,ref=T0,test=T1)