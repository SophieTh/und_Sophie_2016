import numpy as np
import matplotlib.pyplot as plt
from pylab import *
from mpl_toolkits.mplot3d import Axes3D
import scipy.constants as codata
import scipy.integrate as integrate

def trajectory_undulator_reference(K=1.87 , gamma=2544.03131115, lambda_u=0.020, Nb_period=10, Nb_point=10, Beta_et=0.99993):
    N = Nb_period * Nb_point + 1
    ku= 2.0*np.pi/lambda_u
    omega_u = Beta_et*codata.c *ku
    # trajectory =
    #   [t........]
    # 	[ X/c......]
    # 	[ Y/c ......]
    #	[ Z/c ......]
    # 	[ Vx/c .....]
    # 	[ Vy/c .....]
    # 	[ Vz/c .....]
    # 	[ Ax/c .....]
    # 	[ Ay/c .....]
    # 	[ Az/c .....]
    trajectory = np.zeros((10, N))
    # t
    trajectory[0] = np.linspace(-(lambda_u / (codata.c * Beta_et)) * (Nb_period / 2), (lambda_u / (codata.c * Beta_et)) * (Nb_period / 2), N)
    #trajectory[0] = np.linspace(0.0,(lambda_u / (c * Beta_et)) * (Nb_period), N)
    # X et Z en fct de t
    trajectory[3] = Beta_et*trajectory[0] - ((K/gamma)**2) * (1.0/(8.0*ku*codata.c))* np.sin(2.0*omega_u * trajectory[0])
    trajectory[1] = (K/(gamma*ku*codata.c))* np.sin(omega_u * trajectory[0])
    # Vx et Vz en fct de t
    trajectory[6] = Beta_et - ((K/gamma)**2) * (2.0*omega_u/(8.0*ku*codata.c ))*np.cos(2.0*omega_u * trajectory[0])
    trajectory[4] = (K/(gamma*ku*codata.c))*omega_u* np.cos(omega_u * trajectory[0])
    # Ax et Az en fct de t
    trajectory[9] = ((2.0*omega_u*K/gamma)**2) * (1.0/(8.0*ku*codata.c))*np.sin(2.0*omega_u * trajectory[0])
    trajectory[7] = -(K/(gamma*ku*codata.c))*(omega_u**2)* np.sin(omega_u * trajectory[0])
    return trajectory



# hypothesis : B=(0,By,0) and norm(v)=constant
def trajectory_undulator( By=np.zeros(101) ,t=np.zeros(101) ,gamma=2544.03131115, Beta=0.99996, Beta_et=0.99993,Vo=0.0 , Xo=0.0 , Zo=0.0):
    N= len(By)
    # trajectory =
    #   [t........]
    # 	[ X/c......]
    # 	[ Y/c ......]
    #	[ Z/c ......]
    # 	[ Vx/c .....]
    # 	[ Vy/c .....]
    # 	[ Vz/c .....]
    # 	[ Ax/c .....]
    # 	[ Ay/c .....]
    # 	[ Az/c .....]
    trajectory = np.zeros((10,N))
    # t
    trajectory[0] = t
    # Ax(t)
    Xm= codata.e*Beta_et/(gamma*codata.m_e)
    trajectory[7] = -Xm * By
    # Vx et Vz
    for i in range(N):
        trajectory[4][i] = np.trapz(trajectory[7][0:(i + 1)], trajectory[0][0:(i + 1)]) + Vo
        # trajectory[4][i] = integrate.simps(trajectory[7][0:(i + 1)], trajectory[0][0:(i + 1)])
    trajectory[6] = np.sqrt((Beta)**2 - trajectory[4]**2)
    # X et Z
    for i in range(N):
        trajectory[1][i] = np.trapz(trajectory[4][0:(i + 1)], trajectory[0][0:(i + 1)])+Xo
        # trajectory[4][i] = integrate.simps(trajectory[7][0:(i + 1)], trajectory[0][0:(i + 1)])
        trajectory[3][i] = np.trapz(trajectory[6][0:(i + 1)], trajectory[0][0:(i + 1)]) +Zo
        # trajectory[4][i] = integrate.simps(trajectory[7][0:(i + 1)], trajectory[0][0:(i + 1)])
     #Az
    trajectory[9]=-(trajectory[7]*trajectory[4])/trajectory[6]
    #trajectory=np.transpose(trajectory)
    return trajectory


def energy_radiated(omega1=2.53465927101*10**17,trajectory=np.zeros((11,10)) , x=0.00 , y=0.0):
    N = trajectory.shape[1]
    # in meters :
    # R = np.sqrt(x ** 2 + y ** 2 + D ** 2)
    # n_chap = np.array([x, y, D]) / R
    # in radian :
    n_chap = np.array([x, y, 1.0 - 0.5 * (x ** 2 + y ** 2)])
    E = np.full((3,), 0. + 1j * 0., dtype=np.complex)
    integrand = np.full((3,N), 0. + 1j * 0., dtype=np.complex)
    Alpha=trajectory[7]*(n_chap[2]-trajectory[6])-(n_chap[0]-trajectory[4])*trajectory[9]
    Alpha2=np.exp(0. + 1j * omega1 * (trajectory[0] - n_chap[0]*trajectory[1]-n_chap[2]*trajectory[3]))
    integrand[0] += ((-n_chap[1]**2)*trajectory[7]-n_chap[2]*Alpha)*Alpha2
    integrand[1] += n_chap[1]*(n_chap[0]*trajectory[7]+n_chap[2]*trajectory[9])*Alpha2
    integrand[2] += (-n_chap[1]**2*trajectory[9]+n_chap[0]*Alpha)*Alpha2
    integrand *= (1.0 / (1.0 - n_chap[0]*trajectory[4]-n_chap[2]*trajectory[6])) ** 2
    for k in range(3):
        E[k] = np.trapz(integrand[k], trajectory[0])
        # E[k] = integrate.simps(integrand[k], trajectory[0])
    # np.linalg.norm
    return (np.abs(E[0]) ** 2 + np.abs(E[1])** 2 + np.abs(E[2])** 2)



def radiation(K=1.87, E=1.3 * 10 ** 9, lambda_u=0.035, trajectory=np.zeros((11, 10)), D=30.0, X=np.arange(-0.0011, 0.0011, 0.00002), Y=np.arange(-0.0011, 0.0011, 0.00002)):
    res = np.zeros((len(X), len(Y)))
    gamma = E / 0.511e6
    omega1 = ((2.0 * gamma ** 2) / (1.0 + (K ** 2) / 2.0)) * ((2.0 * np.pi * codata.c) / lambda_u)
    # c1 = codata.e ** 2 * omega1 ** 2 / (16 * np.pi ** 3 * codata.epsilon_0 * codata.c )
    # c2 = 1.0 / codata.e  # multiply by number of electrons in 1 A
    # c3 = 2.0*np.pi / (codata.h * omega1)  # divide by e energy (to get number of photons)
    # c4 = 1e-2 / omega1  # to get 1% energy bandwidth  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # c5 = 0.1e-3 * 0.1e-3  # from rad to .1mrad angle bandwidth
    c6= codata.e*1e-10/(8.0*np.pi**2*codata.epsilon_0*codata.c*codata.h)
    for i in range(len(X)):
        for k in range(len(Y)):
            res[i][k] = c6*energy_radiated(omega1=omega1,trajectory=trajectory , x=X[i] , y=Y[k] )
    print(res.max())
    return res