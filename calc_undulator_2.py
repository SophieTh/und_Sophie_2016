import numpy as np
import matplotlib.pyplot as plt
from pylab import *
from mpl_toolkits.mplot3d import Axes3D
import scipy.constants as codata
import scipy.integrate as integrate
from scipy.integrate import ode
from scipy.integrate import odeint
from scipy.interpolate import interp1d


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
    # trajectory *= codata.c
    # trajectory[0] *= (1.0/codata.c)
    return trajectory



# hypothesis : B=(0,By,0) and norm(v)=constant
def trajectory_undulator( By=np.zeros(101) ,gamma=2544.03131115, Beta=0.99996, Beta_et=0.99993,Vo=0.0 , Xo=0.0):
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
    to=(lambda_u / (codata.c * Beta_et)) * (Nb_period / 2)
    trajectory[0] =np.linspace(-to, to, N)

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
        trajectory[3][i] = np.trapz(trajectory[6][0:(i + 1)], trajectory[0][0:(i + 1)]) -lambda_u * (Nb_period / 2)
        # trajectory[4][i] = integrate.simps(trajectory[7][0:(i + 1)], trajectory[0][0:(i + 1)])
     #Az
    trajectory[9]=-(trajectory[7]*trajectory[4])/trajectory[6]

    return trajectory


def f3(y,t,cst,B) :
    return [ -cst*B( y[5]) * y[2],
            0.0,
             cst * B(y[5]) * y[0],
            y[0],
            0.0,
            y[2]]

def enlargement_vector_for_interpolation(Z,By,nb_elarg) :
    dz = Z[1] - Z[0]
    elarg_1 = np.linspace(Z[0] - nb_elarg * dz, Z[0] - dz, nb_elarg)
    elarg_2 = np.linspace(Z[len(Z) - 1] + dz, Z[len(Z) - 1] + nb_elarg * dz, nb_elarg)
    Z = np.concatenate((elarg_1, Z))
    Z = np.concatenate((Z, elarg_2))
    elarg_3 = np.zeros(nb_elarg)
    By = np.concatenate((elarg_3, By))
    By = np.concatenate((By, elarg_3))
    return Z,By




def trajectory_undulator2 (By,Z ,nb_enlarg,Beta,gamma) :
    N=len(Z)
    #   trajectory =
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
    trajectory[0] = Z/codata.c
    Z_centred=Z_centre= Z-(Z[len(Z)-1]/2.0)
    Vo = [0.0, 0.0, Beta * codata.c, 0.0, 0.0, Z_centred[0]]
    Z_centred,By=enlargement_vector_for_interpolation(Z_centred, By, nb_enlarg)
    B = interp1d(Z_centred, By)
    cst=-codata.e/(codata.m_e*gamma)
    res = odeint(f3,Vo,trajectory[0],args=(cst,B))
    print('dim de res')
    print(res.shape)
    res = np.transpose(res)
    trajectory[4] = res[0]
    trajectory[5] = res[1]
    trajectory[6] = res[2]
    trajectory[1] = res[3]
    trajectory[2] = res[4]
    trajectory[3] = res[5]
    trajectory[7] = -cst*B(trajectory[3])*trajectory[6]
    trajectory[9] = cst* B(trajectory[3])*trajectory[4]
    k=1
    while k<10 :
        trajectory[k] *= 1.0/codata.c
        k+=1
    return trajectory


########################### essai
#
# # hypothesis : B=(0,By,0) and norm(v)=constant
def trajectory_undulator3( By=np.zeros(101) ,z=np.zeros(101) ,Nb_period=10,lambda_u=0.035,gamma=2544.03131115, Beta=0.99996, Beta_et=0.99993,Vo=0.0 , Xo=0.0,E=1.3e9):
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
    # Z
    trajectory[3] = z

    # # ax
    # #trajectory[7] = codata.e/(gamma*codata.m_e)*By
    # trajectory[7] = 0.3*codata.c /(E) * By
    # to=-(lambda_u / (codata.c * Beta_et)) * (Nb_period / 2)
    #
    # # dX/dz
    # for i in range(N):
    #     trajectory[4][i] = np.trapz(trajectory[7][0:(i + 1)], trajectory[3][0:(i + 1)])+Vo
    #     # trajectory[3][i] = integrate.simps(trajectory[5][0:(i + 1)], trajectory[1][0:(i + 1)])
    # #trajectory[4]*= (1.0/(Beta_et*codata.c))**2
    # Vx_2=np.sqrt(1.0*trajectory[4]**2)
    # L = np.zeros(N)
    # for i in range(N):
    #     L[i] = np.trapz(Vx_2[0:(i + 1)], trajectory[3][0:(i + 1)])
    # trajectory[0] = to + L / (codata.c * Beta)
    #
    # # dZ/dt
    # trajectory[6][0]=Beta_et*codata.c
    # trajectory[6][N-1]=Beta_et*codata.c
    # k=0
    # while k < N-1:
    #     trajectory[6][k] = (trajectory[3][k+1]-trajectory[3][k])/(trajectory[0][k+1]-trajectory[0][k])
    #     k +=1

    return trajectory


def trajectory_undulator4 (By ,Beta_et,Beta,lambda_u,Nb_period,Bo,gamma,E) :
    N= len(By)
    #   trajectory =
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
    #Vo=[codata.c*np.sqrt(Beta**2-Beta_et**2),0.0,Beta_et*codata.c,0.0,0.0,-lambda_u*(Nb_period/2.0)]
    Vo=[0.0,0.0,Beta*codata.c,0.0,0.0,-lambda_u*(Nb_period/2.0)]
    #Vo = [-0.3*codata.c*lambda_u/(E*2.0*np.pi), 0.0, codata.c*np.sqrt(Beta**2 - Beta_et**2), 0.0, 0.0, -lambda_u * (Nb_period / 2.0)]
    #res = ode(f).set_integrator('dopri5')
    cst=-Bo*codata.e/(codata.m_e*gamma)
    print(cst)
    #cst = 1e1
    #res = ode(f).set_integrator('vode', method='bdf').set_initial_value(Vo, trajectory[0][0]).set_f_params(cst)
    omega=2.0*np.pi/lambda_u
    res = ode(f).set_integrator('dopri5').set_initial_value(Vo, 0.0).set_f_params(cst,omega)
    i=0
    dt=1.0/(N-1)
    while i < N :
        #res.set_f_params(By[i])
        #print(res.t, res.integrate(res.t + dt))
        trajectory[0][i] = res.t
        trajectory[4][i] = res.integrate(res.t + dt)[0]
        trajectory[5][i] = res.integrate(res.t + dt)[1]
        trajectory[6][i] = res.integrate(res.t + dt)[2]
        trajectory[1][i] = res.integrate(res.t + dt)[3]
        trajectory[2][i] = res.integrate(res.t + dt)[4]
        trajectory[3][i] = res.integrate(res.t + dt)[5]
        i +=1
    trajectory[7]=cst*np.sin(omega*trajectory[3])*trajectory[6]
    trajectory[9] = -cst * np.sin(omega * trajectory[3]) * trajectory[4]
    k=1
    while k<10 :
        trajectory[k] *= 1.0/codata.c
        k+=1
    return trajectory

#################################################


def energy_radiated(omega1=2.53465927101*10**17,trajectory=np.zeros((11,10)) , x=0.00 , y=0.0,D=30.0):
    N = trajectory.shape[1]
    # in meters :
    # R = np.sqrt(x ** 2 + y ** 2 + D ** 2)
    # n_chap = np.array([x, y, D]) / R
    # in radian :
    n_chap = np.array([x, y, 1.0 - 0.5 * (x ** 2 + y ** 2)])
    E = np.full((3,), 0. + 1j * 0., dtype=np.complex)
    integrand = np.full((3,N), 0. + 1j * 0., dtype=np.complex)
    Alpha=trajectory[7]*(n_chap[2]-trajectory[6]) - (n_chap[0]-trajectory[4])*trajectory[9]
    Alpha2=np.exp(0. + 1j * omega1 * (trajectory[0] - n_chap[0]*trajectory[1]-n_chap[2]*trajectory[3]))
    integrand[0] += (-(n_chap[1]**2)*trajectory[7]-n_chap[2]*Alpha)*Alpha2
    integrand[1] += n_chap[1]*(n_chap[0]*trajectory[7]+n_chap[2]*trajectory[9])*Alpha2
    integrand[2] += (-(n_chap[1]**2)*trajectory[9]+n_chap[0]*Alpha)*Alpha2
    integrand *= (1.0 / (1.0 - n_chap[0]*trajectory[4]-n_chap[2]*trajectory[6])) ** 2
    for k in range(3):
        E[k] = np.trapz(integrand[k], trajectory[0])
        #E[k] = integrate.simps(integrand[k], trajectory[0])
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
        #print(i)
        for k in range(len(Y)):
            res[i][k] = c6*energy_radiated(omega1=omega1,trajectory=trajectory , x=X[i] , y=Y[k] )
    print(res.max())
    return res