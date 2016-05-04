import numpy as np
import matplotlib.pyplot as plt
from pylab import *
from mpl_toolkits.mplot3d import Axes3D
import scipy.constants as codata


print("c: %g m/s**2"%codata.c)
print("e: %g C"%codata.e)
print("h: %g J s"%codata.h)
print("hbar: %g J s"%codata.hbar)
print("epsilon_0: %g F/m"%codata.epsilon_0)
print("m_e: %g kg"%codata.m_e)


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
    trajectory[3] = Beta_et*trajectory[0] - ((K/gamma)**2) * (1.0/(8.0*ku*codata.c))*np.sin(2.0*omega_u * trajectory[0])
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



#  hypothesis :'FAR FIELD'
def electric_field_undulator(omega1=2.53465927101*10**17, trajectory=np.zeros((11, 10)), x=0.0, y=0.0, D=30.0):
    N = trajectory.shape[1]
    # in meters :
    # R = np.sqrt(x ** 2 + y ** 2 + D ** 2)
    # n_chap = np.array([x, y, D]) / R
    #in radian :
    n_chap = np.array([x,y,1.0-0.5*(x**2 + y**2)])
    trajectory = np.transpose(trajectory)
    integrand = np.full((N, 3), 0. + 1j * 0., dtype=np.complex)
    E = np.full((3,), 0. + 1j * 0., dtype=np.complex)
    for i in range(N):
        integrand[i] = (np.cross(n_chap, np.cross(n_chap - (trajectory[i][4:7]), (trajectory[i][7:10]))) + np.array([0. + 1j * 0., 0. + 1j * 0., 0. + 1j * 0.]))
        integrand[i] *= (1.0 / (1.0 - np.dot(n_chap, trajectory[i][4:7]))) ** 2
        integrand[i] *= np.exp(0. + 1j * omega1 * (trajectory[i][0] - np.dot(n_chap, trajectory[i][1:4])))

    integrand = np.transpose(integrand)
    trajectory=np.transpose(trajectory)
    for k in range(3):
        E[k] = np.trapz( integrand[k],trajectory[0])
    return E


def dI_dOmega(K=1.87 , E = 1.3*10**9 , lambda_u= 0.035 , trajectory=np.zeros((11,10)) , x=0.00 , y= 0.00 ,D=30.0 ) :
    gamma = E /0.511e6
    omega1 = ((2.0 * gamma ** 2) / (1.0 + (K ** 2) / 2.0)) * ((2.0 * np.pi * codata.c) / lambda_u)
    Elec= electric_field_undulator( omega1=omega1 ,trajectory=trajectory,x=x, y=y , D=D )
    norm_2 = (np.abs(Elec[0]) ** 2 + np.abs(Elec[1])** 2 + np.abs(Elec[2])** 2)
    c1 = codata.e ** 2 * omega1 ** 2 / (16 * np.pi ** 3 * codata.epsilon_0 * codata.c )
    c2 = 1.0 / codata.e  # multiply by number of electrons in 1 A
    c3 = 2.0*np.pi / (codata.h * omega1)  # divide by e energy (to get number of photons)
    c4 = 1e-2 / omega1  # to get 1% energy bandwidth  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    c5 = 0.1e-3 * 0.1e-3  # from rad to .1mrad angle bandwidth
    c6= codata.e*1e-10/(8.0*np.pi**2*codata.epsilon_0*codata.c*codata.h)
    return c6*norm_2


def Radiation(K=1.87, E=1.3 *10**9 , lambda_u=0.035, trajectory=np.zeros((11, 10)),D=30.0 ,X=np.arange(-0.0011, 0.0011, 0.00002), Y=np.arange(-0.0011, 0.0011, 0.00002)):
    res = np.zeros((len(X), len(Y)))
    print(' dimention de la matrice dna la fct MAG')
    print(res.shape)
    for i in range(len(X)):
        print(i)
        for k in range(len(Y)):
            res[i][k] = dI_dOmega(K=K, E=E, lambda_u=lambda_u, trajectory=trajectory,x=X[i] , y=Y[k],D=D)
    print(res.max())
    return res


#  hypothesis :'FAR FIELD'
def electric_field_undulator_2(f=4.06e17, trajectory=np.zeros((11, 10)), x=0.0, y=0.0, D=30.0):
    N = trajectory.shape[1]
    # in meters :
    # R = np.sqrt(x ** 2 + y ** 2 + D ** 2)
    # n_chap = np.array([x, y, D]) / R
    #in radian :
    n_chap = np.array([x,y,1.0-0.5*(x**2 + y**2)])
    trajectory = np.transpose(trajectory)
    integrand = np.full((N, 3), 0. + 1j * 0., dtype=np.complex)
    E = np.full((3,), 0. + 1j * 0., dtype=np.complex)
    for i in range(N):
        integrand[i] = (np.cross(n_chap, np.cross(n_chap - (trajectory[i][4:7]), (trajectory[i][7:10]))) + np.array([0. + 1j * 0., 0. + 1j * 0., 0. + 1j * 0.]))
        integrand[i] *= (1.0 / (1.0 - np.dot(n_chap, trajectory[i][4:7]))) ** 2
        integrand[i] *= np.exp(0. + 1j * 2.0*np.pi*f* (trajectory[i][0] - np.dot(n_chap, trajectory[i][1:4])))

    integrand = np.transpose(integrand)
    trajectory=np.transpose(trajectory)
    for k in range(3):
        E[k] = np.trapz( integrand[k],trajectory[0])
    return E


def dI_dOmega_2(K=1.87 , E = 1.3*10**9 , lambda_u= 0.035 , trajectory=np.zeros((11,10)) , x=0.00 , y= 0.00 ,D=30.0 ) :
    gamma = E / mc2
    Beta_et = 1.0 - (1.0 / (2.0 * gamma ** 2)) * (1.0 + (K ** 2) / 2.0)
    #f = ((2.0 * gamma ** 2) / (1.0 + (K ** 2) / 2.0)) * c / lambda_u
    f=2.0*Beta_et*codata.c/((1.0-Beta_et)*lambda_u)
    Elec= electric_field_undulator_2( f=f ,trajectory=trajectory,x=x, y=y , D=D )
    norm_2 = (np.abs(Elec[0]) ** 2 + np.abs(Elec[1])** 2 + np.abs(Elec[2])** 2)
    c1=e**2/(2.0*np.pi*codata.c)
    c2 = 1e3/(f*codata.h*2.0*np.pi*f*codata.e)
    return c1*norm_2


def Radiation_2(K=1.87, E=1.3 *10**9 , lambda_u=0.035, trajectory=np.zeros((11, 10)),D=30.0 ,X=np.arange(-0.0011, 0.0011, 0.00002), Y=np.arange(-0.0011, 0.0011, 0.00002)):
    res = np.zeros((len(X), len(Y)))
    print(' dimention de la matrice dna la fct MAG')
    print(res.shape)
    for i in range(len(X)):
        print(i)
        for k in range(len(Y)):
            res[i][k] = dI_dOmega_2(K=K, E=E, lambda_u=lambda_u, trajectory=trajectory,x=X[i] , y=Y[k],D=D)
    print(res.max())
    return res



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
t = np.linspace(-(lambda_u / (codata.c * Beta_et)) * (Nb_period / 2), (lambda_u / (codata.c * Beta_et)) * (Nb_period / 2), N)
By=Bo*np.sin(omega_u*t)
Vo=( codata.e*Beta_et*Bo/(gamma*codata.m_e*omega_u))*np.cos((-1.0)*Nb_period*np.pi)
Zo = (-lambda_u*Nb_period)/(2.0*codata.c )
Xo=0.0



T = trajectory_undulator_reference(K=K,gamma=gamma, lambda_u=lambda_u, Nb_period=Nb_period, Nb_point=Nb_pts, Beta_et=Beta_et)
T2 = trajectory_undulator(By=By,t=t,gamma=gamma, Beta=Beta, Beta_et=Beta_et,Vo=Vo , Xo=Xo,Zo=Zo)
print("number of points : %d" % (T.shape[1]))
print(T.shape)
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
# plt.plot(T[0],test2)
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
Z2 = Radiation(K=K,E=E,trajectory=T,X=X,Y=Y)
X,Y = np.meshgrid(X,Y)
print('debut du plot')
ax.plot_surface(X,Y, Z2, rstride=1, cstride=1)
ax.set_xlabel("X")
ax.set_ylabel('Y')
ax.set_zlabel("flux")
show()
print(e**2/(2.0*np.pi*c))