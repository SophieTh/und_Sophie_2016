import numpy as np
import matplotlib.pyplot as plt
from pylab import *
from mpl_toolkits.mplot3d import Axes3D

# light speed
c = 2.99792458 * 10 ** 8

# elementary charge
e = 1.602176565 / 10 ** (19)

# masse electron * c2
mc2 = 511 * 10 ** 3

# ...
epsilon0 = 8.85418782 / 10 ** (12)

#masse electron
me= 9.109/10**(31)


# plank
# h=6,62607004/10**(34)  # m*m*kg/s



def trajectory_undulator_reference(K=1.87 , gamma=2544.03131115, lambda_u=0.020, Nb_period=10, Nb_point=10, Beta_et=0.99993):
    N = Nb_period * Nb_point + 1
    ku= 2.0*np.pi/lambda_u
    omega_u = Beta_et*c*ku
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
    trajectory[0] = np.linspace(-(lambda_u / (c * Beta_et)) * (Nb_period / 2), (lambda_u / (c * Beta_et)) * (Nb_period / 2), N)
    #trajectory[0] = np.linspace(0.0,(lambda_u / (c * Beta_et)) * (Nb_period), N)
    # X et Z en fct de t
    trajectory[3] = Beta_et*trajectory[0] - ((K/gamma)**2) * (1.0/(8.0*ku*c))*np.sin(2.0*omega_u * trajectory[0])
    trajectory[1] = (K/(gamma*ku*c))* np.sin(omega_u * trajectory[0])
    # Vx et Vz en fct de t
    trajectory[6] = Beta_et - ((K/gamma)**2) * (2.0*omega_u/(8.0*ku*c))*np.cos(2.0*omega_u * trajectory[0])
    trajectory[4] = (K/(gamma*ku*c))*omega_u* np.cos(omega_u * trajectory[0])
    # Ax et Az en fct de t
    trajectory[9] = ((2.0*omega_u*K/gamma)**2) * (1.0/(8.0*ku*c))*np.sin(2.0*omega_u * trajectory[0])
    trajectory[7] = -(K/(gamma*ku*c))*(omega_u**2)* np.sin(omega_u * trajectory[0])
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
    trajectory = np.zeros((10, N))
    # t
    trajectory[0] = t
    # Ax(t)
    Xm= e*Beta_et/(gamma*me)
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
    trajectory=np.transpose(trajectory)
    return trajectory






#  hypothesis :'FAR FIELD'
def electric_field_undulator(omega1=2.53465927101*10**17, trajectory=np.zeros((11, 10)), x=0.0, y=0.0, D=30.0):
    N = trajectory.shape[0]
    # in meters :
    R = np.sqrt(x ** 2 + y ** 2 + D ** 2)
    n_chap = np.array([x, y, D]) / R
    #in radian :
    #n_chap = np.array([x,y,1.0-0.5*(x**2 + y**2)])
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
    gamma = E / mc2
    omega1 = ((2.0 * gamma ** 2) / (1.0 + (K ** 2) / 2.0)) * ((2.0 * np.pi * c) / lambda_u)
    Elec= electric_field_undulator( omega1=omega1 ,trajectory=trajectory,x=x, y=y , D=D )
    norm_2 = (np.abs(Elec[0]) ** 2 + np.abs(Elec[1])** 2 + np.abs(Elec[2])** 2)
    return ((e**2 * omega1**2)/(16.0*np.pi**3*c*epsilon0))*norm_2


def Radiation(K=1.87, E=1.3 *10**9 , lambda_u=0.035, trajectory=np.zeros((11, 10)),D=30.0 ,X=np.arange(-0.0011, 0.0011, 0.00002), Y=np.arange(-0.0011, 0.0011, 0.00002)):
    res = np.zeros((len(X), len(Y)))
    print(' dimention de la matrice dna la fct MAG')
    print(res.shape)
    for i in range(len(X)):
        print(i)
        for k in range(len(Y)):
            res[i][k] = dI_dOmega(K=K, E=E, lambda_u=lambda_u, trajectory=trajectory,x=X[i] , y=Y[k],D=D)
    return res



K = 1.87
E = 1.3 * 10 ** 9
lambda_u = 0.035
Nb_period = 10
Nb_pts = 20


gamma = E / mc2
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

omega1 = ((2.0 * gamma ** 2) / (1.0 + (K ** 2) / 2.0)) * ((2.0 * np.pi * c) / lambda_u)
print('1st harmonic')
print(omega1)

##
omega_u = Beta_et * c * 2.0 * np.pi / lambda_u
t = np.linspace(-(lambda_u / (c * Beta_et)) * (Nb_period / 2), (lambda_u / (c * Beta_et)) * (Nb_period / 2), N)
By=Bo*np.sin(omega_u*t)
Vo=( e*Beta_et*Bo/(gamma*me*omega_u))*np.cos((-1.0)*Nb_period*np.pi)
Zo = (-lambda_u*Nb_period)/(2.0*c)
Xo=0.0



T = trajectory_undulator_reference(K=K,gamma=gamma, lambda_u=lambda_u, Nb_period=Nb_period, Nb_point=Nb_pts, Beta_et=Beta_et)
T2 = trajectory_undulator(By=By,t=t,gamma=gamma, Beta=Beta, Beta_et=Beta_et,Vo=Vo , Xo=Xo,Zo=Zo)
print("number of points : %d" % (T.shape[0]))
print("number of points : %d" % (T2.shape[0]))

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


# # # 3D plot
fig = figure()
ax = Axes3D(fig)
X=np.arange(0.0, 0.0301, 0.0003)
Y=np.arange(0.0, 0.0301, 0.0003)
Z2 = Radiation(K=K,E=E,trajectory=T,X=X,Y=Y)
X,Y = np.meshgrid(X,Y)
print('debut du plot')
ax.plot_surface(X,Y, Z2, rstride=1, cstride=1, cmap='hot')
ax.set_xlabel("X")
ax.set_ylabel('Y')
ax.set_zlabel("flux")
show()