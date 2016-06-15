import numpy as np
import matplotlib.pyplot as plt
from pylab import *
from mpl_toolkits.mplot3d import Axes3D
import scipy.constants as codata
import scipy.integrate as integrate
from scipy.interpolate import interp1d


from calc_undulator_2 import draw_trajectory,undulator_trajectory, draw_2_trajectory,radiation_single_electron_2D


def trajectory_reference_test(K = 1.87,E = 1.3e9,lambda_u = 0.035,Nb_period = 12,N=1001,Z=np.zeros(1001)) :
    gamma = E / 0.511e6
    Beta = np.sqrt(1.0 - (1.0 / gamma ** 2))
    Beta_et = 1.0 - (1.0 / (2.0 * gamma ** 2)) * (1.0 + (K ** 2) / 2.0)
    trajectory = np.zeros((10, N))
    trajectory[0] = np.linspace(Z[0] / (Beta_et * codata.c), Z[len(Z) - 1] / (Beta_et * codata.c), N)
    Bo = K / (93.4 * lambda_u)
    cst= -codata.e/(gamma*codata.m_e)
    trajectory[1] = (Beta / (cst * Bo)) * np.sin(cst * Bo * trajectory[0])
    trajectory[3] = (-Beta / (cst * Bo*codata.c)) * np.cos(cst * Bo * trajectory[0])
    trajectory[4] = Beta*np.cos(cst * Bo * trajectory[0])
    trajectory[6] = Beta*np.sin(cst * Bo * trajectory[0])
    trajectory[7] = -cst * Bo*Beta * np.sin(cst * Bo * trajectory[0])
    trajectory[9] = cst * Bo*Beta * np.cos(cst * Bo * trajectory[0])
    return trajectory

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

    #print("erreur temps")
    egalite_temp = (ref[0] - test[0])
    #print(egalite_temp.max()/test[0].max())
    erreur01 = np.zeros_like(ref)
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
    # print("proportion erreur.max / ref.max : ")
    # print(' erreur x max=')
    # print((abs(erreur01[1])).max() / (abs(ref[1])).max())
    # print(' erreur z max=')
    # print((abs(erreur01[3])).max() / (abs(ref[3])).max())
    # print('erreur Vx max=')
    # print((abs(erreur01[4])).max() / (abs(ref[4])).max())
    # print(' erreur Vz max=')
    # print((abs(erreur01[4])).max() / (abs(ref[4])).max())
    # print('erreur Ax max=')
    # print((abs(erreur01[6])).max() / (abs(ref[6])).max())
    # print('erreur Az max=')
    # print((abs(erreur01[7])).max() / (abs(ref[7])).max())

    return erreur01

def erreur_of_radiation_for_trajectory(X=np.arange(0.0, 0.0301, 0.0003),Y=np.zeros(101),
                                       ref=np.zeros((10,101)),test=np.zeros((10,101))) :
    Z0, maxref = radiation_single_electron_2D(K=K, E=E, trajectory=ref, X=X, Y=Y, D=30.0)
    Z1, maxtest = radiation_single_electron_2D(K=K, E=E, trajectory=test, X=X, Y=Y, D=30.0)

    u=np.linspace(0.0,len(X),len(X))/len(X)
    plt.plot(u, Z0)
    plt.plot(u, Z1)
    plt.title(" ref and test = f(u) ")
    plt.show()

    plt.plot(u, (abs(Z1 - Z0) / maxref))
    plt.title(" relativ error of Z1 for Z0 ")
    plt.ylabel('|Z1-Z0|/max0')
    plt.xlabel('u')
    plt.show()

    print("erreur radiation max ")
    print((abs(Z1 - Z0)).max()/maxref)

    return abs(Z1-Z0)


def compare_with_reference_traj(K,E,lambda_u,Nb_period,Nb_pts,
                                X=np.arange(0.0, 0.0301, 0.0003), Y=np.zeros(101)) :
    # a magnetic field

    gamma = E /0.511e6
    Beta = np.sqrt(1.0 - (1.0 / gamma ** 2))
    Beta_et = 1.0-(1.0/(2.0*gamma**2))*(1.0+(K**2)/2.0)
    #Beta_et=Beta*(1.0-(K/(2.0*gamma))**2)
    Bo = K / (93.4 * lambda_u)
    N = Nb_period * Nb_pts + 1


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

    print(" --------------------------------- ")
    print('**')
    print("T0 and T1")
    print('**')
    print(" --------------------------------- ")

    #draw_2_trajectory(T0,T1)
    erreur01 =erreur_trajectory(ref=T0,test=T1)
    # draw_trajectory(erreur01)

    erreur_rad_01 = erreur_of_radiation_for_trajectory(X=X, Y=Y, ref=T0, test=T1)

    print(" --------------------------------- ")
    print('**')
    print("T0 and T2")
    print('**')
    print(" --------------------------------- ")

    #draw_2_trajectory(T0,T2)
    erreur02 =erreur_trajectory(ref=T0,test=T2)
    #draw_trajectory(erreur02)

    erreur_rad_02 = erreur_of_radiation_for_trajectory(X=X, Y=Y, ref=T0, test=T2)

    print(" --------------------------------- ")
    print('**')
    print("T1 and T2")
    print('**')
    print(" --------------------------------- ")

    #draw_2_trajectory(T0,T2)
    erreur12 =erreur_trajectory(ref=T1,test=T2)
    #draw_trajectory(erreur02)

    erreur_rad_12=erreur_of_radiation_for_trajectory(X=X,Y=Y,ref=T1,test=T2)

    return erreur01,erreur_rad_01,erreur02,erreur_rad_02,erreur12,erreur_rad_12

def compare_with_reference_traj_test(K, E, lambda_u, Nb_period, Nb_pts):
    gamma = E / 0.511e6
    Beta = np.sqrt(1.0 - (1.0 / gamma ** 2))
    Beta_et = 1.0 - (1.0 / (2.0 * gamma ** 2)) * (1.0 + (K ** 2) / 2.0)
    # Beta_et=Beta*(1.0-(K/(2.0*gamma))**2)
    Bo = K / (93.4 * lambda_u)
    N = Nb_period * Nb_pts + 1
    omega1 = ((2.0 * gamma ** 2) / (1.0 + (K ** 2) / 2.0)) * ((2.0 * np.pi * codata.c) / lambda_u)

    Z = np.linspace(-(2.0*np.pi*gamma*Beta_et*codata.c*codata.m_e/(codata.e*Bo)) * (Nb_period / 2),
                    (2.0*np.pi*gamma*Beta_et*codata.c*codata.m_e/(codata.e*Bo) )* (Nb_period / 2), N)
    By = np.ones_like(Z)*Bo

    Z_By = np.zeros((2, len(Z)))
    Z_By[0] = Z
    Z_By[1] = By

    T0 = trajectory_reference_test(K =K,E =E,lambda_u = lambda_u,Nb_period = Nb_period,N=N,Z=Z)

    T2 = undulator_trajectory(K, E, lambda_u, Nb_period, Nb_pts, Z_By, type_trajectory=2,
                          Vx=T0[4][0] * codata.c, Vz=T0[6][0] * codata.c, Xo=T0[1][0] * codata.c)

    # print(" --------------------------------- ")
    # print('**')
    # print("T0 and T2")
    # print('**')
    # print(" --------------------------------- ")

    # draw_2_trajectory(T0,T2)
    erreur02 = erreur_trajectory(ref=T0, test=T2)
    # draw_trajectory(erreur02)

    draw_trajectory(T2)

    return erreur02




def compare_with_ideal_magnetism(K,E,lambda_u,Nb_period,Nb_pts,
                                 X=np.arange(0.0, 0.0301, 0.0003), Y=np.zeros(101)) :
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

    T1=undulator_trajectory(K,E,lambda_u,Nb_period,Nb_pts,type_trajectory=1,
                        Vx=0.0,Vz=Beta*codata.c,Xo=0.0)
    T2=undulator_trajectory(K,E,lambda_u,Nb_period,Nb_pts,type_trajectory=2,
                            Vx=0.0, Vz=Beta * codata.c, Xo=0.0)

    #draw_2_trajectory(T1,T2)
    erreur12 =erreur_trajectory(ref=T1,test=T2)
    #draw_trajectory(erreur12)

    erreur_rad_12=erreur_of_radiation_for_trajectory(X=X,Y=Y,ref=T1,test=T2)

    return erreur12 , erreur_rad_12




def compare_with_magnetism(K, E, lambda_u, Nb_period, Nb_pts,Z_By,
                                     X=np.arange(0.0, 0.0301, 0.0003), Y=np.zeros(101)):

    gamma = E / 0.511e6
    print('gamma =')
    print(gamma)
    Beta = np.sqrt(1.0 - (1.0 / gamma ** 2))
    print('Beta = ')
    print(Beta)
    Beta_et = 1.0 - (1.0 / (2.0 * gamma ** 2)) * (1.0 + (K ** 2) / 2.0)
    # Beta_et=Beta*(1.0-(K/(2.0*gamma))**2)
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

    T1 = undulator_trajectory(K, E, lambda_u, Nb_period, Nb_pts,Z_By=Z_By, type_trajectory=1,
                      Vx=0.0, Vz=Beta * codata.c, Xo=0.0)
    T2 = undulator_trajectory(K, E, lambda_u, Nb_period, Nb_pts,Z_By=Z_By, type_trajectory=2,
                      Vx=0.0, Vz=Beta * codata.c, Xo=0.0)

    # draw_2_trajectory(T1,T2)
    erreur12 = erreur_trajectory(ref=T1, test=T2)
    # draw_trajectory(erreur12)

    erreur_rad_12 = erreur_of_radiation_for_trajectory(X=X, Y=Y, ref=T1, test=T2)

    return erreur12, erreur_rad_12




########## semi traj #######

def semi_trajectory_reference_test(K = 1.87,E = 1.3e9,lambda_u = 0.035,Nb_period = 12,N=1001) :
    gamma = E / 0.511e6
    Beta = np.sqrt(1.0 - (1.0 / gamma ** 2))
    Beta_et = 1.0 - (1.0 / (2.0 * gamma ** 2)) * (1.0 + (K ** 2) / 2.0)
    trajectory = np.zeros((5, N))
    trajectory[0]=np.linspace(-(lambda_u / (codata.c * Beta_et)) * (Nb_period / 2),
                                              (lambda_u / (codata.c * Beta_et)) * (Nb_period / 2), N)
    Bo = K / (93.4 * lambda_u)
    cst= -codata.e/(gamma*codata.m_e)
    omega=2.0*np.pi*Beta_et*codata.c/lambda_u
    alpha=(cst*Bo/omega)*(1.0-np.cos(omega*trajectory[0]))
    alpha_prim=(cst*Bo*np.sin(omega*trajectory[0]))
    trajectory[1]=-np.sin(alpha)*Beta
    trajectory[2] = np.cos(alpha)*Beta
    trajectory[3] = -Beta*np.cos(alpha)*alpha_prim
    trajectory[1] = -Beta*np.sin(alpha)*alpha_prim

    return trajectory


def semi_trajectory_undulator_from_magnetic_field1_trap(K=1.87, E=1.3e9,lambda_u=0.0035, N=101,Nb_period=12,
                                                  Vx=0.0):

    gamma = E / 0.511e6
    Beta = np.sqrt(1.0 - (1.0 / gamma ** 2))
    Beta_et = 1.0 - (1.0 / (2.0 * gamma ** 2)) * (1.0 + (K ** 2) / 2.0)
    Bo = K / (93.4 * lambda_u)
    omega = 2.0*np.pi*Beta_et*codata.c/lambda_u
    trajectory = np.zeros((5, N))
    trajectory[0]=trajectory[0] = np.linspace(-(lambda_u / (codata.c * Beta_et)) * (Nb_period / 2),
                                              (lambda_u / (codata.c * Beta_et)) * (Nb_period / 2), N)

    By2 = Bo * np.sin(omega * trajectory[0])

    # Ax(t)
    Xm = codata.e * Beta_et / (gamma * codata.m_e)
    trajectory[3] = Xm * By2
    # Vx et Vz
    for i in range(N):
        trajectory[1][i] = np.trapz(trajectory[3][0:(i + 1)], trajectory[0][0:(i + 1)]) + Vx / codata.c
    trajectory[2] = np.sqrt((Beta) ** 2 - trajectory[4] ** 2)
    # X et Z
    # Az
    trajectory[4] = -(trajectory[3] * trajectory[1]) / trajectory[2]

    return trajectory

def semi_trajectory_undulator_from_magnetic_field1_simps(K=1.87, E=1.3e9,lambda_u=0.0035, N=101,Nb_period=12,
                                                  Vx=0.0):
    gamma = E / 0.511e6
    Beta = np.sqrt(1.0 - (1.0 / gamma ** 2))
    Beta_et = 1.0 - (1.0 / (2.0 * gamma ** 2)) * (1.0 + (K ** 2) / 2.0)
    Bo = K / (93.4 * lambda_u)
    omega = 2.0*np.pi*Beta_et*codata.c/lambda_u
    trajectory = np.zeros((5, N))
    trajectory[0] = trajectory[0] = np.linspace(-(lambda_u / (codata.c * Beta_et)) * (Nb_period / 2),
                                                (lambda_u / (codata.c * Beta_et)) * (Nb_period / 2), N)

    By2 = Bo*np.sin(omega*trajectory[0])

    # Ax(t)
    Xm = codata.e * Beta_et / (gamma * codata.m_e)
    trajectory[3] = Xm * By2
    # Vx et Vz
    for i in range(N):
        trajectory[1][i] = integrate.simps(trajectory[3][0:(i + 1)], trajectory[0][0:(i + 1)]) + Vx / codata.c
    trajectory[2] = np.sqrt((Beta) ** 2 - trajectory[4] ** 2)
    # X et Z
    # Az
    trajectory[4] = -(trajectory[3] * trajectory[1]) / trajectory[2]

    return trajectory

def erreur_2_semi_traj(ref,test):
    erreur =  np.zeros((ref.shape[0]-1,ref.shape[1]))
    for k in range(ref.shape[0]-1) :
        erreur[k]= np.abs(ref[k+1]-test[k+1])
    return erreur

#########################################

def compare_method_calc_trajectory(K = 1.87,E = 1.3e9,lambda_u = 0.035,Nb_period = 12) :
    erreur_max01= np.zeros((10, 100))
    erreur_max02 = np.zeros((10, 100))
    erreur_max12 = np.zeros((10, 100))
    for i in range(100) :
        #print(i)
        Nb_pts=i+1
        erreur_max01[0][i]=Nb_period*Nb_pts+1
        erreur_max02[0][i] = Nb_period * Nb_pts + 1
        erreur_max12[0][i] = Nb_period * Nb_pts + 1
        res01, res02, res12=compare_with_reference_traj_test(K, E, lambda_u, Nb_period, Nb_pts)
        for k in range(9) :
            erreur_max01[k+1][i]=  res01[k+1].max()
            erreur_max02[k+1][i] = res02[k+1].max()
            erreur_max12[k+1][i] = res12[k+1].max()

    return erreur_max01 , erreur_max02, erreur_max12




######################################

####                            #####
    # ************************ #
####                            #####

#######################################
K = 1.87
E = 1.3e9
lambda_u = 0.035
Nb_period = 12
Nb_pts = 100


X=np.arange(0.0, 0.0301, 0.0003)
Y=np.zeros_like(X)

# compare_with_reference_traj(K,E,lambda_u,Nb_period,Nb_pts,X=X,Y=Y)
#
# compare_with_ideal_magnetism(K,E,lambda_u,Nb_period,Nb_pts,X=X,Y=Y)
#
# reference=np.load("x_ray_booklet_field.npz")
# #print(reference.keys())
# Z=reference['ct']
# Z -= (Z[len(Z)-1])/2.0
# By2=reference['B_y']
#
# Z_By=np.zeros((2,len(Z)))
# Z_By[0]=Z
# Z_By[1]=By2

##compare_with_magnetism(K, E, lambda_u, Nb_period, Nb_pts, Z_By=Z_By,X=X, Y=Y)


compare_with_reference_traj_test(K, E, lambda_u, Nb_period, Nb_pts)

# resultat02=compare_method_calc_trajectory(K = K,E = E,lambda_u = lambda_u,Nb_period = Nb_period)
# print(resultat01.shape)
#
# draw_trajectory(resultat02)



