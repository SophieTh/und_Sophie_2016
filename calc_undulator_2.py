import numpy as np
import matplotlib.pyplot as plt
from pylab import *
from mpl_toolkits.mplot3d import Axes3D
import scipy.constants as codata
import scipy.integrate as integrate
from scipy.integrate import ode
from scipy.integrate import odeint
from scipy.interpolate import interp1d



#calculate a theorical trajectory in an undulator
def analytical_trajectory_undulator(K=1.87 ,E=1.3e9,lambda_u=0.020, Nb_period=10, Nb_point=10):
    N = Nb_period * Nb_point + 1
    ku= 2.0*np.pi/lambda_u
    gamma = E / 0.511e6
    Beta = np.sqrt(1.0 - (1.0 / gamma ** 2))
    Beta_et = Beta * (1.0 - (K / (2.0 * gamma)) ** 2)
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
    trajectory[1] = (K/(gamma*Beta_et*ku*codata.c))* np.sin(omega_u * trajectory[0])
    # Vx et Vz en fct de t
    trajectory[6] = Beta_et - ((K/gamma)**2) *(1.0/4.0)*np.cos(2.0*omega_u * trajectory[0])
    trajectory[4] = (K/(gamma*ku*Beta_et*codata.c))*omega_u* np.cos(omega_u * trajectory[0])
    # Ax et Az en fct de t
    trajectory[9] = ((omega_u*K/gamma)**2) * (1.0/(2.0*ku*codata.c))*np.sin(2.0*omega_u * trajectory[0])
    trajectory[7] = -(K/(gamma*ku*Beta_et*codata.c))*(omega_u**2)* np.sin(omega_u * trajectory[0])
    # Bo = K / (93.4 * lambda_u)
    # By = -Bo * np.sin((2.0 * np.pi / lambda_u) *codata.c*trajectory[0])
    # plt.plot(codata.c*trajectory[0], By)
    # plt.show()
    return trajectory

#simulate the magnatic field in a undulator :
def creation_magnetic_field(K,lambda_u,Nb_period,z) :
    Bo = K / (93.4 * lambda_u)
    By = -Bo * np.sin((2.0 * np.pi / lambda_u) * z)
    # Hamming windowing
    windpar=1.0/(2.0*Nb_period)
    zmin = z.min()
    apo1 = zmin + windpar
    apo2 = z.max() - windpar
    wind = np.ones(len(z))
    for i in range(len(z)):
        if z[i] <= apo1:
            wind[i] *= 1.08 - (.54+0.46*np.cos(np.pi*(z[i]-zmin)/windpar))
        if z[i] >= apo2:
            wind[i] *=   1.08 - (.54-0.46*np.cos(np.pi*(z[i]-apo2)/windpar))
    By *= wind
    return By


def fct_ODE(y,t,cst,B) :
    return [ -cst*B( y[5]) * y[2],
            0.0,
             cst * B(y[5]) * y[0],
            y[0],
            0.0,
            y[2]]

# enelarge 2 vector for the interpolation of the magnetic field use in trajectory_undulator_from_magnetic_field2
def enlargement_vector_for_interpolation(Z,By,nb_enlarg=0) :
    dz = Z[1] - Z[0]
    enlarg_1 = np.linspace(Z[0] - nb_enlarg * dz, Z[0] - dz, nb_enlarg)
    enlarg_2 = np.linspace(Z[len(Z) - 1] + dz, Z[len(Z) - 1] + nb_enlarg * dz, nb_enlarg)
    Z = np.concatenate((enlarg_1, Z))
    Z = np.concatenate((Z, enlarg_2))
    enlarg_3 = np.zeros(nb_enlarg)
    By = np.concatenate((enlarg_3, By))
    By = np.concatenate((By, enlarg_3))
    return Z,By

# electron's trajectory in a PLANE undulator that mean :  B=(0,By,0)
# other hypothesis norm(v)=constant
def trajectory_undulator_from_magnetic_field1( By=np.zeros(101),Z= np.zeros(101),K=1.87,E=1.3e9,N=101,
                                               Vx=0.0,Xo=0.0,Zo=-0.175):
    gamma = E / 0.511e6
    Beta = np.sqrt(1.0 - (1.0 / gamma ** 2))
    Beta_et = Beta * (1.0 - (K / (2.0 * gamma)) ** 2)
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
    trajectory[0] = np.linspace(Z[0] / (codata.c * Beta_et),Z[len(Z)-1]/ (codata.c * Beta_et),N)
    #trajectory[0] = Z/(codata.c ) ???

    # Ax(t)
    Xm= codata.e*Beta_et/(gamma*codata.m_e)
    trajectory[7] = Xm * By
    # Vx et Vz
    for i in range(N):
        trajectory[4][i] = np.trapz(trajectory[7][0:(i + 1)], trajectory[0][0:(i + 1)]) + Vx/codata.c
    trajectory[6] = np.sqrt((Beta)**2 - trajectory[4]**2)
    # X et Z
    for i in range(N):
        trajectory[1][i] = np.trapz(trajectory[4][0:(i + 1)], trajectory[0][0:(i + 1)]) + Xo/codata.c
        trajectory[3][i] = np.trapz(trajectory[6][0:(i + 1)], trajectory[0][0:(i + 1)]) + Zo/codata.c
     #Az
    trajectory[9]=-(trajectory[7]*trajectory[4])/trajectory[6]

    return trajectory

#electron's trajectory in a PLANE undulator that mean :  B=(0,By,0)
def trajectory_undulator_from_magnetic_field2(By,Z,E,N,nb_enlarg,Vx=0.0,Vz=0.999997*codata.c,Xo=0.0,Zo=-0.175) :
    gamma = E / 0.511e6
    Beta = np.sqrt(1.0 - (1.0 / gamma ** 2))
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
    trajectory[0] = np.linspace(Z[0]/codata.c,Z[len(Z)-1]/codata.c, N)
    Vo = [Vx, 0.0, Vz, Xo, 0.0, Z[0]]
    Z,By=enlargement_vector_for_interpolation(Z, By, nb_enlarg)
    B = interp1d(Z, By)
    print(" min et max de B")
    print(B(Z).min())
    print(B(Z).max())
    cst=-codata.e/(codata.m_e*gamma)
    print("cst ODE")
    print(cst)
    res = odeint(fct_ODE,Vo,trajectory[0],args=(cst,B),full_output=True)
    traj=res[0]
    info=res[1]
    print("1 : nonstiff problems, Adams . 2: stiff problem, BDF")
    print(info.get('mused'))
    traj = np.transpose(traj)
    trajectory[4] = traj[0]
    trajectory[5] = traj[1]
    trajectory[6] = traj[2]
    trajectory[1] = traj[3]
    trajectory[2] = traj[4]
    trajectory[3] = traj[5]
    trajectory[7] = -cst * B(trajectory[3]) * trajectory[6]
    trajectory[9] = cst  * B(trajectory[3]) * trajectory[4]
    k=1
    while k<10 :
        trajectory[k] *= 1.0/codata.c
        k+=1
    return trajectory


r"""
    undulator_trajectory(K,E,lambda_u,Nb_period,Nb_point,Z_By=None,type_trajectory=1)

  PURPOSE:
 	This procedure calculates the electron trajectory in a undulator

  INPUTS:
    Undulator's parameters :
        K
        E : energy (eV)
        lambda_u : period of the undulator (m)
    Trajectory's parameters :
        Nb_points : number of points for one period
        Z_By : numpy.array(2,npoints)   (CASE A)
                Z_By[0] = the coordinate along the undulator's length
                Z_By[1] = a numpy.array which is the magnatic field mesured in fonction of Z
               None                     (CASE B)
                The fonction create a theoritical magnatic fields
    type_trajectory : int
        1 : calculate by integration
        2 : calculte with an ODE
        other : a theorytical trajectory is return

  OUTPUT:
 	trajectory = numpy.array(10,npoints) the trajectory in fonction of the observator's time and divide by c ,
 	                (the light speed).
 	    trajectory[0]= times ( NOT DIVIDE BY C )
 	    trajectory[1] = X/c
 	    trajectory[2] = Y/c
 	    trajectory[3] = Z/c
 	    trajectory[4] = Vx/c
 	    trajectory[5] = Vy/c
 	    trajectory[6] = Vz/c
 	    trajectory[7] = Ax/c
 	    trajectory[8] = Ay/c
 	    trajectory[9] = Az/c

    REMARQUE :
        CASE A :
            - Z_By[0] must be centred around 0.0
            - With the case A and the trajectory_type 1 the number of points put in argument is not considerate
                the number npoints is the length of the magnetic field given


;-
"""
def undulator_trajectory(K,E,lambda_u,Nb_period=10,Nb_point=100,Z_By=None,type_trajectory=0,
                         Vx=0.0,Vz=0.999997*codata.c,Xo=0.0) :

    if (type_trajectory == 1 or type_trajectory == 2) :

        # CASE B
        if Z_By == None :
            Z = np.linspace(-lambda_u * (Nb_period / 2),lambda_u * (Nb_period / 2), Nb_period * Nb_point + 1)
            By= creation_magnetic_field(K,lambda_u,Nb_period,Z)

        # CASE A
        else :
            Z= Z_By[0]
            By=Z_By[1]
        # plt.plot(Z, By)
        # plt.show()

        if (type_trajectory == 1) :

            trajectory = trajectory_undulator_from_magnetic_field1(By=By, Z=Z,K=K,E=E,N=len(By),Vx=Vx,Xo=Xo,Zo=Z[0])

        else :
            trajectory = trajectory_undulator_from_magnetic_field2(By=By, Z=Z,E=E,N=Nb_period * Nb_point + 1
                                                                   ,nb_enlarg=np.floor(len(Z)/(0.1*Nb_point)),
                                                                   Vx=Vx,Vz=Vz,Xo=Xo,Zo=Z[0])

    else :
        trajectory = analytical_trajectory_undulator(K=K, lambda_u=lambda_u, Nb_period=Nb_period, Nb_point=Nb_point)

    return trajectory


#draw all coordinate of the trajectory in function of the time
def draw_trajectory( trajectory = np.zeros((10,101))) :
    plt.plot(trajectory[0],trajectory[1])
    plt.title(" X = f(t) ")
    plt.xlabel('t')
    plt.ylabel('x')
    plt.show()

    plt.plot(trajectory[0],trajectory[3])
    plt.title(" Z  = f(t) ")
    plt.xlabel('t')
    plt.ylabel('Z')
    plt.show()

    print('average of Vz  =')
    ## vrai valeur de Beta_et
    Beta_et= np.sum(trajectory[6])/len(trajectory[6])
    print(Beta_et)

    Z = Beta_et*trajectory[0]
    plt.plot(trajectory[0],trajectory[3]-Z)
    plt.title(" Z - Beta* t = f(t) ")
    plt.xlabel('t')
    plt.ylabel('Z - Beta* t')
    plt.show()

    plt.plot(trajectory[0], trajectory[4])
    plt.title(" Vx = f(t) ")
    plt.xlabel('t')
    plt.ylabel('Vx')
    plt.show()

    plt.plot(trajectory[0], trajectory[6])
    plt.title(" Vz = f(t) ")
    plt.xlabel('t')
    plt.ylabel('Vz')
    plt.show()

    plt.plot(trajectory[0], trajectory[6]-Beta_et)
    plt.title(" Vz = f(t) ")
    plt.xlabel('t')
    plt.ylabel('Vz')
    plt.show()

    plt.plot(trajectory[0], trajectory[7])
    plt.title(" Ax = f(t) ")
    plt.xlabel('t')
    plt.ylabel('Ax')
    plt.show()

    plt.plot(trajectory[0], trajectory[9])
    plt.title(" Az = f(t) ")
    plt.xlabel('t')
    plt.ylabel('Az')
    plt.show()

#draw all coordinate of the 2 trajectories in function of the time
# It must have the same shape and the same fistr vector wich represent time.
def draw_2_trajectory( trajectory1 = np.zeros((10,101)),trajectory2 = np.zeros((10,101))) :

    plt.plot(trajectory1[0],trajectory2[0])
    plt.title(" times 1 = f(times2) ")
    plt.xlabel('times1')
    plt.ylabel('times2')
    plt.show()


    plt.plot(trajectory1[0],trajectory1[1])
    plt.plot(trajectory2[0],trajectory2[1])
    plt.title(" X = f(t) ")
    plt.xlabel('t')
    plt.ylabel('x')
    plt.show()

    plt.plot(trajectory1[0],trajectory1[3])
    plt.plot(trajectory2[0],trajectory2[3])
    plt.title(" Z  = f(t) ")
    plt.xlabel('t')
    plt.ylabel('Z')
    plt.show()

    plt.plot(trajectory1[0],trajectory1[4])
    plt.plot(trajectory2[0],trajectory2[4])
    plt.title(" Vx = f(t) ")
    plt.xlabel('t')
    plt.ylabel('Vx')
    plt.show()

    plt.plot(trajectory1[0],trajectory1[6])
    plt.plot(trajectory2[0],trajectory2[6])
    plt.title(" Vz = f(t) ")
    plt.xlabel('t')
    plt.ylabel('Vz')
    plt.show()

    plt.plot(trajectory1[0],trajectory1[7])
    plt.plot(trajectory2[0],trajectory2[7])
    plt.title(" Ax = f(t) ")
    plt.xlabel('t')
    plt.ylabel('Ax')
    plt.show()

    plt.plot(trajectory1[0],trajectory1[9])
    plt.plot(trajectory2[0],trajectory2[9])
    plt.title(" Az = f(t) ")
    plt.xlabel('t')
    plt.ylabel('Az')
    plt.show()


#exact equation for the energy radiated
#warning !!!!!!! far field approximation
def energy_radiated(omega=2.53465927101*10**17,gamma=2000.0,trajectory=np.zeros((11,10)) , x=0.00 , y=0.0, D=None):
    N = trajectory.shape[1]

    # far field !
    if D == None:
        # in radian :
        n_chap = np.array([x, y, 1.0 - 0.5 * (x ** 2 + y ** 2)])
        R=100.0
    else:
        # in meters :
        n_chap =np.array([x,y,D])
        R=norm(n_chap)
        n_chap /= R

    E = np.full((3,), 0. + 1j * 0., dtype=np.complex)
    integrand = np.full((3,N), 0. + 1j * 0., dtype=np.complex)
    Alpha=trajectory[7]*(n_chap[2]-trajectory[6]) - (n_chap[0]-trajectory[4])*trajectory[9]
    Alpha2=np.exp(0. + 1j * omega * (trajectory[0] - n_chap[0]*trajectory[1]-n_chap[2]*trajectory[3]))
    Alpha1=(1.0 / (1.0 - n_chap[0]*trajectory[4]-n_chap[2]*trajectory[6])) ** 2
    Alpha3= codata.c/(gamma**2 * R)
    integrand[0] += (-(n_chap[1]**2)*trajectory[7]-n_chap[2]*Alpha
                     + Alpha3 *(n_chap[0]-trajectory[4]))*Alpha2*Alpha1
    integrand[1] += (n_chap[1]*(n_chap[0]*trajectory[7]+n_chap[2]*trajectory[9])
                            + Alpha3 * (n_chap[1] - trajectory[5]))*Alpha2*Alpha1
    integrand[2] += (-(n_chap[1]**2)*trajectory[9]+n_chap[0]*Alpha
                     + Alpha3 *(n_chap[2]-trajectory[6]))*Alpha2*Alpha1
    for k in range(3):
        E[k] = np.trapz(integrand[k], trajectory[0])
    return (np.abs(E[0]) ** 2 + np.abs(E[1])** 2 + np.abs(E[2])** 2)

# energy radiated without the the far filed approxiamation :
#doesn't work for the moment....
def energy_radiated_no_far_field(omega=2.53465927101*10**17,gamma=2000.0,trajectory=np.zeros((11,10)) , x=0.00 , y=0.0, D=100.0):
    N = trajectory.shape[1]
    n_chap =np.array([x/codata.c-trajectory[1],y/codata.c-trajectory[2],D/codata.c-trajectory[3]])
    R=np.zeros(n_chap.shape[1])
    for i in range(n_chap.shape[1]) :
        R[i] = norm(n_chap[:,i])
        n_chap[:,i] /= R[i]

    E = np.full((3,), 0. + 1j * 0., dtype=np.complex)
    integrand = np.full((3,N), 0. + 1j * 0., dtype=np.complex)
    Alpha=trajectory[7]*(n_chap[2]-trajectory[6]) - (n_chap[0]-trajectory[4])*trajectory[9]
    Alpha2=np.exp(0. + 1j * omega * (trajectory[0] + R))
    #Alpha2 = np.exp(0. + 1j * omega * (trajectory[0] - n_chap[0] * trajectory[1] - n_chap[2] * trajectory[3]))
    Alpha1=(1.0 / (1.0 - n_chap[0]*trajectory[4]-n_chap[2]*trajectory[6])) ** 2
    Alpha3= codata.c/(gamma**2 * R)
    integrand[0] += (-(n_chap[1]**2)*trajectory[7]-n_chap[2]*Alpha
                     + Alpha3 *(n_chap[0]-trajectory[4]))*Alpha2*Alpha1
    integrand[1] += (n_chap[1]*(n_chap[0]*trajectory[7]+n_chap[2]*trajectory[9])
                            + Alpha3 * (n_chap[1] - trajectory[5]))*Alpha2*Alpha1
    integrand[2] += (-(n_chap[1]**2)*trajectory[9]+n_chap[0]*Alpha
                     + Alpha3 *(n_chap[2]-trajectory[6]))*Alpha2*Alpha1
    for k in range(3):
        E[k] = np.trapz(integrand[k], trajectory[0])
    return (np.abs(E[0]) ** 2 + np.abs(E[1])** 2 + np.abs(E[2])** 2)

# approximation of the energy radiated by an electron in a PLANE undulator
# warning : trajectory is the trajectory difine like the function "undulator trajectory" before :
def energy_radiated_initial(omega=2.53465927101*10**17,trajectory=np.zeros((11,10)) , x=0.00 , y=0.0, D=None):
    N = trajectory.shape[1]

    if D == None:
        # in radian :
        n_chap = np.array([x, y, 1.0 - 0.5 * (x ** 2 + y ** 2)])
    else:
        # in meters :
        R = np.sqrt(x ** 2 + y ** 2 + D ** 2)
        n_chap = np.array([x, y, D]) / R

    E = np.full((3,), 0. + 1j * 0., dtype=np.complex)
    integrand = np.full((3,N), 0. + 1j * 0., dtype=np.complex)
    Alpha=trajectory[7]*(n_chap[2]-trajectory[6]) - (n_chap[0]-trajectory[4])*trajectory[9]
    Alpha2=np.exp(0. + 1j * omega * (trajectory[0] - n_chap[0]*trajectory[1]-n_chap[2]*trajectory[3]))
    Alpha1=(1.0 / (1.0 - n_chap[0]*trajectory[4]-n_chap[2]*trajectory[6])) ** 2
    integrand[0] += (-(n_chap[1]**2)*trajectory[7]-n_chap[2]*Alpha)*Alpha2*Alpha1
    integrand[1] += n_chap[1]*(n_chap[0]*trajectory[7]+n_chap[2]*trajectory[9])*Alpha2*Alpha1
    integrand[2] += (-(n_chap[1]**2)*trajectory[9]+n_chap[0]*Alpha)*Alpha2*Alpha1
    for k in range(3):
        E[k] = np.trapz(integrand[k], trajectory[0])
    return (np.abs(E[0]) ** 2 + np.abs(E[1])** 2 + np.abs(E[2])** 2)

# Photon's flow all over a screen situate at distance D of an undulator
# warning : trajectory is the trajectory difine like the function "undulator trajectory" before :
def radiation_single_electron(K=1.87, E=1.3 * 10 ** 9, lambda_u=0.035, trajectory=np.zeros((11, 10)), D=None, omega=None,
              X=np.arange(-0.0011, 0.0011, 0.00002), Y=np.arange(-0.0011, 0.0011, 0.00002)):

    gamma = E / 0.511e6

    if omega == None:  # use the first harmonic at the resonance
        omega = ((2.0 * gamma ** 2) / (1.0 + (K ** 2) / 2.0)) * ((2.0 * np.pi * codata.c) / lambda_u)


    # c1 = codata.e ** 2 * omega1 ** 2 / (16 * np.pi ** 3 * codata.epsilon_0 * codata.c )
    # c2 = 1.0 / codata.e  # multiply by number of electrons in 1 A
    # c3 = 2.0*np.pi / (codata.h * omega1)  # divide by e energy (to get number of photons)
    # c4 = 1e-2 / omega1  # to get 1% energy bandwidth  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # c5 = 0.1e-3 * 0.1e-3  # from rad to .1mrad angle bandwidth
    c6= codata.e*1e-10/(8.0*np.pi**2*codata.epsilon_0*codata.c*codata.h)
    #c6 = codata.e**2 / (16.0* np.pi ** 3 * codata.epsilon_0 * codata.c)

    # res = np.zeros((len(X), len(Y)))
    # for i in range(len(X)):
    #     #print(i)
    #     for k in range(len(Y)):
    #         res[i][k] = c6*energy_radiated(omega=omega,trajectory=trajectory , x=X[i] , y=Y[k], D=D )

    if X.size != Y.size:
        raise Exception("X and Y dimensions must be equal.")

    res = np.zeros_like(X)
    shape1 = res.shape

    X = X.flatten()
    Y = Y.flatten()
    res = res.flatten()
    print(len(X))

    for i in range(len(X)):
        res[i] = c6*energy_radiated_initial(omega=omega,trajectory=trajectory , x=X[i] , y=Y[i], D=D )
        #res[i]= c6*energy_radiated_no_far_field(omega=omega,trajectory=trajectory, x=X[i], y=Y[i], D=D)
        #res[i] = c6 * energy_radiated(omega=omega,gamma=gamma, trajectory=trajectory, x=X[i], y=Y[i], D=D)

    X = X.reshape(shape1)
    Y = Y.reshape(shape1)
    res = res.reshape(shape1)
    print("radiation max ")
    print(res.max())
    return res


# approximation of the energy radiated by an electron in a PLANE undulator
# Use the second formula write in JACKSON's book .... there is a factor's difference : 5.203e21 with the other formula
# warning : trajectory is the trajectory difine like the function "undulator trajectory" before :
def energy_radiated2(omega=2.53465927101*10**17,gamma=2000.0,trajectory=np.zeros((11,10)) , x=0.00 , y=0.0, D=None):
    N = trajectory.shape[1]

    if D == None:
        # in radian :
        n_chap = np.array([x, y, 1.0 - 0.5 * (x ** 2 + y ** 2)])
    else:
        # in meters :
        R = np.sqrt(x ** 2 + y ** 2 + D ** 2)
        n_chap = np.array([x, y, D]) / R

    E = np.full((3,), 0. + 1j * 0., dtype=np.complex)
    integrand = np.full((3,N), 0. + 1j * 0., dtype=np.complex)
    Alpha=trajectory[7]*n_chap[2] - n_chap[0]*trajectory[9]
    Alpha2 = np.exp(0. + 1j * omega * (trajectory[0]- n_chap[0] * trajectory[1] - n_chap[2] * trajectory[3]))
    integrand[0] += (-(n_chap[1]**2)*trajectory[7]-n_chap[2]*Alpha)*Alpha2
    integrand[1] += n_chap[1]*(n_chap[0]*trajectory[7]+n_chap[2]*trajectory[9])*Alpha2
    integrand[2] += (-(n_chap[1]**2)*trajectory[9]+n_chap[0]*Alpha)*Alpha2
    for k in range(3):
        E[k] = np.trapz(integrand[k], trajectory[0])
    E *= 1j*(-omega)

    Alpha3= -trajectory[0][0] + n_chap[2]*trajectory[3][0]
    Beta=trajectory[6][0]
    terme_bord = np.full((3,), 0. + 1j * 0., dtype=np.complex)
    terme_bord[0] =n_chap[2]*n_chap[0]*Beta +0.0*1j
    terme_bord[1]=n_chap[2]*n_chap[1]*Beta+0.0*1j
    terme_bord[2]=Beta*(-n_chap[0]**2 - n_chap[1]**2)+0.0*1j
    terme_bord *= 2j*np.sin(omega*Alpha3)/(1.0-n_chap[2]*trajectory[6][0])
    E += terme_bord

    return (np.abs(E[0]) ** 2 + np.abs(E[1])** 2 + np.abs(E[2])** 2)


# Same that before bur call the energy radiated 2
def radiation_single_electron2(K=1.87, E=1.3 * 10 ** 9, lambda_u=0.035, trajectory=np.zeros((11, 10)), D=None, omega=None,
              X=np.arange(-0.0011, 0.0011, 0.00002), Y=np.arange(-0.0011, 0.0011, 0.00002)):

    gamma = E / 0.511e6

    if omega == None:  # use the first harmonic at the resonance
        omega = ((2.0 * gamma ** 2) / (1.0 + (K ** 2) / 2.0)) * ((2.0 * np.pi * codata.c) / lambda_u)


    if X.size != Y.size:
        raise Exception("X and Y dimensions must be equal.")

    res = np.zeros_like(X)
    shape1 = res.shape

    X = X.flatten()
    Y = Y.flatten()
    res = res.flatten()

    shape2 = res.shape

    c6 = codata.e * 1e-10 / (8.0 * np.pi ** 2 * codata.epsilon_0 * codata.c * codata.h)
    #c6 = codata.e ** 2 / (16.0 * np.pi ** 3 * codata.epsilon_0 * codata.c)
    for i in range(len(X)):
        res[i] = c6*energy_radiated2(omega=omega,gamma=gamma,trajectory=trajectory , x=X[i] , y=Y[i], D=D )

    X = X.reshape(shape1)
    Y = Y.reshape(shape1)
    res = res.reshape(shape1)
    print("radiation max ")
    print(res.max())
    return res


