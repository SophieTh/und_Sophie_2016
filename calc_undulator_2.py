import numpy as np
import matplotlib.pyplot as plt
from pylab import *
from mpl_toolkits.mplot3d import Axes3D
import scipy.constants as codata
import scipy.integrate as integrate
from scipy.integrate import ode
from scipy.integrate import odeint
from scipy.interpolate import interp1d


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
    trajectory[1] = (K/(gamma*ku*codata.c))* np.sin(omega_u * trajectory[0])
    # Vx et Vz en fct de t
    trajectory[6] = Beta_et - ((K/gamma)**2) * (2.0*omega_u/(8.0*ku*codata.c ))*np.cos(2.0*omega_u * trajectory[0])
    trajectory[4] = (K/(gamma*ku*codata.c))*omega_u* np.cos(omega_u * trajectory[0])
    # Ax et Az en fct de t
    trajectory[9] = ((2.0*omega_u*K/gamma)**2) * (1.0/(8.0*ku*codata.c))*np.sin(2.0*omega_u * trajectory[0])
    trajectory[7] = -(K/(gamma*ku*codata.c))*(omega_u**2)* np.sin(omega_u * trajectory[0])
    # Bo = K / (93.4 * lambda_u)
    # By = -Bo * np.sin((2.0 * np.pi / lambda_u) *codata.c*trajectory[0])
    # plt.plot(codata.c*trajectory[0], By)
    # plt.show()
    return trajectory



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


# hypothesis : B=(0,By,0) and norm(v)=constant
def trajectory_undulator_from_magnetic_field1( By=np.zeros(101),Z= np.zeros(101),K=1.87,E=1.3e9,N=101,Vo=0.0,Xo=0.0,Zo=-0.175):
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
        trajectory[4][i] = np.trapz(trajectory[7][0:(i + 1)], trajectory[0][0:(i + 1)]) + Vo/codata.c
    trajectory[6] = np.sqrt((Beta)**2 - trajectory[4]**2)
    # X et Z
    for i in range(N):
        trajectory[1][i] = np.trapz(trajectory[4][0:(i + 1)], trajectory[0][0:(i + 1)]) + Xo/codata.c
        trajectory[3][i] = np.trapz(trajectory[6][0:(i + 1)], trajectory[0][0:(i + 1)]) + Zo/codata.c
     #Az
    trajectory[9]=-(trajectory[7]*trajectory[4])/trajectory[6]

    return trajectory



def trajectory_undulator_from_magnetic_field2(By,Z,E,N,nb_enlarg) :
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
    Vo = [0.0, 0.0, Beta * codata.c, 0.0, 0.0, Z[0]]
    Z,By=enlargement_vector_for_interpolation(Z, By, nb_enlarg)
    B = interp1d(Z, By)
    cst=-codata.e/(codata.m_e*gamma)
    res = odeint(fct_ODE,Vo,trajectory[0],args=(cst,B))
    res = np.transpose(res)
    trajectory[4] = res[0]
    trajectory[5] = res[1]
    trajectory[6] = res[2]
    trajectory[1] = res[3]
    trajectory[2] = res[4]
    trajectory[3] = res[5]
    trajectory[7] = -cst * B(trajectory[3]) * trajectory[6]
    trajectory[9] = cst* B(trajectory[3]) * trajectory[4]
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
def undulator_trajectory(K,E,lambda_u,Nb_period=10,Nb_point=100,Z_By=None,type_trajectory=1) :

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

            trajectory = trajectory_undulator_from_magnetic_field1(By=By, Z=Z,K=K,E=E,N=len(By),Vo=0.0,Xo=0.0,Zo=Z[0])

        else :
            trajectory = trajectory_undulator_from_magnetic_field2(By=By, Z=Z,E=E,N=Nb_period * Nb_point + 1
                                                                   ,nb_enlarg=np.floor(len(Z)/10))

    else :
        trajectory = analytical_trajectory_undulator(K=K, lambda_u=lambda_u, Nb_period=Nb_period, Nb_point=Nb_point)

    return trajectory



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

def energy_radiated(omega=2.53465927101*10**17,trajectory=np.zeros((11,10)) , x=0.00 , y=0.0, D=None):
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
    integrand[0] += (-(n_chap[1]**2)*trajectory[7]-n_chap[2]*Alpha)*Alpha2
    integrand[1] += n_chap[1]*(n_chap[0]*trajectory[7]+n_chap[2]*trajectory[9])*Alpha2
    integrand[2] += (-(n_chap[1]**2)*trajectory[9]+n_chap[0]*Alpha)*Alpha2
    integrand *= (1.0 / (1.0 - n_chap[0]*trajectory[4]-n_chap[2]*trajectory[6])) ** 2
    for k in range(3):
        E[k] = np.trapz(integrand[k], trajectory[0])
    return (np.abs(E[0]) ** 2 + np.abs(E[1])** 2 + np.abs(E[2])** 2)





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

    shape2 = res.shape

    for i in range(len(X)):
        res[i] = c6*energy_radiated(omega=omega,trajectory=trajectory , x=X[i] , y=Y[i], D=D )

    X = X.reshape(shape1)
    Y = Y.reshape(shape1)
    res = res.reshape(shape1)
    print("radiation max ")
    print(res.max())
    return res


def energy_radiated2(omega=2.53465927101*10**17,trajectory=np.zeros((11,10)) , x=0.00 , y=0.0, D=None):
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
    Alpha2 = np.exp(0. + 1j * omega * (trajectory[0] - n_chap[0] * trajectory[1] - n_chap[2] * trajectory[3]))
    integrand[0] += (-(n_chap[1]**2)*trajectory[7]-n_chap[2]*Alpha)*Alpha2
    integrand[1] += n_chap[1]*(n_chap[0]*trajectory[7]+n_chap[2]*trajectory[9])*Alpha2
    integrand[2] += (-(n_chap[1]**2)*trajectory[9]+n_chap[0]*Alpha)*Alpha2
    for k in range(3):
        E[k] = np.trapz(integrand[k], trajectory[0])
    return (np.abs(E[0]) ** 2 + np.abs(E[1])** 2 + np.abs(E[2])** 2)





def radiation_single_electron2(K=1.87, E=1.3 * 10 ** 9, lambda_u=0.035, trajectory=np.zeros((11, 10)), D=None, omega=None,
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

    shape2 = res.shape

    for i in range(len(X)):
        res[i] = c6*energy_radiated2(omega=omega,trajectory=trajectory , x=X[i] , y=Y[i], D=D )

    X = X.reshape(shape1)
    Y = Y.reshape(shape1)
    res = res.reshape(shape1)
    print("radiation max ")
    print(res.max())
    return res


########################### essai


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


