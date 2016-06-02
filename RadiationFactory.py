import numpy as np
import scipy.constants as codata
import scipy.integrate as integrate
from Radiation import Radiation

RADIATION_METHOD_AUTOMATIC=0
RADIATION_METHOD_NEAR_FIELD=1
RADIATION_METHOD_FARFIELD=2
RADIATION_METHOD_APPROX_FARFIELD=3


class RadiationFactory(object):
    def __init__(self,method,omega):
        self.omega=omega
        self.method=method

    def copy(self):
        return RadiationFactory( method=self.method,omega=self.omega)


    # Photon's flow all over a screen situate at distance D of an undulator
    def calculate_radiation_intensity(self, trajectory,undulator,distance=None, X_arrays=None, Y_arrays=None):
        # c1 = codata.e ** 2 * omega1 ** 2 / (16 * np.pi ** 3 * codata.epsilon_0 * codata.c )
        # c2 = I / codata.e  # multiply by number of electrons in 1 A
        # c3 = 2.0*np.pi / (codata.h * omega1)  # divide by e energy (to get number of photons)
        # c4 = 1e-2 / omega1  # to get 1% energy bandwidth  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        # c5 = 0.1e-3 * 0.1e-3  # from rad to .1mrad angle bandwidth
        c6 = codata.e *undulator.I* 1e-10 / (8.0 * np.pi ** 2 * codata.epsilon_0 * codata.c * codata.h)
        if X_arrays== None or Y_arrays == None:
            print('ok')
            if distance == None:
                distance = 100
            X = np.arange(0.0, (distance * 1.01) * 1e-3, distance * 1e-5)
            Y = np.arange(0.0, (distance * 1.01) * 1e-3, distance * 1e-5)
        else :
            X=X_arrays.copy()
            Y=Y_arrays.copy()
        X, Y = np.meshgrid(X, Y)

        if X.size != Y.size:
            raise Exception("X and Y dimensions must be equal.")

        res = np.zeros_like(X)
        shape1 = res.shape

        X = X.flatten()
        Y = Y.flatten()
        res = res.flatten()

        if self.method == RADIATION_METHOD_APPROX_FARFIELD:
            for i in range(len(X)):
                res[i] = c6 * self.energy_radiated_approximation_and_farfield(trajectory=trajectory,distance=distance, x=X[i], y=Y[i])
        else:
            if self.method == RADIATION_METHOD_FARFIELD:
                for i in range(len(X)):
                    res[i] = c6 * self.energy_radiated_farfield(trajectory=trajectory,distance=distance,
                                                                      x=X[i], y=Y[i])
            else:
                if self.method == RADIATION_METHOD_NEAR_FIELD:
                    for i in range(len(X)):
                        res[i] = c6 * self.energy_radiated_near_field(trajectory=trajectory,distance=distance,
                                                                         x=X[i], y=Y[i])
                else:
                    print('coucou')
        res = res.reshape(shape1)
        return res

    # Photon's flow all over a screen situate at distance D of an undulator
    def create_for_single_electron(self, trajectory,undulator,distance,X=None,Y=None):
        map=self.calculate_radiation_intensity(trajectory=trajectory,undulator=undulator,
                                               distance=distance,X_arrays=X,Y_arrays=Y)
        radiation= Radiation(map=map,X=X,Y=Y,distance=distance)
        return radiation

    # approximation of the energy radiated by an electron in a PLANE undulator
    # warning : trajectory is the trajectory difine like the function "undulator trajectory" before :
    def energy_radiated_approximation_and_farfield(self,trajectory, x, y, distance):
        # N = trajectory.shape[1]
        N = trajectory.nb_points()
        if distance == None:
            # in radian :
            n_chap = np.array([x, y, 1.0 - 0.5 * (x ** 2 + y ** 2)])
        #in meters :
        else :
            R = np.sqrt(x ** 2 + y ** 2 + distance ** 2)
            n_chap = np.array([x, y, distance]) / R

        E = np.full((3,), 0. + 1j * 0., dtype=np.complex)
        integrand = np.full((3, N), 0. + 1j * 0., dtype=np.complex)
        A1 =(n_chap[0]*trajectory.a_x +n_chap[1]*trajectory.a_y +n_chap[2]*trajectory.a_z )
        A2=(n_chap[0]*(n_chap[0]-trajectory.v_x)+n_chap[1]*(n_chap[1]-trajectory.v_y)
            +n_chap[2]*(n_chap[2]-trajectory.v_z))
        Alpha2 = np.exp(
            0. + 1j * self.omega * (trajectory.t - n_chap[0] * trajectory.x
                                    -n_chap[1] * trajectory.y- n_chap[2] * trajectory.z))
        Alpha1 = (1.0 / (1.0 - n_chap[0] * trajectory.v_x
                         -n_chap[1] * trajectory.v_y- n_chap[2] * trajectory.v_z)) ** 2

        integrand[0] += (A1*(n_chap[0]-trajectory.v_x)-A2*trajectory.a_x) * Alpha2 * Alpha1
        integrand[1] += (A1 * (n_chap[1] - trajectory.v_y) - A2 * trajectory.a_y) * Alpha2 * Alpha1
        integrand[2] += (A1*(n_chap[2]-trajectory.v_z)-A2*trajectory.a_z) * Alpha2 * Alpha1

        for k in range(3):
            # E[k] = np.trapz(integrand[k], self.trajectory.t)
            E[k] = integrate.simps(integrand[k], trajectory.t)
        return (np.abs(E[0]) ** 2 + np.abs(E[1]) ** 2 + np.abs(E[2]) ** 2)

    # exact equation for the energy radiated
    # warning !!!!!!! far field approximation
    def energy_radiated_farfield(self,trajectory, x, y,distance):
        N = trajectory.nb_points()
        gamma=self.undulator.E/0.511e6
        # in meters :
        n_chap = np.array([x, y, distance])
        R = np.norm(n_chap)
        n_chap /= R

        E = np.full((3,), 0. + 1j * 0., dtype=np.complex)
        integrand = np.full((3, N), 0. + 1j * 0., dtype=np.complex)
        A1 = (n_chap[0] * trajectory.a_x + n_chap[1] * trajectory.a_y + n_chap[2] * trajectory.a_z)
        A2 = (n_chap[0] * (n_chap[0] - trajectory.v_x) + n_chap[1] * (n_chap[1] - trajectory.v_y)
              + n_chap[2] * (n_chap[2] - trajectory.v_z))
        Alpha2 = np.exp(
            0. + 1j * self.omega * (trajectory.t - n_chap[0] * trajectory.x
                                    -n_chap[1] * trajectory.y- n_chap[2] * trajectory.z))
        Alpha1 = (1.0 / (1.0 - n_chap[0] * trajectory.v_x
                         -n_chap[1] * trajectory.v_y- n_chap[2] * trajectory.v_z)) ** 2
        Alpha3 = codata.c / (gamma ** 2 * R)
        integrand[0] += ((A1*(n_chap[0]-trajectory.v_x)-A2*trajectory.a_x)
                         + Alpha3 * (n_chap[0] - trajectory.v_x)) * Alpha2 * Alpha1
        integrand[1] += ((A1 * (n_chap[1] - trajectory.v_y) - A2 * trajectory.a_y)
                         + Alpha3 * (n_chap[1] - trajectory.v_y)) * Alpha2 * Alpha1
        integrand[2] += ((A1*(n_chap[2]-trajectory.v_z)-A2*trajectory.a_z)
                         + Alpha3 * (n_chap[2] - trajectory.v_z)) * Alpha2 * Alpha1
        for k in range(3):
            # E[k] = np.trapz(integrand[k], trajectory[0])
            E[k] = integrate.simps(integrand[k], trajectory.t)
        return (np.abs(E[0]) ** 2 + np.abs(E[1]) ** 2 + np.abs(E[2]) ** 2)

    # energy radiated without the the far filed approxiamation
    def energy_radiated_near_field(self,undulator,trajectory, x, y,distance):
        N = trajectory.nb_points()
        gamma = undulator.E / 0.511e6
        n_chap = np.array([x - trajectory.x * codata.c, y - trajectory.y * codata.c, distance - trajectory.z * codata.c])
        R = np.zeros(n_chap.shape[1])
        for i in range(n_chap.shape[1]):
            R[i] = np.norm(n_chap[:, i])
            n_chap[:, i] /= R[i]

        E = np.full((3,), 0. + 1j * 0., dtype=np.complex)
        integrand = np.full((3, N), 0. + 1j * 0., dtype=np.complex)
        A1 = (n_chap[0] * trajectory.a_x + n_chap[1] * trajectory.a_y + n_chap[2] * trajectory.a_z)
        A2 = (n_chap[0] * (n_chap[0] - trajectory.v_x) + n_chap[1] * (n_chap[1] - trajectory.v_y)
              + n_chap[2] * (n_chap[2] - trajectory.v_z))
        Alpha2 = np.exp(
            0. + 1j * self.omega * (trajectory.t - n_chap[0] * trajectory.x
                                    -n_chap[1] * trajectory.y- n_chap[2] * trajectory.z))
        Alpha1 = (1.0 / (1.0 - n_chap[0] * trajectory.v_x
                         -n_chap[1] * trajectory.v_y- n_chap[2] * trajectory.v_z)) ** 2
        Alpha3 = codata.c / (gamma ** 2 * R)
        integrand[0] += ((A1 * (n_chap[0] - trajectory.v_x) - A2 * trajectory.a_x)
                         + Alpha3 * (n_chap[0] - trajectory.v_x)
                         ) * Alpha2 * Alpha1
        integrand[1] += ((A1 * (n_chap[1] - trajectory.v_y) - A2 * trajectory.a_y)
                         + Alpha3 * (n_chap[1] - trajectory.v_y)
                         ) * Alpha2 * Alpha1
        integrand[2] += ((A1 * (n_chap[2] - trajectory.v_z) - A2 * trajectory.a_z)
                         + Alpha3 * (n_chap[2] - trajectory.v_z)
                         ) * Alpha2 * Alpha1
        for k in range(3):
            # E[k] = np.trapz(integrand[k], trajectory[0])
            E[k] = integrate.simps(integrand[k], trajectory.t)
        return (np.abs(E[0]) ** 2 + np.abs(E[1]) ** 2 + np.abs(E[2]) ** 2)




