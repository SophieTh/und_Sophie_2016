import numpy as np
from pylab import *
import scipy.constants as codata
import scipy.integrate as integrate
from scipy.integrate import odeint
from scipy.interpolate import interp1d
import Trajectory
from Radiation import Radiation

RADIATION_METHOD_AUTOMATIC=0
RADIATION_METHOD_FARFIELD=1
RADIATION_METHOD_APPROX=2
RADIATION_METHOD_NO_FARFIELD=3


class UndulatorRadiationFactory(object):
    def __init__(self,Trajectory,method,omega=None):
        self.trajectory=Trajectory
        self.omega=omega
        self.method=method

    # Photon's flow all over a screen situate at distance D of an undulator
    def calculate_radiation_intensity(self, distance, X_arrays=None, Y_arrays=None):

        gamma = self.trajectory.parameters.E / 0.511e6

        if self.omega == None:  # use the first harmonic at the resonance
            self.omega = ((2.0 * gamma ** 2) / (1.0 + (self.trajectory.parameters.K ** 2) / 2.0)) * (
            (2.0 * np.pi * codata.c)
            / self.trajectory.parameters.lambda_u)
        # c1 = codata.e ** 2 * omega1 ** 2 / (16 * np.pi ** 3 * codata.epsilon_0 * codata.c )
        # c2 = 1.0 / codata.e  # multiply by number of electrons in 1 A
        # c3 = 2.0*np.pi / (codata.h * omega1)  # divide by e energy (to get number of photons)
        # c4 = 1e-2 / omega1  # to get 1% energy bandwidth  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        # c5 = 0.1e-3 * 0.1e-3  # from rad to .1mrad angle bandwidth
        c6 = codata.e * 1e-10 / (8.0 * np.pi ** 2 * codata.epsilon_0 * codata.c * codata.h)

        if X_arrays== None or Y_arrays == None:
            print('ok')
            if distance == None:
                distance = 100
            X = np.arange(0.0, (distance * 1.01) * 1e-3, distance * 1e-5)
            Y = np.arange(0.0, (distance * 1.01) * 1e-3, distance * 1e-5)
        else :
            X=X_arrays.copy()
            Y=Y_arrays.copy()
        X, Y = meshgrid(X, Y)

        if X.size != Y.size:
            raise Exception("X and Y dimensions must be equal.")

        res = np.zeros_like(X)
        shape1 = res.shape

        X = X.flatten()
        Y = Y.flatten()
        res = res.flatten()

        if self.method == RADIATION_METHOD_FARFIELD:
            for i in range(len(X)):
                res[i] = c6 * self.energy_radiated_fast_and_farfield(distance=distance, x=X[i], y=Y[i])
        else:
            if self.method == RADIATION_METHOD_APPROX:
                for i in range(len(X)):
                    res[i] = c6 * self.energy_radiated_approxiamation(distance=distance,
                                                                      gamma=gamma, x=X[i], y=Y[i])
            else:
                if self.method == RADIATION_METHOD_NO_FARFIELD:
                    for i in range(len(X)):
                        res[i] = c6 * self.energy_radiated_no_far_field(distance=distance,
                                                                        gamma=gamma, x=X[i], y=Y[i])
                else:
                    print('coucou')
        res = res.reshape(shape1)
        return res

    # Photon's flow all over a screen situate at distance D of an undulator
    def create_for_single_electron(self, distance,X=None,Y=None):
        map=self.calculate_radiation_intensity(distance=distance,X_arrays=X,Y_arrays=Y)
        radiation= Radiation(map,X,Y,distance,self)
        return radiation

        # approximation of the energy radiated by an electron in a PLANE undulator
        # warning : trajectory is the trajectory difine like the function "undulator trajectory" before :
    def energy_radiated_fast_and_farfield(self, x, y,distance):
            #N = trajectory.shape[1]
            N = self.trajectory.parameters.Nb_pts
            # if self.distance == None:
            #     # in radian :
            #     n_chap = np.array([x, y, 1.0 - 0.5 * (x ** 2 + y ** 2)])
            # in meters :
            R = np.sqrt(x ** 2 + y ** 2 + distance ** 2)
            n_chap = np.array([x, y, distance]) / R

            E = np.full((3,), 0. + 1j * 0., dtype=np.complex)
            integrand = np.full((3, N), 0. + 1j * 0., dtype=np.complex)
            Alpha = self.trajectory.a_x * (n_chap[2] - self.trajectory.v_z) - (n_chap[0] - self.trajectory.v_x) * self.trajectory.a_z
            Alpha2 = np.exp(0. + 1j * self.omega * (self.trajectory.t - n_chap[0] * self.trajectory.x - n_chap[2] * self.trajectory.z))
            Alpha1 = (1.0 / (1.0 - n_chap[0] * self.trajectory.v_x - n_chap[2] * self.trajectory.v_z)) ** 2
            integrand[0] += (-(n_chap[1] ** 2) * self.trajectory.a_x - n_chap[2] * Alpha) * Alpha2 * Alpha1
            integrand[1] += n_chap[1] * (n_chap[0] * self.trajectory.a_x + n_chap[2] * self.trajectory.a_z
                                         ) * Alpha2 * Alpha1
            integrand[2] += (-(n_chap[1] ** 2) * self.trajectory.a_z + n_chap[0] * Alpha) * Alpha2 * Alpha1
            for k in range(3):
                #E[k] = np.trapz(integrand[k], self.trajectory.t)
                E[k] = integrate.simps(integrand[k], self.trajectory.t)
            return (np.abs(E[0]) ** 2 + np.abs(E[1]) ** 2 + np.abs(E[2]) ** 2)

    # exact equation for the energy radiated
    # warning !!!!!!! far field approximation
    def energy_radiated_approxiamation(self,gamma, x, y,distance):
        N = self.trajectory.parameters.Nb_pts
        # in meters :
        n_chap = np.array([x, y, distance])
        R = norm(n_chap)
        n_chap /= R

        E = np.full((3,), 0. + 1j * 0., dtype=np.complex)
        integrand = np.full((3, N), 0. + 1j * 0., dtype=np.complex)
        Alpha = self.trajectory.a_x * (n_chap[2] - self.trajectory.v_z) - (n_chap[0] - self.trajectory.v_x
                                                                           ) * self.trajectory.a_z
        Alpha2 = np.exp(0. + 1j * self.omega * (self.trajectory.t - n_chap[0] * self.trajectory.x
                                           - n_chap[2] * self.trajectory.z))
        Alpha1 = (1.0 / (1.0 - n_chap[0] * self.trajectory.v_x - n_chap[2] * self.trajectory.v_z)) ** 2
        Alpha3 = codata.c / (gamma ** 2 * R)
        integrand[0] += (-(n_chap[1] ** 2) * self.trajectory.a_x - n_chap[2] * Alpha
                         + Alpha3 * (n_chap[0] - self.trajectory.v_x)) * Alpha2 * Alpha1
        integrand[1] += (n_chap[1] * (n_chap[0] * self.trajectory.a_x + n_chap[2] * self.trajectory.a_z)
                         + Alpha3 * (n_chap[1] - self.trajectory.v_y)) * Alpha2 * Alpha1
        integrand[2] += (-(n_chap[1] ** 2) * self.trajectory.a_z + n_chap[0] * Alpha
                         + Alpha3 * (n_chap[2] - self.trajectory.v_z)) * Alpha2 * Alpha1
        for k in range(3):
            # E[k] = np.trapz(integrand[k], trajectory[0])
            E[k] = integrate.simps(integrand[k], self.trajectory.t)
        return (np.abs(E[0]) ** 2 + np.abs(E[1]) ** 2 + np.abs(E[2]) ** 2)

    # energy radiated without the the far filed approxiamation
    def energy_radiated_no_far_field(self,gamma, x, y,distance):
        #N = trajectory.shape[1]
        N = self.trajectory.parameters.Nb_pts
        n_chap = np.array([x - self.trajectory.x * codata.c, y - self.trajectory.y * codata.c, distance - self.trajectory.z * codata.c])
        R = np.zeros(n_chap.shape[1])
        for i in range(n_chap.shape[1]):
            R[i] = norm(n_chap[:, i])
            n_chap[:, i] /= R[i]

        E = np.full((3,), 0. + 1j * 0., dtype=np.complex)
        integrand = np.full((3, N), 0. + 1j * 0., dtype=np.complex)
        Alpha = self.trajectory.a_x * (n_chap[2] - self.trajectory.v_z) - (n_chap[0] - self.trajectory.v_x
                                                                           ) * self.trajectory.a_z
        Alpha2 = np.exp(0. + 1j * self.omega * (self.trajectory.t + R / codata.c))
        Alpha1 = (1.0 / (1.0 - n_chap[0] * self.trajectory.v_x - n_chap[2] * self.trajectory.v_z)) ** 2
        Alpha3 = codata.c / (gamma ** 2 * R)
        integrand[0] += (-(n_chap[1] ** 2) * self.trajectory.a_x - n_chap[2] * Alpha
                         + Alpha3 * (n_chap[0] - self.trajectory.v_x)
                         ) * Alpha2 * Alpha1
        integrand[1] += (n_chap[1] * (n_chap[0] * self.trajectory.a_x + n_chap[2] * self.trajectory.a_z)
                         + Alpha3 * (n_chap[1] - self.trajectory.v_y)
                         ) * Alpha2 * Alpha1
        integrand[2] += (-(n_chap[1] ** 2) * self.trajectory.a_z + n_chap[0] * Alpha
                         + Alpha3 * (n_chap[2] - self.trajectory.v_z)
                         ) * Alpha2 * Alpha1
        for k in range(3):
            # E[k] = np.trapz(integrand[k], trajectory[0])
            E[k] = integrate.simps(integrand[k], self.trajectory.t)
        return (np.abs(E[0]) ** 2 + np.abs(E[1]) ** 2 + np.abs(E[2]) ** 2)


