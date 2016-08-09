import numpy as np
from pySRU.Radiation import Radiation
import scipy.constants as codata
import scipy.integrate as integrate
from pySRU.TrajectoryFactory import TrajectoryFactory, TRAJECTORY_METHOD_ANALYTIC , TRAJECTORY_METHOD_ODE



RADIATION_METHOD_AUTOMATIC=0
RADIATION_METHOD_NEAR_FIELD=1
RADIATION_METHOD_FARFIELD=2
RADIATION_METHOD_APPROX=3
RADIATION_METHOD_APPROX_FARFIELD=4
eV_to_J=1.602176487e-19


class RadiationFactory(object):
    def __init__(self,method,omega,Nb_pts,formula=None):
        self.omega=omega
        if formula==None :
            formula =1
        self.formula=formula
        self.method=method
        if Nb_pts==None :
            Nb_pts=101
        self.Nb_pts=Nb_pts


    def copy(self):
        return RadiationFactory( method=self.method,omega=self.omega,formula=self.formula,Nb_pts=self.Nb_pts)

    def energy_eV(self):
        return self.omega*codata.hbar/eV_to_J

    # Photon's flow all over a screen situate at distance D of an undulator
    def create_for_one_relativistic_electron(self, trajectory, source, XY_are_list=False, distance=None, X=None, Y=None):
        if X == None or Y == None:
            print('calculate X and Y array')
            if distance == None:
                if self.method == RADIATION_METHOD_APPROX_FARFIELD :
                    xy_max = np.tan(1. / source.Lorentz_factor())
                else:
                    distance = source.choose_distance_automatic(2)
                    xy_max = np.tan(1. / source.Lorentz_factor()) * distance
            else:
                xy_max = np.tan(1. / source.Lorentz_factor()) * distance
            X = np.linspace(0.0, xy_max, self.Nb_pts)
            Y = np.linspace(0.0, xy_max, self.Nb_pts)
        else:
            if type(X) != np.ndarray:
                X = np.linspace(0.0, X, self.Nb_pts)
            if type(Y) != np.ndarray:
                Y = np.linspace(0.0, Y, self.Nb_pts)
                # Nb_pts=len(X)

        if not XY_are_list :
            X, Y = np.meshgrid(X, Y)

        intensity = self.calculate_radiation_intensity(trajectory=trajectory, source=source,
                                                       distance=distance, X_array=X, Y_array=Y)

        radiation = Radiation(intensity=intensity, X=X, Y=Y, distance=distance)

        return radiation



    # Photon's flow all over a screen situate at distance D of an undulator
    def calculate_radiation_intensity(self, trajectory, source, X_array, Y_array, distance=None):
        # c1 = codata.e ** 2 * omega1 ** 2 / (16 * np.pi ** 3 * codata.epsilon_0 * codata.c )
        # c2 = I / codata.e  # multiply by number of electrons in 1 A
        # c3 = 2.0*np.pi / (codata.h * omega1)  # divide by e energy (to get number of photons)
        # c4 = 1e-3 / omega1  # to get 1% energy bandwidth  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        # c5 = 1e-3 * 1e-3  # from rad to .1mrad angle bandwidth
        c6 = codata.e * source.I_current() * 1e-9 / (8.0 * np.pi ** 2 * codata.epsilon_0 * codata.c * codata.h)
        if X_array.size != Y_array.size:
            raise Exception("X and Y dimensions must be equal.")

        res = np.zeros_like(X_array)
        shape1 = res.shape
        X = X_array.flatten()
        Y = Y_array.flatten()
        res = res.flatten()
        gamma = source.Lorentz_factor()
        if self.method == RADIATION_METHOD_APPROX_FARFIELD:
            if self.formula ==1 :
                #print('APPROX AND FARFIELD OK')
                for i in range(len(X)):
                    res[i] = c6 * self.energy_radiated_approximation_and_farfield(trajectory=trajectory,
                                                                                  distance=distance, x=X[i], y=Y[i])
            else :
                for i in range(len(X)):
                    res[i] = c6 * self.energy_radiated_approximation_and_farfield2(trajectory=trajectory,
                                                                                   distance=distance, x=X[i], y=Y[i])
        else:
            if self.method == RADIATION_METHOD_FARFIELD:
                if self.formula == 1:
                    for i in range(len(X)):
                        res[i] = c6 * self.energy_radiated_farfield(trajectory=trajectory, gamma=gamma,
                                                                    distance=distance, x=X[i], y=Y[i])
                else:
                    for i in range(len(X)):
                        res[i] = c6 * self.energy_radiated_farfield2(trajectory=trajectory, distance=distance,
                                                                     gamma=gamma, x=X[i], y=Y[i])
            else:
                if self.method == RADIATION_METHOD_APPROX:
                    if self.formula == 1:
                        for i in range(len(X)):
                            res[i] = c6 * self.energy_radiated_approx(trajectory=trajectory, distance=distance,
                                                                          gamma=gamma, x=X[i], y=Y[i])
                    else:
                        for i in range(len(X)):
                            res[i] = c6 * self.energy_radiated_approx2(trajectory=trajectory, distance=distance,
                                                                           gamma=gamma, x=X[i], y=Y[i])
                else : # Nearfield
                    if self.formula == 1:
                        for i in range(len(X)):
                            res[i] = c6 * self.energy_radiated_near_field(trajectory=trajectory, distance=distance,
                                                                      gamma=gamma, x=X[i], y=Y[i])
                    else:
                        for i in range(len(X)):
                            res[i] = c6 * self.energy_radiated_near_field2(trajectory=trajectory, distance=distance,
                                                                       gamma=gamma, x=X[i], y=Y[i])
        res = res.reshape(shape1)
        return res

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

        E = np.zeros((3,), dtype=np.complex)
        integrand = np.zeros((3, N), dtype=np.complex)
        A1 = (n_chap[0] * trajectory.a_x + n_chap[1] * trajectory.a_y + n_chap[2] * trajectory.a_z)
        A2 = (n_chap[0] * (n_chap[0] - trajectory.v_x) + n_chap[1] * (n_chap[1] - trajectory.v_y)
                + n_chap[2] * (n_chap[2] - trajectory.v_z))
        Alpha2 = np.exp(
            0. + 1j * self.omega * (trajectory.t - n_chap[0] * trajectory.x
                                        - n_chap[1] * trajectory.y - n_chap[2] * trajectory.z))
        Alpha1 = (1.0 / (1.0 - n_chap[0] * trajectory.v_x
                             - n_chap[1] * trajectory.v_y - n_chap[2] * trajectory.v_z)) ** 2

        integrand[0] += (A1 * (n_chap[0] - trajectory.v_x) - A2 * trajectory.a_x) * Alpha2 * Alpha1
        integrand[1] += (A1 * (n_chap[1] - trajectory.v_y) - A2 * trajectory.a_y) * Alpha2 * Alpha1
        integrand[2] += (A1 * (n_chap[2] - trajectory.v_z) - A2 * trajectory.a_z) * Alpha2 * Alpha1
        for k in range(3):
            # E[k] = np.trapz(integrand[k], self.trajectory.t)
            E[k] = np.trapz(integrand[k], trajectory.t)
        return (np.abs(E[0]) ** 2 + np.abs(E[1]) ** 2 + np.abs(E[2]) ** 2)

    def energy_radiated_approximation_and_farfield2(self, trajectory, x, y, distance):
        # N = trajectory.shape[1]
        N = trajectory.nb_points()
        if distance == None:
            # in radian :
            n_chap = np.array([x, y, 1.0 - 0.5 * (x ** 2 + y ** 2)])
        # in meters :
        else:
            R = np.sqrt(x ** 2 + y ** 2 + distance ** 2)
            n_chap = np.array([x, y, distance]) / R

        E = np.zeros((3,), dtype=np.complex)
        integrand = np.zeros((3, N), dtype=np.complex)
        A1 = (n_chap[1] * trajectory.v_z - n_chap[2] * trajectory.v_y)
        A2 = (-n_chap[0] * trajectory.v_z + n_chap[2] * trajectory.v_x)
        A3 = (n_chap[0] * trajectory.v_y - n_chap[1] * trajectory.v_x)
        Alpha2 = np.exp(
            0. + 1j * self.omega * (trajectory.t - n_chap[0] * trajectory.x
                                    - n_chap[1] * trajectory.y - n_chap[2] * trajectory.z))


        integrand[0] += ( n_chap[1]*A3 - n_chap[2]*A2) * Alpha2
        integrand[1] += (- n_chap[0]*A3 + n_chap[2]*A1) * Alpha2
        integrand[2] += ( n_chap[0]*A2 - n_chap[1]*A1) * Alpha2

        for k in range(3):
            # E[k] = np.trapz(integrand[k], self.trajectory.t)
            E[k] = integrate.simps(integrand[k], trajectory.t)
        E *= self.omega * 1j

        Alpha3 = (np.exp(1j*self.omega*(trajectory.t[-1] - n_chap[2] * trajectory.z[-1]))-
                  np.exp(1j*self.omega*(trajectory.t[0] - n_chap[2] * trajectory.z[0])))
        terme_bord = np.full((3,), 0. + 1j * 0., dtype=np.complex)
        terme_bord[0] += n_chap[0] * n_chap[2] * trajectory.v_z[0]
        terme_bord[1] += n_chap[1] * n_chap[2] * trajectory.v_z[0]
        terme_bord[2] += trajectory.v_z[0] * (n_chap[2]**2 - 1.0)
        terme_bord *= Alpha3 / (1.0 - n_chap[2] * trajectory.v_z[0])
        E += terme_bord

        return (np.abs(E[0]) ** 2 + np.abs(E[1]) ** 2 + np.abs(E[2]) ** 2)

    def energy_radiated_farfield(self, trajectory, x, y,gamma, distance):
        N = trajectory.nb_points()
        if distance==None :
            distance=200

        R = np.sqrt(x ** 2 + y ** 2 + distance ** 2)
        n_chap = np.array([x, y, distance]) / R

        E = np.zeros((3,), dtype=np.complex)
        integrand = np.zeros((3, N), dtype=np.complex)
        A1 = (n_chap[0] * trajectory.a_x + n_chap[1] * trajectory.a_y + n_chap[2] * trajectory.a_z)
        A2 = (n_chap[0] * (n_chap[0] - trajectory.v_x) + n_chap[1] * (n_chap[1] - trajectory.v_y)
                + n_chap[2] * (n_chap[2] - trajectory.v_z))
        Alpha2 = np.exp(
            0. + 1j * self.omega * (trajectory.t - n_chap[0] * trajectory.x
                                        - n_chap[1] * trajectory.y - n_chap[2] * trajectory.z))
        Alpha1 = (1.0 / (1.0 - n_chap[0] * trajectory.v_x
                             - n_chap[1] * trajectory.v_y - n_chap[2] * trajectory.v_z)) ** 2
        cst=codata.c/(R*gamma**2)
        integrand[0] += (A1 * (n_chap[0] - trajectory.v_x) - A2 * trajectory.a_x
                         + cst * (n_chap[0] - trajectory.v_x) ) * Alpha2 * Alpha1
        integrand[1] += (A1 * (n_chap[1] - trajectory.v_y) - A2 * trajectory.a_y
                         + cst * (n_chap[1] - trajectory.v_y)  ) * Alpha2 * Alpha1
        integrand[2] += (A1 * (n_chap[2] - trajectory.v_z) - A2 * trajectory.a_z
                         + cst * (n_chap[2] - trajectory.v_z)  ) * Alpha2 * Alpha1
        for k in range(3):
                # E[k] = np.trapz(integrand[k], self.trajectory.t)
            E[k] = np.trapz(integrand[k], trajectory.t)
        return (np.abs(E[0]) ** 2 + np.abs(E[1]) ** 2 + np.abs(E[2]) ** 2)

    def energy_radiated_farfield2(self, trajectory, x, y,gamma, distance):
        # N = trajectory.shape[1]
        N = trajectory.nb_points()
        if distance == None:
            # in radian :
            n_chap = np.array([x, y, 1.0 - 0.5 * (x ** 2 + y ** 2)])
            R=10000.0
        # in meters :
        else:
            R = np.sqrt(x ** 2 + y ** 2 + distance ** 2)
            n_chap = np.array([x, y, distance]) / R

        E = np.zeros((3,), dtype=np.complex)
        integrand = np.zeros((3, N), dtype=np.complex)
        A1 = (n_chap[1] * trajectory.v_z - n_chap[2] * trajectory.v_y)
        A2 = (-n_chap[0] * trajectory.v_z + n_chap[2] * trajectory.v_x)
        A3 = (n_chap[0] * trajectory.v_y - n_chap[1] * trajectory.v_x)
        Alpha1 = (1.0 / (1.0 - n_chap[0] * trajectory.v_x
                         - n_chap[1] * trajectory.v_y - n_chap[2] * trajectory.v_z)) ** 2
        Alpha2 = np.exp(
            0. + 1j * self.omega * (trajectory.t - n_chap[0] * trajectory.x
                                    - n_chap[1] * trajectory.y - n_chap[2] * trajectory.z))
        cst = codata.c / ((gamma ** 2) * R)
        integrand[0] += ((n_chap[1] * A3 - n_chap[2] * A2)*self.omega * 1j
                         + cst * (n_chap[0] - trajectory.v_x) * Alpha1) * Alpha2
        integrand[1] += ((- n_chap[0] * A3 + n_chap[2] * A1)*self.omega * 1j
                         + cst * (n_chap[1] - trajectory.v_y) * Alpha1) * Alpha2
        integrand[2] += ((n_chap[0] * A2 - n_chap[1] * A1)*self.omega * 1j
                         + cst * (n_chap[2] - trajectory.v_z) * Alpha1) * Alpha2

        for k in range(3):
            E[k] = np.trapz(integrand[k], trajectory.t)
            #E[k] = integrate.simps(integrand[k], trajectory.t)


        Alpha3 = (np.exp(1j * self.omega * (trajectory.t[-1] - n_chap[2] * trajectory.z[-1])) -
                  np.exp(1j * self.omega * (trajectory.t[0] - n_chap[2] * trajectory.z[0])))
        terme_bord = np.full((3,), 0. + 1j * 0., dtype=np.complex)
        terme_bord[0] += n_chap[0] * n_chap[2] * trajectory.v_z[0]
        terme_bord[1] += n_chap[1] * n_chap[2] * trajectory.v_z[0]
        terme_bord[2] += trajectory.v_z[0] * (n_chap[2] ** 2 - 1.0)
        terme_bord *= Alpha3 / (1.0 - n_chap[2] * trajectory.v_z[0])
        E += terme_bord

        return (np.abs(E[0]) ** 2 + np.abs(E[1]) ** 2 + np.abs(E[2]) ** 2)

    # exact equation for the energy radiated
    # warning !!!!!!! far field approximation
    def energy_radiated_approx(self,trajectory,gamma, x, y,distance):
        N = trajectory.nb_points()
        n_chap = np.array([x - trajectory.x * codata.c, y - trajectory.y * codata.c, distance - trajectory.z * codata.c])
        R = np.zeros(n_chap.shape[1])
        for i in range(n_chap.shape[1]):
            R[i] = np.linalg.norm(n_chap[:, i])
            n_chap[:, i] /= R[i]

        E = np.zeros((3,), dtype=np.complex)
        integrand = np.zeros((3, N), dtype=np.complex)
        A1 = (n_chap[0] * trajectory.a_x + n_chap[1] * trajectory.a_y + n_chap[2] * trajectory.a_z)
        A2 = (n_chap[0] * (n_chap[0] - trajectory.v_x) + n_chap[1] * (n_chap[1] - trajectory.v_y)
              + n_chap[2] * (n_chap[2] - trajectory.v_z))
        Alpha2 = np.exp(
            0. + 1j * self.omega * (trajectory.t - n_chap[0] * trajectory.x
                                    - n_chap[1] * trajectory.y - n_chap[2] * trajectory.z))
        Alpha1 = (1.0 / (1.0 - n_chap[0] * trajectory.v_x
                         - n_chap[1] * trajectory.v_y - n_chap[2] * trajectory.v_z)) ** 2
        integrand[0] += ((A1 * (n_chap[0] - trajectory.v_x) - A2 * trajectory.a_x)
                         ) * Alpha2 * Alpha1
        integrand[1] += ((A1 * (n_chap[1] - trajectory.v_y) - A2 * trajectory.a_y)
                         ) * Alpha2 * Alpha1
        integrand[2] += ((A1 * (n_chap[2] - trajectory.v_z) - A2 * trajectory.a_z)
                         ) * Alpha2 * Alpha1
        for k in range(3):
            E[k] = np.trapz(integrand[k], trajectory.t)
            #E[k] = integrate.simps(integrand[k], trajectory.t)
        return (np.abs(E[0]) ** 2 + np.abs(E[1]) ** 2 + np.abs(E[2]) ** 2)

    def energy_radiated_approx2(self, trajectory, gamma, x, y, distance):
        N = trajectory.nb_points()
        n_chap = np.array([x - trajectory.x * codata.c, y - trajectory.y * codata.c, distance - trajectory.z * codata.c])
        R = np.zeros(n_chap.shape[1])
        for i in range(n_chap.shape[1]):
            R[i] = np.linalg.norm(n_chap[:, i])
            n_chap[:, i] /= R[i]

        E = np.zeros((3,), dtype=np.complex)
        integrand = np.zeros((3, N), dtype=np.complex)

        A1 = (n_chap[1] * trajectory.v_z - n_chap[2] * trajectory.v_y)
        A2 = (-n_chap[0] * trajectory.v_z + n_chap[2] * trajectory.v_x)
        A3 = (n_chap[0] * trajectory.v_y - n_chap[1] * trajectory.v_x)
        Alpha2 = np.exp( 0. + 1j * self.omega * (trajectory.t +R/codata.c))

        integrand[0] += (n_chap[1] * A3 - n_chap[2] * A2) * Alpha2
        integrand[1] += (- n_chap[0] * A3 + n_chap[2] * A1) * Alpha2
        integrand[2] += (n_chap[0] * A2 - n_chap[1] * A1) * Alpha2
        for k in range(3):
            E[k] = np.trapz(integrand[k], trajectory.t)
            #E[k] = integrate.simps(integrand[k], trajectory.t)
        E *= self.omega * 1j

        Alpha3 = (np.exp(1j * self.omega * (trajectory.t[-1] - n_chap[2][-1] * trajectory.z[-1])) -
                  np.exp(1j * self.omega * (trajectory.t[0] - n_chap[2][0] * trajectory.z[0])))
        terme_bord = np.full((3,), 0. + 1j * 0., dtype=np.complex)
        terme_bord[0] += (n_chap[0][0] * n_chap[2][0] * trajectory.v_z[0]
                            )*Alpha3 / (1.0 - n_chap[2][0] * trajectory.v_z[0])
        terme_bord[1] += (n_chap[1][0] * n_chap[2][0] * trajectory.v_z[0]
                    )*Alpha3 / (1.0 - n_chap[2][0] * trajectory.v_z[0])
        terme_bord[2] += (trajectory.v_z[0] * (n_chap[2][0] ** 2 - 1.0)
                        )*Alpha3 / (1.0 - n_chap[2][0] * trajectory.v_z[0])
        E += terme_bord

        return (np.abs(E[0]) ** 2 + np.abs(E[1]) ** 2 + np.abs(E[2]) ** 2)

    # energy radiated without the the far filed approxiamation
    def energy_radiated_near_field(self,gamma,trajectory, x, y,distance):
        N = trajectory.nb_points()
        n_chap = np.array([x - trajectory.x * codata.c, y - trajectory.y * codata.c, distance - trajectory.z * codata.c])
        R = np.zeros(n_chap.shape[1])
        for i in range(n_chap.shape[1]):
            R[i] = np.linalg.norm(n_chap[:, i])
            n_chap[:, i] /= R[i]

        E = np.zeros((3,), dtype=np.complex)
        integrand = np.zeros((3, N), dtype=np.complex)
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
            E[k] = np.trapz(integrand[k], trajectory.t)
        return (np.abs(E[0]) ** 2 + np.abs(E[1]) ** 2 + np.abs(E[2]) ** 2)

    def energy_radiated_near_field2(self, trajectory, gamma,x, y, distance):
        # N = trajectory.shape[1]
        N = trajectory.nb_points()
        n_chap = np.array(
            [x - trajectory.x * codata.c, y - trajectory.y * codata.c, distance - trajectory.z * codata.c])
        R = np.zeros(n_chap.shape[1])
        for i in range(n_chap.shape[1]):
            R[i] = np.linalg.norm(n_chap[:, i])
            n_chap[:, i] /= R[i]

        E = np.zeros((3,), dtype=np.complex)
        integrand = np.zeros((3, N), dtype=np.complex)

        A1 = (n_chap[1] * trajectory.v_z - n_chap[2] * trajectory.v_y)
        A2 = (-n_chap[0] * trajectory.v_z + n_chap[2] * trajectory.v_x)
        A3 = (n_chap[0] * trajectory.v_y - n_chap[1] * trajectory.v_x)
        Alpha1 = np.exp(
            0. + 1j * self.omega * (trajectory.t - n_chap[0] * trajectory.x
                                    - n_chap[1] * trajectory.y - n_chap[2] * trajectory.z))
        Alpha2 = codata.c / (gamma ** 2 * R)
        integrand[0] += ((n_chap[1] * A3 - n_chap[2] * A2) * self.omega * 1j
                         + Alpha2 * (n_chap[0] - trajectory.v_x)
                         ) * Alpha1
        integrand[1] += ((-n_chap[0] * A3 + n_chap[2] * A1) * self.omega * 1j
                         + Alpha2 * (n_chap[1] - trajectory.v_y)
                         ) * Alpha1
        integrand[2] += ((n_chap[0] * A2 - n_chap[1] * A1) * self.omega * 1j
                         + Alpha2 * (n_chap[2] - trajectory.v_z)
                         ) * Alpha1
        for k in range(3):
            # E[k] = np.trapz(integrand[k], self.trajectory.t)
            E[k] = integrate.simps(integrand[k], trajectory.t)

        terme_bord = np.full((3), 0. + 1j * 0., dtype=np.complex)
        unit_z = np.array([0.0, 0.0, 1.0])
        terme_bord += ((n_chap[2][-1] * n_chap[:, -1] - unit_z) *
                       np.exp(1j * self.omega * (trajectory.t[-1] - n_chap[2][-1] * trajectory.z[-1])) / (
                           1.0 - n_chap[2][-1] * trajectory.v_z[-1]
                       ))
        terme_bord -= ((n_chap[2][0] * n_chap[:, 0] - unit_z) *
                       np.exp(1j * self.omega * (trajectory.t[-1] - n_chap[2][0] * trajectory.z[0])) / (
                           1.0 - n_chap[2][0] * trajectory.v_z[0]
                       ))
        terme_bord *= trajectory.v_z[0]
        E += terme_bord

        return (np.abs(E[0]) ** 2 + np.abs(E[1]) ** 2 + np.abs(E[2]) ** 2)

    def get_method(self):
        if self.method == RADIATION_METHOD_NEAR_FIELD:
            method = 'Near field calculation'

        elif self.method == RADIATION_METHOD_APPROX:
            method = ' Near field and approximate calculation '

        elif self.method == RADIATION_METHOD_FARFIELD:
            method = ' Far field calculation '

        elif self.method == RADIATION_METHOD_APPROX_FARFIELD:
            method = '  Far field and approximation calculation '
        else :
            method=' Problem unknow method'
        return method

    def print_parameters(self):
        print('Radiation ')
        print('    method : %s' %self.get_method())
        print('    number of point in each direction : %d' %self.Nb_pts)
        print('    energy of the emission  (eV): %f' %self.energy_eV())


def Exemple_FARFIELD():
    from MagneticStructureUndulatorPlane import MagneticStructureUndulatorPlane as Undulator
    from ElectronBeam import ElectronBeam
    from SourceUndulatorPlane import SourceUndulatorPlane

    undulator_test = Undulator(K=1.87, period_length=0.035, length=0.035 * 14)
    electron_beam_test = ElectronBeam(Electron_energy=1.3, I_current=1.0)
    magnetic_field_test = undulator_test.create_magnetic_field()
    source_test = SourceUndulatorPlane(undulator=undulator_test,
                         electron_beam=electron_beam_test,
                         magnetic_field=magnetic_field_test)

    traj = TrajectoryFactory(Nb_pts=2000, method=TRAJECTORY_METHOD_ANALYTIC).create_from_source(source_test)

    Rad = RadiationFactory(omega=source_test.harmonic_frequency(1), method=RADIATION_METHOD_APPROX_FARFIELD, Nb_pts=101
                            ).create_for_one_relativistic_electron(trajectory=traj, source=source_test)

    print('Screen distance :')
    print(Rad.distance)

    print("screen shape ")
    print(Rad.intensity.shape)

    print('X max :')
    print(Rad.X.max())

    print('Y max :')
    print(Rad.Y.max())

    print('intensity max ()')
    print(Rad.max())

    print('plot')
    Rad.plot()


def Exemple_NEARFIELD():
    from MagneticStructureUndulatorPlane import MagneticStructureUndulatorPlane as Undulator
    from ElectronBeam import ElectronBeam
    from SourceUndulatorPlane import SourceUndulatorPlane

    undulator_test = Undulator(K=1.87, period_length=0.035, length=0.035 * 14)
    electron_beam_test = ElectronBeam(Electron_energy=1.3, I_current=1.0)
    magnetic_field_test = undulator_test.create_magnetic_field()
    source_test = SourceUndulatorPlane(undulator=undulator_test,
                         electron_beam=electron_beam_test,
                         magnetic_field=magnetic_field_test)

    traj = TrajectoryFactory(Nb_pts=2000, method=TRAJECTORY_METHOD_ODE).create_from_source(source_test)

    Rad = RadiationFactory(omega=source_test.harmonic_frequency(1), method=RADIATION_METHOD_NEAR_FIELD, Nb_pts=101
                           ).create_for_one_relativistic_electron(trajectory=traj, source=source_test)

    print('Screen distance :')
    print(Rad.distance)

    print("screen shape ")
    print(Rad.intensity.shape)

    print('X max :')
    print(Rad.X.max())

    print('Y max :')
    print(Rad.Y.max())

    print('intensity max ()')
    print(Rad.max())

    print('plot')
    Rad.plot()

if __name__ == "__main__" :
    pass
    Exemple_FARFIELD()
    #Exemple_NEARFIELD()

