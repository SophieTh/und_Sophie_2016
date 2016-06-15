import numpy as np
import time
import scipy.constants as codata
import scipy.integrate as integrate
from Radiation import Radiation
from TrajectoryFactory import TrajectoryFactory, TRAJECTORY_METHOD_ANALYTIC
from UndulatorParameter import UndulatorParameters as Undulator

RADIATION_METHOD_AUTOMATIC=0
RADIATION_METHOD_NEAR_FIELD=1
RADIATION_METHOD_FARFIELD=2
RADIATION_METHOD_APPROX_FARFIELD=3


class RadiationFactory(object):
    def __init__(self,method,omega,Nb_pts):
        self.omega=omega
        self.method=method
        self.Nb_pts=Nb_pts


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
            print('calculate X and Y array')
            if distance == None:
                distance = 100
            #xy_max=np.sqrt(undulator.lambda_u*2.0*undulator.L)/(2.0*np.pi)
            xy_max = np.tan(1.0/undulator.gamma())*distance
            X_arrays = np.linspace(0.0, xy_max, self.Nb_pts)
            Y_arrays = np.linspace(0.0, xy_max, self.Nb_pts)

        X = X_arrays.copy()
        Y = Y_arrays.copy()
        X, Y = np.meshgrid(X, Y)

        if X.size != Y.size:
            raise Exception("X and Y dimensions must be equal.")

        res = np.zeros_like(X)
        shape1 = res.shape

        X = X.flatten()
        Y = Y.flatten()
        res = res.flatten()
        gamma=undulator.E/0.511e6
        if self.method == RADIATION_METHOD_APPROX_FARFIELD:
            for i in range(len(X)):
                # res[i] = c6 * self.energy_radiated_approximation_and_farfield(trajectory=trajectory,
                #                                                                distance=distance, x=X[i], y=Y[i])
                res[i] = c6 * self.energy_radiated_approximation_and_farfield(trajectory=trajectory,
                                                                      distance=distance, x=X[i], y=Y[i])
        else:
            if self.method == RADIATION_METHOD_FARFIELD:
                for i in range(len(X)):
                    res[i] = c6 * self.energy_radiated_farfield(trajectory=trajectory,distance=distance,
                                                                gamma=gamma,  x=X[i], y=Y[i])
            else:
                if self.method == RADIATION_METHOD_NEAR_FIELD:
                    print("near field method")
                    for i in range(len(X)):
                        res[i] = c6 * self.energy_radiated_near_field(trajectory=trajectory,distance=distance,
                                                                      gamma=gamma, x=X[i], y=Y[i])
                else:
                    print('coucou')
        res = res.reshape(shape1)
        return res

    # Photon's flow all over a screen situate at distance D of an undulator
    def create_for_single_electron(self, trajectory,undulator,distance=None,X=None,Y=None):
        if X== None or Y == None:
            print('create X and Y array')
            if distance == None:
                distance = 100
            xy_max=np.sqrt(undulator.lambda_u*2.0*undulator.L)/(2.0*np.pi)
            X = np.linspace(0.0, xy_max, self.Nb_pts)
            Y = np.linspace(0.0, xy_max, self.Nb_pts)

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

        E = np.zeros((3,), dtype=np.complex)
        integrand = np.zeros((3, N), dtype=np.complex)
        A1 =(n_chap[0]*trajectory.a_x +n_chap[1]*trajectory.a_y +n_chap[2]*trajectory.a_z )
        A2=(n_chap[0]*(n_chap[0]-trajectory.v_x)+n_chap[1]*(n_chap[1]-trajectory.v_y)
            +n_chap[2]*(n_chap[2]-trajectory.v_z))
        Alpha2 = np.exp(
            0. + 1j * self.omega * (trajectory.t - n_chap[0] * trajectory.x
                                    -n_chap[1] * trajectory.y- n_chap[2] * trajectory.z))
        Alpha1 = (1.0 / (1.0 - n_chap[0] * trajectory.v_x
                         -n_chap[1] * trajectory.v_y- n_chap[2] * trajectory.v_z)) ** 2

        integrand[0] += (A1*(n_chap[0]-trajectory.v_x)-A2*trajectory.a_x) * Alpha2*Alpha1
        integrand[1] += (A1 * (n_chap[1] - trajectory.v_y) - A2 * trajectory.a_y) * Alpha2*Alpha1
        integrand[2] += (A1*(n_chap[2]-trajectory.v_z)-A2*trajectory.a_z) * Alpha2*Alpha1

        for k in range(3):
            # E[k] = np.trapz(integrand[k], self.trajectory.t)
            E[k] = integrate.simps(integrand[k], trajectory.t)
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

    # exact equation for the energy radiated
    # warning !!!!!!! far field approximation
    def energy_radiated_farfield(self,trajectory,gamma, x, y,distance):
        N = trajectory.nb_points()
        # in meters :
        n_chap = np.array([x, y, distance])
        R = np.sqrt(x ** 2 + y ** 2 + distance ** 2)
        n_chap /= R

        E = np.zeros((3,), dtype=np.complex)
        integrand = np.zeros((3, N), dtype=np.complex)
        A1 = (n_chap[0] * trajectory.a_x + n_chap[1] * trajectory.a_y + n_chap[2] * trajectory.a_z)
        A2 = (n_chap[0] * (n_chap[0] - trajectory.v_x) + n_chap[1] * (n_chap[1] - trajectory.v_y)
              + n_chap[2] * (n_chap[2] - trajectory.v_z))
        Alpha2 = np.exp(
            0. + 1j * self.omega * (trajectory.t - n_chap[0] * trajectory.x
                                    -n_chap[1] * trajectory.y- n_chap[2] * trajectory.z))
        # Alpha2 = np.exp(
        #        0. + 1j * self.omega * (trajectory.t + R/codata.c))
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

    def energy_radiated_farfield2(self, trajectory, gamma, x, y, distance):
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

        Alpha1 = np.exp(
            0. + 1j * self.omega * (trajectory.t - n_chap[0] * trajectory.x
                                    - n_chap[1] * trajectory.y - n_chap[2] * trajectory.z))
        Alpha2 = codata.c / (gamma ** 2 * R)
        Alpha3 = (1.0 / (1.0 - n_chap[0] * trajectory.v_x
                         - n_chap[1] * trajectory.v_y - n_chap[2] * trajectory.v_z)) ** 2
        integrand[0] += (n_chap[1] * A3 - n_chap[2] * A2
                         + Alpha2 * (n_chap[0] - trajectory.v_x)*Alpha3
                         ) * Alpha1
        integrand[1] += (-n_chap[0] * A3 + n_chap[2] * A1
                         + Alpha2 * (n_chap[1] - trajectory.v_y)*Alpha3
                         ) * Alpha1
        integrand[2] += (n_chap[0] * A2 - n_chap[1] * A1
                         + Alpha2 * (n_chap[2] - trajectory.v_z)*Alpha3
                         ) * Alpha1
        for k in range(3):
            # E[k] = np.trapz(integrand[k], self.trajectory.t)
            E[k] = integrate.simps(integrand[k], trajectory.t)
        E *= self.omega * 1j

        Alpha3 = (np.exp(1j * self.omega * (trajectory.t[-1] - n_chap[2] * trajectory.z[-1])) -
                  np.exp(1j * self.omega * (trajectory.t[0] - n_chap[2] * trajectory.z[0])))
        terme_bord = np.full((3,), 0. + 1j * 0., dtype=np.complex)
        terme_bord[0] += n_chap[0] * n_chap[2] * trajectory.v_z[0]
        terme_bord[1] += n_chap[1] * n_chap[2] * trajectory.v_z[0]
        terme_bord[2] += trajectory.v_z[0] * (n_chap[2] ** 2 - 1.0)
        terme_bord *= Alpha3 / (1.0 - n_chap[2] * trajectory.v_z[0])
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
            E[k] = integrate.simps(integrand[k], trajectory.t)
        return (np.abs(E[0]) ** 2 + np.abs(E[1]) ** 2 + np.abs(E[2]) ** 2)

    # attention ici c'est faux
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
        integrand[0] += (n_chap[1] * A3 - n_chap[2] * A2
                         + Alpha2 * (n_chap[0] - trajectory.v_x)
                         ) * Alpha1
        integrand[1] += (-n_chap[0] * A3 + n_chap[2] * A1
                         + Alpha2 * (n_chap[1] - trajectory.v_y)
                         ) * Alpha1
        integrand[2] += (n_chap[0] * A2 - n_chap[1] * A1
                         + Alpha2 * (n_chap[2] - trajectory.v_z)
                         ) * Alpha1
        for k in range(3):
            # E[k] = np.trapz(integrand[k], self.trajectory.t)
            E[k] = integrate.simps(integrand[k], trajectory.t)
        E *= self.omega * 1j

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


if __name__ == "__main__" :
    und_test = Undulator(K=1.87, E=1.3e9, lambda_u=0.035, L=0.035 * 12, I=1.0)
    traj_test=TrajectoryFactory(Nb_pts=1001,method=TRAJECTORY_METHOD_ANALYTIC).create_for_plane_undulator_ideal(
                                                                                           undulator=und_test)

    # start_time=time.time()
    rad=RadiationFactory(omega=und_test.omega1(),method=RADIATION_METHOD_APPROX_FARFIELD,Nb_pts=101
                         ).create_for_single_electron(trajectory=traj_test, undulator=und_test)
    print(rad.X)

    # interval = time.time() - start_time
    # print("interval temps :")
    # print(interval)
    # distance=100
    # X = np.linspace(0.0, distance*1.01e-3, 101)
    # Y = np.linspace(0.0,distance*1.01e-3, 101)
    # start_time=time.time()
    # rad=RadiationFactory(omega=5.0*und_test.omega1(),method=RADIATION_METHOD_APPROX_FARFIELD).create_for_single_electron(
    #                     trajectory=traj_test, undulator=und_test)
    # interval = time.time() - start_time
    # print("interval temps :")
    # print(interval)
    # print(rad.max())
    rad.draw()

