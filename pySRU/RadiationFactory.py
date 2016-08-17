import numpy as np
from pySRU.Radiation import Radiation
import scipy.constants as codata
import scipy.integrate as integrate
from pySRU.TrajectoryFactory import TrajectoryFactory, TRAJECTORY_METHOD_ANALYTIC , TRAJECTORY_METHOD_ODE




RADIATION_METHOD_NEAR_FIELD=0
RADIATION_METHOD_APPROX=1
RADIATION_METHOD_APPROX_FARFIELD=2
eV_to_J=1.602176487e-19


class RadiationFactory(object):
    def __init__(self,method,omega,Nb_pts=101):
        self.omega=omega
        self.method=method
        #TODO useful ?
        self.Nb_pts=Nb_pts


    def copy(self):
        return RadiationFactory( method=self.method,omega=self.omega,Nb_pts=self.Nb_pts)

    def energy_eV(self):
        return self.omega*codata.hbar/eV_to_J

    # Photon's flow all over a screen situate at distance D of an undulator
    def create_for_one_relativistic_electron(self, trajectory, source, XY_are_list=False, distance=None, X=None, Y=None):
        if X is None or Y is None:
            print('calculate X and Y array')
            if distance == None:
                if self.method==RADIATION_METHOD_APPROX_FARFIELD :
                    xy_max = np.tan(1./ source.Lorentz_factor())
                else :
                    distance=source.choose_distance_automatic(2)
                    xy_max = np.tan(1. / source.Lorentz_factor())*distance
            else :
                xy_max = np.tan(1. / source.Lorentz_factor())*distance
            X = np.linspace(0.0, xy_max, self.Nb_pts)
            Y = np.linspace(0.0, xy_max, self.Nb_pts)
        else :
            if type(X) != np.ndarray:

                X = np.linspace(0.0, X, self.Nb_pts)
            if type(Y) != np.ndarray:
                Y = np.linspace(0.0, Y, self.Nb_pts)
            #Nb_pts=len(X)

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
        # c4 = 1e-3 / omega1  # to get .1% energy bandwidth  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        # c5 = 1e-3 * 1e-3  # from rad to .1mrad angle bandwidth



        c6 = codata.e * source.I_current() * 1e-9 / (8.0 * np.pi ** 2 * codata.epsilon_0 * codata.c * codata.h)

        # to guarantee the flux in phot/mm2 reduces with 1/r2
        if distance != None:
            c6 /= distance**2

        if X_array.size != Y_array.size:
            raise Exception("X and Y dimensions must be equal.")

        res = np.zeros_like(X_array)
        shape1 = res.shape
        X = X_array.flatten()
        Y = Y_array.flatten()
        res = res.flatten()
        gamma = source.Lorentz_factor()
        if self.method == RADIATION_METHOD_NEAR_FIELD:
            for i in range(len(X)):
                res[i] = c6 * self.energy_radiated_near_field(trajectory=trajectory,gamma=gamma,
                                                                                  distance=distance, x=X[i], y=Y[i])
        elif self.method == RADIATION_METHOD_APPROX:
            for i in range(len(X)):
                res[i] = c6 * self.energy_radiated_approx(trajectory=trajectory,
                                                                                  distance=distance, x=X[i], y=Y[i])
        else :
            for i in range(len(X)):
                res[i] = c6 * self.energy_radiated_approximation_and_farfield(trajectory=trajectory,
                                                            distance=distance,x=X[i], y=Y[i])
        res = res.reshape(shape1)
        return res

    # on neglige le deuxieme terme et on considere n comme constant
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



    # on neglige le deuxieme terme
    def energy_radiated_approx(self,trajectory, x, y,distance):
        N = trajectory.nb_points()
        n_chap = np.array([x - trajectory.x * codata.c, y - trajectory.y * codata.c, distance - trajectory.z * codata.c])
        # R = np.zeros(n_chap.shape[1])
        # for i in range(n_chap.shape[1]):
        #     R[i] = np.linalg.norm(n_chap[:, i])
        #     n_chap[:, i] /= R[i]
        R = np.sqrt( n_chap[0]**2 + n_chap[1]**2 + n_chap[2]**2 )
        n_chap[0,:] /= R
        n_chap[1,:] /= R
        n_chap[2,:] /= R

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


    # formule general
    def energy_radiated_near_field(self, gamma, trajectory, x, y, distance):
        N = trajectory.nb_points()
        n_chap = np.array(
            [x - trajectory.x * codata.c, y - trajectory.y * codata.c, distance - trajectory.z * codata.c])
        # R = np.zeros(n_chap.shape[1])
        # for i in range(n_chap.shape[1]):
        #     R[i] = np.linalg.norm(n_chap[:, i])
        #     n_chap[:, i] /= R[i]

        R = np.sqrt( n_chap[0]**2 + n_chap[1]**2 + n_chap[2]**2 )
        n_chap[0,:] /= R
        n_chap[1,:] /= R
        n_chap[2,:] /= R

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




    def get_method(self):
        if self.method == RADIATION_METHOD_NEAR_FIELD:
            method = 'Near field calculation'

        elif self.method == RADIATION_METHOD_APPROX_FARFIELD:
            method = '  Far field and approximation calculation '
        else :
            method=' Problem unknow method'
        return method

    def print_parameters(self):
        print('Radiation ')
        print('    method: %s' %self.get_method())
        print('    number of points in each direction  %d' %self.Nb_pts)
        print('    energy of the emission: %f eV, omega: %f Hz' %(self.energy_eV(),self.omega))


def Exemple_FARFIELD():
    from pySRU.MagneticStructureUndulatorPlane import MagneticStructureUndulatorPlane as Undulator
    from pySRU.ElectronBeam import ElectronBeam
    from pySRU.SourceUndulatorPlane import SourceUndulatorPlane
    from pySRU.TrajectoryFactory import TrajectoryFactory,TRAJECTORY_METHOD_ODE,TRAJECTORY_METHOD_ANALYTIC
    import time

    undulator_test = Undulator(K=1.87, period_length=0.035, length=0.035 * 14)
    electron_beam_test = ElectronBeam(Electron_energy=1.3, I_current=1.0)
    magnetic_field_test = undulator_test.create_magnetic_field()

    magnetic_field_test.plot_z(0,0,np.linspace(-0.035 * (14+8) / 2,0.035 * (14+8) / 2,500))

    source_test = SourceUndulatorPlane(undulator=undulator_test,
                         electron_beam=electron_beam_test,
                         magnetic_field=magnetic_field_test)

    traj = TrajectoryFactory(Nb_pts=2000, method=TRAJECTORY_METHOD_ODE).create_from_source(source_test)

    t0 = time.time()
    Rad = RadiationFactory(omega=source_test.harmonic_frequency(1), method=RADIATION_METHOD_APPROX_FARFIELD, Nb_pts=101
                            ).create_for_one_relativistic_electron(trajectory=traj, source=source_test)
    print("Elapsed time in RadiationFactory: ",time.time()-t0)

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
    Rad.plot(title="FAR FIELD")


def Exemple_NEARFIELD():
    from pySRU.MagneticStructureUndulatorPlane import MagneticStructureUndulatorPlane as Undulator
    from pySRU.ElectronBeam import ElectronBeam
    from pySRU.SourceUndulatorPlane import SourceUndulatorPlane
    import time

    undulator_test = Undulator(K=1.87, period_length=0.035, length=0.035 * 14)
    electron_beam_test = ElectronBeam(Electron_energy=1.3, I_current=1.0)
    magnetic_field_test = undulator_test.create_magnetic_field()
    source_test = SourceUndulatorPlane(undulator=undulator_test,
                         electron_beam=electron_beam_test,
                         magnetic_field=magnetic_field_test)

    traj = TrajectoryFactory(Nb_pts=2000, method=TRAJECTORY_METHOD_ODE).create_from_source(source_test)

    t0 = time.time()
    Rad = RadiationFactory(omega=source_test.harmonic_frequency(1), method=RADIATION_METHOD_NEAR_FIELD, Nb_pts=101
                           ).create_for_one_relativistic_electron(trajectory=traj, source=source_test,distance=100)
    print("Elapsed time in RadiationFactory: ",time.time()-t0)

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

def Exemple_APPROX():
    from pySRU.MagneticStructureUndulatorPlane import MagneticStructureUndulatorPlane as Undulator
    from pySRU.ElectronBeam import ElectronBeam
    from pySRU.SourceUndulatorPlane import SourceUndulatorPlane
    import time

    undulator_test = Undulator(K=1.87, period_length=0.035, length=0.035 * 14)
    electron_beam_test = ElectronBeam(Electron_energy=1.3, I_current=1.0)
    magnetic_field_test = undulator_test.create_magnetic_field()
    source_test = SourceUndulatorPlane(undulator=undulator_test,
                         electron_beam=electron_beam_test,
                         magnetic_field=magnetic_field_test)

    traj = TrajectoryFactory(Nb_pts=2000, method=TRAJECTORY_METHOD_ODE).create_from_source(source_test)

    t0 = time.time()
    Rad = RadiationFactory(omega=source_test.harmonic_frequency(1), method=RADIATION_METHOD_APPROX, Nb_pts=101
                           ).create_for_one_relativistic_electron(trajectory=traj, source=source_test,distance=100)
    print("Elapsed time in RadiationFactory: ",time.time()-t0)

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
    # Exemple_FARFIELD()
    # Exemple_APPROX()
    Exemple_NEARFIELD()



