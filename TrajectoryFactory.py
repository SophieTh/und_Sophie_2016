import numpy as np
from pylab import *
import scipy.constants as codata
import scipy.integrate as integrate
from scipy.integrate import odeint
from scipy.interpolate import interp1d
from Trajectory import Trajectory

TRAJECTORY_METHOD_ANALYTIC=0
TRAJECTORY_METHOD_ODE=1
TRAJECTORY_METHOD_INTEGRATION=2


def fct_ODE_plane_undulator(y, t, cst, B):
    return [-cst * B(y[5]) * y[2],
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


class TrajectoryFactory(object):
    def __init__(self,K,E,lambda_u,Nb_period,Nb_pts,method):
        self.K=K
        self.E=E
        self.lambda_u=lambda_u
        self.Nb_period=Nb_period
        self.Nb_pts=Nb_pts
        self.method=method

    # calculate a theorical trajectory in an undulator
    def analytical_trajectory_undulator(self):
        N= self.Nb_pts
        ku = 2.0 * np.pi / self.lambda_u
        gamma = self.E / 0.511e6
        Beta = np.sqrt(1.0 - (1.0 / gamma ** 2))
        Beta_et = 1.0 - (1.0 / (2.0 * gamma ** 2)) * (1.0 + (self.K ** 2) / 2.0)
        omega_u = Beta_et * codata.c * ku
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
        trajectory[0] = np.linspace(-(self.lambda_u / (codata.c * Beta_et)) * (self.Nb_period / 2),
                                    (self.lambda_u / (codata.c * Beta_et)) * (self.Nb_period / 2), N)
        # X et Z en fct de t
        trajectory[3] = Beta_et * trajectory[0] - ((self.K / gamma) ** 2) * (1.0 / (8.0 * ku * codata.c)) * np.sin(
            2.0 * omega_u * trajectory[0])
        trajectory[1] = (self.K / (gamma * Beta_et * ku * codata.c)) * np.sin(omega_u * trajectory[0])
        # Vx et Vz en fct de t
        trajectory[6] = Beta_et - ((self.K / gamma) ** 2) * (1.0 / 4.0) * np.cos(2.0 * omega_u * trajectory[0])
        trajectory[4] = (self.K / (gamma * ku * Beta_et * codata.c)) * omega_u * np.cos(omega_u * trajectory[0])
        # Ax et Az en fct de t
        trajectory[9] = ((omega_u * self.K / gamma) ** 2) * (1.0 / (2.0 * ku * codata.c)) * np.sin(
            2.0 * omega_u * trajectory[0])
        trajectory[7] = -(self.K / (gamma * ku * Beta_et * codata.c)) * (omega_u ** 2) * np.sin(omega_u * trajectory[0])
        return trajectory

    # simulate the magnatic field in a undulator :
    def creation_magnetic_field_plane_undulator(self,z):
        Bo = self.K / (93.4 * self.lambda_u)
        By = -Bo * np.sin((2.0 * np.pi / self.lambda_u) * z)
        # Hamming windowing
        windpar = 1.0 / (2.0 * self.Nb_period)
        zmin = z.min()
        apo1 = zmin + windpar
        apo2 = z.max() - windpar
        wind = np.ones(len(z))
        for i in range(len(z)):
            if z[i] <= apo1:
                wind[i] *= 1.08 - (.54 + 0.46 * np.cos(np.pi * (z[i] - zmin) / windpar))
            if z[i] >= apo2:
                wind[i] *= 1.08 - (.54 - 0.46 * np.cos(np.pi * (z[i] - apo2) / windpar))
        By *= wind
        return By

    # electron's trajectory in a PLANE undulator that mean :  B=(0,By,0)
    # other hypothesis norm(v)=constant
    def trajectory_undulator_from_magnetic_field_method1(self,By, Z,Vx, Xo, Zo, nb_enlarg):
        gamma = self.E / 0.511e6
        Beta = np.sqrt(1.0 - (1.0 / gamma ** 2))
        Beta_et = 1.0 - (1.0 / (2.0 * gamma ** 2)) * (1.0 + (self.K ** 2) / 2.0)
        N=self.Nb_pts
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
        trajectory[0] = np.linspace(Z[0] / (Beta_et * codata.c), Z[len(Z) - 1] / (Beta_et * codata.c), N)

        Z, By = enlargement_vector_for_interpolation(Z, By, nb_enlarg)
        B = interp1d(Z, By)
        By2 = B(Beta_et * codata.c * trajectory[0])
        # Ax(t)
        Xm = codata.e * Beta_et / (gamma * codata.m_e)
        trajectory[7] = Xm * By2
        # Vx et Vz
        for i in range(N):
            trajectory[4][i] = integrate.simps(trajectory[7][0:(i + 1)], trajectory[0][0:(i + 1)]) + Vx / codata.c
        trajectory[6] = np.sqrt((Beta) ** 2 - trajectory[4] ** 2)
        # X et Z
        for i in range(N):
            trajectory[1][i] = integrate.simps(trajectory[4][0:(i + 1)], trajectory[0][0:(i + 1)]) + Xo / codata.c
            trajectory[3][i] = integrate.simps(trajectory[6][0:(i + 1)], trajectory[0][0:(i + 1)]) + Zo / codata.c
            # Az
        trajectory[9] = -(trajectory[7] * trajectory[4]) / trajectory[6]

        return trajectory

    # electron's trajectory in a PLANE undulator that mean :  B=(0,By,0)
    def trajectory_undulator_from_magnetic_field_method2(self,By, Z, nb_enlarg, Vx, Vz, Xo):
        gamma = self.E / 0.511e6
        Beta = np.sqrt(1.0 - (1.0 / gamma ** 2))
        Beta_et = 1.0 - (1.0 / (2.0 * gamma ** 2)) * (1.0 + (self.K ** 2) / 2.0)
        N=self.Nb_pts
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
        trajectory = np.zeros((10, N))
        trajectory[0] = np.linspace(Z[0] / (Beta_et * codata.c), Z[- 1] / (Beta_et * codata.c), N)
        boundary_condition = [Vx, 0.0, Vz, Xo, 0.0, Z[0]]
        Z, By = enlargement_vector_for_interpolation(Z, By, nb_enlarg * 100)
        B = interp1d(Z, By)
        cst = -codata.e / (codata.m_e * gamma)
        res = odeint(fct_ODE_plane_undulator, boundary_condition, trajectory[0], args=(cst, B), full_output=True)
        traj = res[0]
        info = res[1]
        # print("1 : nonstiff problems, Adams . 2: stiff problem, BDF")
        # print(info.get('mused'))
        traj = np.transpose(traj)
        trajectory[4] = traj[0]
        trajectory[5] = traj[1]
        trajectory[6] = traj[2]
        trajectory[1] = traj[3]
        trajectory[2] = traj[4]
        trajectory[3] = traj[5]
        trajectory[7] = -cst * B(trajectory[3]) * trajectory[6]
        trajectory[9] = cst * B(trajectory[3]) * trajectory[4]
        k = 1
        while k < 10:
            trajectory[k] *= 1.0 / codata.c
            k += 1
        return trajectory

    def create_for_plane_undulator(self,Vx=None,Vz=None,Xo=None,Z_By=None):
        if (self.method == 1 or self.method == 2):
            if Z_By == None:
                Z = np.linspace(-self.lambda_u * (self.Nb_period / 2), self.lambda_u * (self.Nb_period / 2),
                                self.Nb_pts)
                B_y = self.creation_magnetic_field_plane_undulator(Z)
            else:
                Z = Z_By[0]
                B_y = Z_By[1]

            if (Vx==None or Vz==None or Xo==None) :
                Vx=0.0
                Vz=np.sqrt(1.0 - (1.0 / ( self.E /0.511e6)** 2))*codata.c # = Beta*c=sqrt(1-1/gamma**2)*c
                Xo=0.0

            if self.method == 1:
                T= self.trajectory_undulator_from_magnetic_field_method1(By=B_y, Z=Z,
                                                                         nb_enlarg=np.floor(len(Z) / (0.1 * len(B_y))),
                                                                             Vx=Vx, Xo=Xo, Zo=Z[0])

            else:
                T= self.trajectory_undulator_from_magnetic_field_method2(By=B_y, Z=Z,
                                                                    nb_enlarg=np.floor(len(Z) / (0.1 * self.Nb_pts)),
                                                                       Vx=Vx, Vz=Vz, Xo=Xo)

        else:
            T = self.analytical_trajectory_undulator()

        trajectory = Trajectory(T[0],T[1],T[2],T[3],T[4],T[5],T[6],T[7],T[8],T[9],self)

        return trajectory

    def create_for_plane_undulator_ideal(self):
        T = self.analytical_trajectory_undulator()
        trajectory = Trajectory(T[0], T[1], T[2], T[3], T[4], T[5], T[6], T[7], T[8], T[9], self)
        return trajectory

