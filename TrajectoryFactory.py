import numpy as np
import scipy.constants as codata
import scipy.integrate as integrate
from scipy.integrate import odeint
from scipy.interpolate import interp1d
from Trajectory import Trajectory
from UndulatorParameter import UndulatorParameters as Undulator


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

def fct_ODE_undulator(y, t, cst, Bx,By,Bz):
    return [cst * (Bz(y[5]) * y[1] - By(y[5]) * y[2]),
            cst * (Bx(y[5]) * y[2] - Bz(y[5]) * y[0]),
            cst * (By(y[5]) * y[0] - Bx(y[5]) * y[1]),
            y[0],
            y[1],
            y[2]]



'''
initial condition : [Vx,Vy,Vz,x,y,z]
'''
class TrajectoryFactory(object):
    def __init__(self,Nb_pts,method,initial_condition=None):
        self.Nb_pts = Nb_pts
        self.method = method
        self.initial_condition=initial_condition

    def copy(self):
        if self.initial_condition==None :
            cond=None
        else :
            cond= self.initial_condition.copy()
        return TrajectoryFactory(Nb_pts=self.Nb_pts,method=self.method,
                                 initial_condition=cond)

    # calculate a theorical trajectory in an undulator
    def analytical_trajectory_plane_undulator(self,undulator):
        N= self.Nb_pts
        ku = 2.0 * np.pi / undulator.lambda_u
        gamma = undulator.E / 0.511e6
        Beta = np.sqrt(1.0 - (1.0 / gamma ** 2))
        Beta_et = 1.0 - (1.0 / (2.0 * gamma ** 2)) * (1.0 + (undulator.K ** 2) / 2.0)
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
        trajectory[0] = np.linspace(-undulator.L / (2.0 * codata.c * Beta_et ),
                                     undulator.L / (2.0 * codata.c * Beta_et ), N)
        # X et Z en fct de t
        trajectory[3] = Beta_et * trajectory[0] - ((undulator.K / gamma) ** 2) * (1.0 / (8.0 * ku * codata.c)) * np.sin(
            2.0 * omega_u * trajectory[0])
        trajectory[1] = (undulator.K / (gamma * Beta_et * ku * codata.c)) * np.sin(omega_u * trajectory[0])
        # Vx et Vz en fct de t
        trajectory[6] = Beta_et - ((undulator.K / gamma) ** 2) * (1.0 / 4.0) * np.cos(2.0 * omega_u * trajectory[0])
        trajectory[4] = (undulator.K / (gamma * ku * Beta_et * codata.c)) * omega_u * np.cos(omega_u * trajectory[0])
        # Ax et Az en fct de t
        trajectory[9] = ((omega_u * undulator.K / gamma) ** 2) * (1.0 / (2.0 * ku * codata.c)) * np.sin(
            2.0 * omega_u * trajectory[0])
        trajectory[7] = -(undulator.K / (gamma * ku * Beta_et * codata.c)) * (omega_u ** 2) * np.sin(omega_u * trajectory[0])
        return trajectory


    #method aui ne marche pas
    # electron's trajectory in a PLANE undulator that mean :  B=(0,By,0)
    # other hypothesis norm(v)=constant
    def trajectory_undulator_from_magnetic_field_method1(self,undulator, B):
        gamma = undulator.E / 0.511e6
        Beta = np.sqrt(1.0 - (1.0 / gamma ** 2))
        Beta_et = 1.0 - (1.0 / (2.0 * gamma ** 2)) * (1.0 + (undulator.K ** 2) / 2.0)
        N=self.Nb_pts
        Z=B.z
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
        trajectory[0] = np.linspace(Z[0] / (Beta_et * codata.c), Z[- 1] / (Beta_et * codata.c), N)

        By2 = B.By(Z)
        # Ax(t)
        Xm = codata.e * Beta_et / (gamma * codata.m_e)
        trajectory[7] = Xm * By2
        # Vx et Vz
        for i in range(N):
            trajectory[4][i] = integrate.simps(trajectory[7][0:(i + 1)], trajectory[0][0:(i + 1)]) \
            + self.initial_condition[0]/ codata.c
        trajectory[6] = np.sqrt((Beta) ** 2 - trajectory[4] ** 2)
        # X et Z
        for i in range(N):
            trajectory[1][i] = integrate.simps(trajectory[4][0:(i + 1)], trajectory[0][0:(i + 1)]) \
                               + self.initial_condition[3]/ codata.c
            trajectory[3][i] = integrate.simps(trajectory[6][0:(i + 1)], trajectory[0][0:(i + 1)]) \
                               +  self.initial_condition[5]/ codata.c

            # Az
        trajectory[9] = -(trajectory[7] * trajectory[4]) / trajectory[6]

        return trajectory

    # electron's trajectory in a PLANE undulator that mean :  B=(0,By,0)
    def trajectory_undulator_from_magnetic_field_method2(self,undulator, B):
        gamma = undulator.gamma()
        Z = B.z
        Beta = np.sqrt(1.0 - (1.0 / gamma ** 2))
        Beta_et = 1.0 - (1.0 / (2.0 * gamma ** 2)) * (1.0 + (undulator.K ** 2) / 2.0)
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
        cst = -codata.e / (codata.m_e * gamma)
        # res = odeint(fct_ODE_undulator, self.initial_condition, trajectory[0],
        #              args=(cst, B.Bx,B.By,B.Bz), full_output=True)
        res = odeint(fct_ODE_plane_undulator,self.initial_condition, trajectory[0],
                     args=(cst,B.By), full_output=True)
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
        trajectory[7] = -cst * B.By(trajectory[3]) * trajectory[6]
        trajectory[9] = cst * B.By(trajectory[3]) * trajectory[4]
        k = 1
        while k < 10:
            trajectory[k] *= 1.0 / codata.c
            k += 1
        return trajectory


    def calculate_trajectory(self,undulator,B):
        if (self.method == TRAJECTORY_METHOD_ODE or self.method == TRAJECTORY_METHOD_INTEGRATION):
            if self.method == TRAJECTORY_METHOD_INTEGRATION:
                T = self.trajectory_undulator_from_magnetic_field_method1(undulator=undulator, B=B)


            else: # method=TRAJECTORY_METHOD_ODE
                T = self.trajectory_undulator_from_magnetic_field_method2(undulator=undulator, B=B)

        else:
            T = self.analytical_trajectory_plane_undulator(undulator=undulator)
        return T


    def create_for_plane_undulator(self,undulator,B):
        if (self.method == 1 or self.method == 2):
            if (self.initial_condition==None) :
                self.initial_condition=np.array([0.0,0.0,np.sqrt(1.0 - (1.0 / ( undulator.E /0.511e6)** 2))*codata.c,
                                                 0.0,0.0,B.z[0]])
            T=self.calculate_trajectory(undulator=undulator,B=B)

        else:
            T = self.analytical_trajectory_plane_undulator(undulator=undulator)

        trajectory = Trajectory(T[0],T[1],T[2],T[3],T[4],T[5],T[6],T[7],T[8],T[9])

        return trajectory


    def create_for_plane_undulator_ideal(self,undulator):
        T = self.analytical_trajectory_plane_undulator(undulator)
        trajectory = Trajectory(T[0], T[1], T[2], T[3], T[4], T[5], T[6], T[7], T[8], T[9])
        return trajectory






if __name__ == "__main__" :
    und_test = Undulator(K=1.87, E=1.3e9, lambda_u=0.035, L=0.035 * 12, I=1.0)
    traj_test=TrajectoryFactory(Nb_pts=201,method=TRAJECTORY_METHOD_ANALYTIC).create_for_plane_undulator_ideal(
                                                                                                undulator=und_test)
    traj_test.draw()