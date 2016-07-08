import numpy as np
import scipy.constants as codata
import scipy.integrate as integrate
from scipy.integrate import odeint
from pySRU.Trajectory import Trajectory
from pySRU.TrajectoryAnalitic import TrajectoryAnalitic
from pySRU.TrajectoryArray import TrajectoryArray
from pySRU.Parameter import Parameter,PLANE_UNDULATOR,BENDING_MAGNET
from pySRU.ParameterPlaneUndulator import ParameterPlaneUndulator as Undulator


TRAJECTORY_METHOD_ANALYTIC=0
TRAJECTORY_METHOD_ODE=1
TRAJECTORY_METHOD_INTEGRATION=2


def fct_ODE_plane_undulator(y, t, cst, B):
    return [-cst * B(y[5],y[4]) * y[2],
            0.0,
            cst * B(y[5],y[4]) * y[0],
            y[0],
            0.0,
            y[2]]

def fct_ODE_undulator(y, t, cst, Bx,By,Bz):

    return [cst * (Bz(z=y[5],y=y[4],x=y[3]) * y[1] - By(z=y[5],y=y[4],x=y[3]) * y[2]),
            cst * (Bx(z=y[5],y=y[4],x=y[3]) * y[2] - Bz(z=y[5],y=y[4],x=y[3]) * y[0]),
            cst * (By(z=y[5],y=y[4],x=y[3]) * y[0] - Bx(z=y[5],y=y[4],x=y[3]) * y[1]),
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
        #
        # t
        time = np.linspace(-undulator.L / (2.0 * codata.c * Beta_et),
                                     undulator.L / (2.0 * codata.c * Beta_et), N)
        #utiliser omegat!
        # # X et Z en fct de t
        z = (lambda t :Beta_et * t + ((undulator.K / gamma) ** 2) * (1.0 / (8.0 * omega_u)
                                    ) * np.sin( 2.0 * omega_u*t))
        x =(lambda t :(-(undulator.K / (gamma * omega_u)) * np.cos(omega_u*t)))
        # # Vx et Vz en fct de t
        v_z = (lambda  t :Beta_et + ((undulator.K / gamma) ** 2) * (1.0 / 4.0) * np.cos(2.0 *omega_u*t))
        v_x=(lambda t : (undulator.K / (gamma )) * np.sin(omega_u*t))
        # # Ax et Az en fct de t
        a_z=(lambda t :-omega_u *( undulator.K / gamma) ** 2 * 0.5 * np.sin( 2.0 * omega_u*t))
        # #trajectory[7] = (undulator.K / (gamma * ku * Beta_et * codata.c)) * (omega_u ** 2) * np.cos(omegat)
        a_x=(lambda  t : (undulator.K / (gamma )) * (omega_u ) * np.cos(omega_u*t))
        fct_null=(lambda t : 0.0*t)
        return TrajectoryAnalitic(t=time,x=x,y=fct_null,z=z,v_x=v_x,v_y=fct_null,v_z=v_z,
                                  a_x=a_x,a_y=fct_null,a_z=a_z)

    # # a change !!
    def analytical_trajectory_cst_magnf(self,bending_magnet):
        #time= np.linspace(-undulator.L / (2.0 * codata.c * Beta_et),
                                     #undulator.L / (2.0 * codata.c * Beta_et), N)
        time = np.linspace(bending_magnet.Zo_analitic(),-bending_magnet.Zo_analitic(),
                           self.Nb_pts) / ( bending_magnet.Beta_et()* codata.c)
        f = -bending_magnet.Bo * codata.e / (bending_magnet.gamma() * codata.m_e)
        vz_0 = bending_magnet.Beta()
        x_0 = -bending_magnet.Beta() / f
        x = (lambda t:-vz_0 * (1.0 / f) * np.cos(f * t)+vz_0/f + x_0)
        y = (lambda t:vz_0 * 0.0 * t)
        z = (lambda t:vz_0 * (1.0 / f) * np.sin(f * t))
        vx = (lambda t:vz_0 * np.sin(f * t))
        vy = (lambda t:vz_0 * 0.0 * t)
        vz = (lambda t:vz_0 * np.cos(f * t))
        ax =(lambda t: vz_0 * f * np.cos(f * t))
        ay = (lambda t:vz_0 * 0.0 * t)
        az = (lambda t:-vz_0* f * np.sin(f * t))
        return TrajectoryAnalitic(t=time,x=x,y=y,z=z,v_x=vx,v_y=vy,v_z=vz,a_x=ax,a_y=ay,a_z=az)
    # # a change !!


    # ### cree speciale ment pour un test
    # def analytical_trajectory_cst_magnf(self,Bo,gamma,t,vz_0,x_0):
    #     T=self.copy()
    #     f=Bo*codata.e/(gamma*codata.m_e)
    #     T.x = (lambda t:-vz_0 * (1.0 / f) * np.cos(f * t) +vz_0/f + x_0)
    #     T.y = (lambda t: vz_0 * (0.0) * t)
    #     T.z = (lambda t: vz_0 * (1.0 / f) * np.sin(f * t))
    #     T.v_x = (lambda t: vz_0 * np.sin(f * t))
    #     T.v_y = (lambda t: vz_0 * 0.0 * t)
    #     T.v_z = (lambda t: vz_0 * np.cos(f * t))
    #     T.a_x = (lambda t:vz_0 * f * np.cos(f * t))
    #     T.a_y = (lambda t:vz_0 * 0.0 * t)
    #     T.a_z = (lambda t:-vz_0* f * np.sin(f * t))
    #     return T


    # electron's trajectory in a PLANE undulator that mean :  B=(0,By,0)
    # other hypothesis norm(v)=constant
    def trajectory_from_magnetic_field_method_INT(self, parameter, B):
        gamma=parameter.gamma()
        Beta = parameter.Beta()
        Beta_et = parameter.Beta_et()
        N=self.Nb_pts
        Z=np.linspace(B.z[0],B.z[-1],N)
        if type(B.y)== np.ndarray :
            Y = np.linspace(B.y[0], B.y[-1], N)
        else :
            Y=B.y
        if type(B.x) == np.ndarray:
            X = np.linspace(B.x[0], B.x[-1], N)
        else:
            X = B.x
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
        trajectory[0] = Z/(Beta_et*codata.c)
        # Ax
        Xm = codata.e *Beta_et/ (gamma * codata.m_e )
        trajectory[7] = Xm * B.By(Z,Y,X)
        # Vx et Vz
        for i in range(N):
            #np.trapz
            # trajectory[4][i] = integrate.simps(trajectory[7][0:(i + 1)], trajectory[0][0:(i + 1)]) \
            trajectory[4][i] = integrate.simps(trajectory[7][0:(i + 1)], trajectory[0][0:(i + 1)]) \
                                       + self.initial_condition[0]/ codata.c
        trajectory[6] = np.sqrt((Beta) ** 2 - trajectory[4] ** 2)
        # X et Z
        for i in range(N):
            #trajectory[1][i] = integrate.simps(trajectory[4][0:(i + 1)], trajectory[0][0:(i + 1)]) \
            trajectory[1][i] =  np.trapz(trajectory[4][0:(i + 1)], trajectory[0][0:(i + 1)]) \
                               + self.initial_condition[3]/ codata.c
            #trajectory[3][i] = integrate.simps(trajectory[6][0:(i + 1)], trajectory[0][0:(i + 1)]) \
            trajectory[3][i] = np.trapz(trajectory[6][0:(i + 1)], trajectory[0][0:(i + 1)]) \
                               +  self.initial_condition[5]/ codata.c

            # Az
        trajectory[9] = -(trajectory[7] * trajectory[4]) / trajectory[6]

        return trajectory

    # electron's trajectory in a PLANE undulator that mean :  B=(0,By,0)
    def trajectory_from_magnetic_field_method_ODE(self, parameter, B):
        gamma = parameter.gamma()
        Beta_et=parameter.Beta_et()
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
        trajectory = np.zeros((10,self.Nb_pts))
        trajectory[0] = np.linspace(B.z[0] / (Beta_et * codata.c),
                                    B.z[-1] / (Beta_et * codata.c), self.Nb_pts)
        #trajectory[0] = np.linspace(Z[0] / (Beta_et * codata.c), Z[- 1] / (Beta_et * codata.c), N)
        cst = -codata.e / (codata.m_e * gamma)
        res = odeint(fct_ODE_undulator,self.initial_condition, trajectory[0],
                     args=(cst,B.Bx,B.By,B.Bz), full_output=True)
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
        trajectory[7] = - cst * B.By(trajectory[3], trajectory[2],trajectory[1]) * trajectory[6]
        trajectory[9] = cst * B.By(trajectory[3], trajectory[2],trajectory[1]) * trajectory[4]
        k = 1
        while k < 10:
            trajectory[k] *= 1.0 / codata.c
            k += 1
        return trajectory

    def calculate_trajectory(self, parameter, B):
        if (self.method == TRAJECTORY_METHOD_ODE or self.method == TRAJECTORY_METHOD_INTEGRATION):
            if self.method == TRAJECTORY_METHOD_INTEGRATION:
                T = self.trajectory_from_magnetic_field_method_INT(parameter=parameter, B=B)


            else: # method=TRAJECTORY_METHOD_ODE
                T = self.trajectory_from_magnetic_field_method_ODE(parameter=parameter, B=B)

        else:
            if parameter.type_magnet==PLANE_UNDULATOR:
                T = self.analytical_trajectory_plane_undulator(undulator=parameter)
            else :
                T = self.analytical_trajectory_cst_magnf(undulator=parameter,B=B)
        return T

    def create_for_parameter(self,parameter,B):
        if (self.method == TRAJECTORY_METHOD_INTEGRATION or self.method == TRAJECTORY_METHOD_ODE):
            if (self.initial_condition==None) :
                self.initial_condition=np.array([0.0,0.0,np.sqrt(1.0 - (1.0 / ( parameter.E /0.511e6)** 2))*codata.c,
                                                 0.0,0.0,B.z[0]])
                #print(self.initial_condition)
            T=self.calculate_trajectory(parameter=parameter,B=B)
            trajectory = TrajectoryArray(T[0], T[1], T[2], T[3], T[4], T[5], T[6], T[7], T[8], T[9])
        else:
            if parameter.type_magnet==PLANE_UNDULATOR:
                trajectory = self.analytical_trajectory_plane_undulator(undulator=parameter)
            else :
                trajectory = self.analytical_trajectory_cst_magnf(bending_magnet=parameter)
            t0=trajectory.t[0]
            self.initial_condition = np.array([trajectory.v_x(t0), trajectory.v_y(t0),
                                    trajectory.v_z(t0), trajectory.x(t0),
                                    trajectory.y(t0), trajectory.z(t0)])
            self.initial_condition *= codata.c
        # a changer TODO
        if (type(trajectory)==TrajectoryAnalitic) :
            trajectory=trajectory.convert()

        return trajectory

    def create_for_plane_undulator_ideal(self,undulator):
        trajectory = self.analytical_trajectory_plane_undulator(undulator=undulator)
        self.initial_condition = np.array([trajectory.v_x(trajectory.t[0]), trajectory.v_y(trajectory.t[0]),
                                           trajectory.v_z(trajectory.t[0]), trajectory.x(trajectory.t[0]),
                                           trajectory.y(trajectory.t[0]), trajectory.z(trajectory.t[0])])
        self.method=TRAJECTORY_METHOD_ANALYTIC
        return trajectory

    def create_for_cst_magnetic_field_ideal(self, parameter):
        T = self.analytical_trajectory_cst_magnf(bending_magnet=parameter)
        return T

    def create_from_array(self,array):
        if array.shape[0] != 10 :
            raise Exception('this array can not be convert in Trajectory')
        return TrajectoryArray(t=array[0],x=array[1],y=array[2],z=array[3],v_x=array[4],v_y=array[5],
                    v_z = array[6],a_x=array[7],a_y=array[8],a_z=array[9])




if __name__ == "__main__" :
    und_test = Undulator(K=1.87, E=1.3e9, lambda_u=0.035, L=0.035 * 12, I=1.0)
    traj_test=TrajectoryFactory(Nb_pts=201,method=TRAJECTORY_METHOD_ANALYTIC).create_for_plane_undulator_ideal(
                                                                                                undulator=und_test)
    Beta_et = 1.0 - (1.0 / (2.0 * und_test.gamma() ** 2)) * (1.0 + (und_test.K ** 2) / 2.0)

    ku=2.0*np.pi/und_test.lambda_u

    xo = und_test.K / (und_test.gamma() * Beta_et * ku)
    zo=Beta_et*codata.c*und_test.L / (2.0 * codata.c * Beta_et)
    vzo=Beta_et*codata.c - ((und_test.K / und_test.gamma()) ** 2)

    initial_condition=np.array([0.0,0.0,vzo,xo,0.0,zo])
    traj_test2 = TrajectoryFactory(Nb_pts=201, method=TRAJECTORY_METHOD_INTEGRATION,
                                   initial_condition=initial_condition).create_for_plane_undulator_ideal(
                                    undulator=und_test)

    traj_test.plot_2_trajectory(traj_test2)
    print(traj_test.z.max())