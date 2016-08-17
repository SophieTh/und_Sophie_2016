import numpy as np
import scipy.constants as codata
from scipy.integrate import odeint
from pySRU.Trajectory import Trajectory
from pySRU.SourceUndulatorPlane import SourceUndulatorPlane,PLANE_UNDULATOR,BENDING_MAGNET
from pySRU.MagneticStructureUndulatorPlane import MagneticStructureUndulatorPlane as Undulator
from pySRU.ElectronBeam import ElectronBeam


TRAJECTORY_METHOD_ANALYTIC=0
TRAJECTORY_METHOD_ODE=1






def fct_ODE_magnetic_field(y, t, cst, Bx,By,Bz):
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
        if self.initial_condition is None :
            cond=None
        else :
            cond= self.initial_condition.copy()
        return TrajectoryFactory(Nb_pts=self.Nb_pts,method=self.method,
                                 initial_condition=cond)

    def choise_initial_condition(self,source):
        if self.method != TRAJECTORY_METHOD_ANALYTIC:
            self.initial_condition = source.choose_initial_contidion_automatic()

    # calculate a theorical trajectory in an undulator
    def analytical_trajectory_plane_undulator(self,undulator):
        ku = 2.0 * np.pi / undulator.magnetic_structure.period_length
        gamma = undulator.Lorentz_factor()
        Beta_et = undulator.average_z_speed_in_undulator()
        K=undulator.magnetic_structure.K
        omega_u = Beta_et * codata.c * ku

        t = undulator.analytical_times_vector(Nb_pts=self.Nb_pts)
        ## x and z
        z = Beta_et * t + ((K / gamma) ** 2) * (1.0 / (8.0 * omega_u)) * np.sin( 2.0 * omega_u*t)
        x = (-(K / (gamma * omega_u)) * np.cos(omega_u*t))
        # # Vx and Vz
        v_z = Beta_et + ((K / gamma) ** 2) * (1.0 / 4.0) * np.cos(2.0 *omega_u*t)
        v_x= (K / (gamma )) * np.sin(omega_u*t)
        # # Ax and Az
        a_z=-omega_u *(K / gamma) ** 2 * 0.5 * np.sin( 2.0 * omega_u*t)
        a_x= (K / (gamma )) * (omega_u ) * np.cos(omega_u*t)
        # y
        y=0.0*t
        v_y=y
        a_y=y
        return Trajectory(t=t,x=x,y=y,z=z,v_x=v_x,v_y=v_y,v_z=v_z,
                                  a_x=a_x,a_y=a_y,a_z=a_z)



    # #TODO changer pour mettre avec des conditions initial ?
    def analytical_trajectory_cst_magnf(self, bending_magnet):
        Bo= bending_magnet.magnetic_structure.Bo
        #ro=bending_magnet.magnetic_structure.radius_curvature(bending_magnet.E())
        omega_p=Bo*codata.e/(codata.m_e*bending_magnet.Lorentz_factor())
        vz_0 = bending_magnet.electron_speed()

        t=bending_magnet.analytical_times_vector(Nb_pts=self.Nb_pts)
        Zo = bending_magnet.magnetic_structure.L / 2.
        #t=np.linspace(0.0,2.*Zo/(vz_0*codata.c))
        to=t[0]
        x = (vz_0/omega_p)*(1.-np.cos(omega_p*(t-to)))
        y = 0.0 * (t)
        z = (vz_0/omega_p)*np.sin(omega_p*(t-to))-Zo/codata.c
        vx = vz_0 * np.sin(omega_p*(t-to))
        vy = vz_0 * 0.0 * (t)
        vz = vz_0 * np.cos(omega_p*(t-to))
        ax = vz_0 * omega_p * np.cos(omega_p*(t-to))
        ay = vz_0 * 0.0 * (t)
        az = -vz_0* omega_p * np.sin(omega_p*(t-to))

        return Trajectory(t=t,x=x,y=y,z=z,v_x=vx,v_y=vy,v_z=vz,a_x=ax,a_y=ay,a_z=az)




    def trajectory_from_magnetic_field_method_ODE(self, source):
        gamma = source.Lorentz_factor()

        B=source.magnetic_field
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
        trajectory[0] = source.construct_times_vector(initial_contition=self.initial_condition,Nb_pts=self.Nb_pts)

        cst = -codata.e / (codata.m_e * gamma)
        initial_condition_for_ODE=self.initial_condition
        res = odeint(fct_ODE_magnetic_field,initial_condition_for_ODE, trajectory[0],
                     args=(cst,B.Bx,B.By,B.Bz), full_output=True)
        traj = res[0]
        info = res[1]
        # print("1 : nonstiff problems, Adams . 2: stiff problem, BDF")
        #print(info.get('nst'))
        traj = np.transpose(traj)
        trajectory[4] = traj[0]
        trajectory[5] = traj[1]
        trajectory[6] = traj[2]
        trajectory[1] = traj[3]
        trajectory[2] = traj[4]
        trajectory[3] = traj[5]
        trajectory[7] = - cst * B.By(trajectory[3], trajectory[2],trajectory[1]) * trajectory[6]
        trajectory[9] = cst * B.By(trajectory[3], trajectory[2],trajectory[1]) * trajectory[4]
        T=self.create_from_array(trajectory)
        T.multiply_by((1.0/codata.c))
        return T

    def create_from_source(self, source):
        if (self.method == TRAJECTORY_METHOD_ODE):
            if (self.initial_condition is None):
                self.choise_initial_condition(source=source)
            trajectory = self.trajectory_from_magnetic_field_method_ODE(source=source)
        else:
            if source.magnet_type()==PLANE_UNDULATOR:
                trajectory = self.analytical_trajectory_plane_undulator(undulator=source)
            else :
                trajectory = self.analytical_trajectory_cst_magnf(bending_magnet=source)
            self.initial_condition = np.array([trajectory.v_x[0], trajectory.v_y[0], trajectory.v_z[0],
                                                   trajectory.x[0], trajectory.y[0], trajectory.z[0]])
            self.initial_condition *= codata.c

        return trajectory

    def create_from_array(self,array):
        if array.shape[0] != 10 :
            raise Exception('this array can not be convert in Trajectory')
        return Trajectory(t=array[0],x=array[1],y=array[2],z=array[3],v_x=array[4],v_y=array[5],
                    v_z = array[6],a_x=array[7],a_y=array[8],a_z=array[9])

    def get_method(self):
        if self.method==TRAJECTORY_METHOD_ANALYTIC :
            method='Analitical Trajectory'

        elif self.method == TRAJECTORY_METHOD_ODE:
            method = ' ODE solution '

        else :# self.method == TRAJECTORY_METHOD_Integration:
            method = ' Trajectory from integration on the magnetic field'
        return method


    def print_parameters(self):
        print("Trajectory ")
        print('    method : %s' %self.get_method())
        print('    number of points : %d' %self.Nb_pts)
        print('    initial position (x,y,z) : ')
        print(self.initial_condition[3:6])
        print('    initial velocity (x,y,z)  ')
        print(self.initial_condition[0:3])


if __name__ == "__main__" :
    undulator_test = Undulator(K=1.87, period_length=0.035, length=0.035 * 14)
    electron_beam_test = ElectronBeam(Electron_energy=1.3e9, I_current=1.0)

    source_test=SourceUndulatorPlane(undulator=undulator_test,electron_beam=electron_beam_test)


    print('Create trajectory with autamatic choice of initial condition and automatic magnetic field')

    trajectory_fact_ODE = TrajectoryFactory(Nb_pts=20000, method=TRAJECTORY_METHOD_ODE)
    trajectory1=trajectory_fact_ODE.create_from_source(source_test)
    print(' ')
    print('trajectory 1 created with ODE method')
    print(trajectory_fact_ODE.print_parameters())
    print('initial condition =')
    print(trajectory_fact_ODE.initial_condition * (1. / codata.c))
    trajectory1.plot_3D(title="traj 1: autamatic initial condition and field")

    trajectory_fact_ANALITIC = TrajectoryFactory(Nb_pts=20000, method=TRAJECTORY_METHOD_ANALYTIC)
    trajectory2=trajectory_fact_ANALITIC.create_from_source(source_test)
    print(' ')
    print('trajectory 2 create with ANALYTIC method')
    print(trajectory_fact_ANALITIC.print_parameters())
    print('initial condition =')
    print(trajectory_fact_ANALITIC.initial_condition*(1./codata.c))
    trajectory2.plot_3D(title="traj 2: ANALYTICAL method")


    print(' ')
    print('initial condition can be impose:')
    print('trajectory 1 modified with initial condition of trajectory2')
    trajectory_fact_ODE.initial_condition=trajectory_fact_ANALITIC.initial_condition
    #print(all(trajectory_fact_ODE.initial_condition==trajectory_fact_ANALITIC.initial_condition))
    trajectory1=trajectory_fact_ODE.create_from_source(source_test)
    trajectory1.plot_3D("traj 1 modified with intitial conditions from traj 2")

    print(' ')
    print('Now trajectory 1 and 2 have the same vector time, we make the difference : ')
    diff =trajectory1.difference_with(trajectory2)
    diff.plot(title="Difference")