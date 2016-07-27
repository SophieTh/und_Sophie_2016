import numpy as np
import scipy.constants as codata
import scipy.integrate as integrate
from scipy.integrate import odeint
from pySRU.Trajectory import Trajectory
from pySRU.Source import Source,PLANE_UNDULATOR,BENDING_MAGNET
from pySRU.MagneticStructureUndulatorPlane import MagneticStructureUndulatorPlane as Undulator
from pySRU.ElectronBeam import ElectronBeam


TRAJECTORY_METHOD_ANALYTIC=0
TRAJECTORY_METHOD_ODE=1
TRAJECTORY_METHOD_INTEGRATION=2


def fct_ODE_magnetic_field(y, t, cst, Bx,By,Bz):
    return [cst * (Bz(z=y[5]*codata.c,y=y[4]*codata.c,x=y[3]*codata.c) * y[1]
                   - By(z=y[5]*codata.c,y=y[4]*codata.c,x=y[3]*codata.c) * y[2]),
            cst * (Bx(z=y[5]*codata.c,y=y[4]*codata.c,x=y[3]*codata.c) * y[2]
                   - Bz(z=y[5]*codata.c,y=y[4]*codata.c,x=y[3]*codata.c) * y[0]),
            cst * (By(z=y[5]*codata.c,y=y[4]*codata.c,x=y[3]*codata.c) * y[0]
                   - Bx(z=y[5]*codata.c,y=y[4]*codata.c,x=y[3]*codata.c) * y[1]),
            y[0],
            y[1],
            y[2]]

def fct_ODE_magnetic_field2(y, t, cst, Bx,By,Bz):
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

    def choise_initial_condition(self,source):
        if self.method != TRAJECTORY_METHOD_ANALYTIC:
            self.initial_condition = np.array([0.0, 0.0, source.Beta() * codata.c,
                                                      0.0, 0.0, source.Zo_symetry()])

    # calculate a theorical trajectory in an undulator
    def analytical_trajectory_plane_undulator(self,undulator):
        N= self.Nb_pts
        ku = 2.0 * np.pi / undulator.lambda_u()
        gamma = undulator.gamma()
        Beta_et = undulator.Beta_et()
        K=undulator.K()
        omega_u = Beta_et * codata.c * ku

        t = np.linspace(undulator.Zo_analitic(),-undulator.Zo_analitic(), N)*(1./(Beta_et*codata.c))
        z = Beta_et * t + ((K / gamma) ** 2) * (1.0 / (8.0 * omega_u)) * np.sin( 2.0 * omega_u*t)
        x = (-(K / (gamma * omega_u)) * np.cos(omega_u*t))
        # # Vx and Vz
        v_z = Beta_et + ((K / gamma) ** 2) * (1.0 / 4.0) * np.cos(2.0 *omega_u*t)
        v_x= (K / (gamma )) * np.sin(omega_u*t)
        # # Ax and Az
        a_z=-omega_u *(K / gamma) ** 2 * 0.5 * np.sin( 2.0 * omega_u*t)
        a_x= (K / (gamma )) * (omega_u ) * np.cos(omega_u*t)
        y=0.0*t
        v_y=y
        a_y=y
        return Trajectory(t=t,x=x,y=y,z=z,v_x=v_x,v_y=v_y,v_z=v_z,
                                  a_x=a_x,a_y=a_y,a_z=a_z)

    def analytical_trajectory_cst_magnf(self,source):
        t = np.linspace(source.Zo_analitic(),-source.Zo_analitic(),
                           self.Nb_pts) / ( source.Beta_et()* codata.c)
        f = -source.Bo() * codata.e / (source.gamma() * codata.m_e)
        vz_0 = self.initial_condition[2]/codata.c
        t0 = source.Zo_analitic()/(codata.c*source.Beta_et())
        x = -vz_0 * (1.0 / f) * np.cos(f * (t-t0))+vz_0/f+self.initial_condition[3]/codata.c
        y = vz_0 * 0.0 * (t-t0)+self.initial_condition[4]/codata.c
        z = vz_0 * (1.0 / f) * np.sin(f * (t-t0))+self.initial_condition[5]/codata.c
        vx = vz_0 * np.sin(f * (t-t0))+self.initial_condition[0]/codata.c
        vy = vz_0 * 0.0 * (t-t0)+self.initial_condition[1]/codata.c
        vz = vz_0 * np.cos(f * (t-t0))+self.initial_condition[2]/codata.c
        ax = vz_0 * f * np.cos(f * (t-t0))
        ay = vz_0 * 0.0 * (t-t0)
        az = -vz_0* f * np.sin(f * (t-t0))
        return Trajectory(t=t,x=x,y=y,z=z,v_x=vx,v_y=vy,v_z=vz,a_x=ax,a_y=ay,a_z=az)
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
    def trajectory_from_magnetic_field_method_INT(self,source):
        gamma=source.gamma()
        Beta = source.Beta()
        Beta_et = source.Beta_et()
        N=self.Nb_pts
        B=source.magnetic_field
        # faire un truc create Z or t pour INT et ODE
        Zo=self.initial_condition[5]
        Z=np.linspace(Zo,-Zo,N)
        Y=0.0
        X=0.0
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
            trajectory[4][i] = integrate.simps(trajectory[7][0:(i + 1)], trajectory[0][0:(i + 1)],even='first') \
                                       + self.initial_condition[0]/ codata.c
        trajectory[6] = np.sqrt((Beta) ** 2 - trajectory[4] ** 2)
        # X et Z
        for i in range(N):
            #trajectory[1][i] = integrate.simps(trajectory[4][0:(i + 1)], trajectory[0][0:(i + 1)]) \
            trajectory[1][i] =  integrate.simps(trajectory[4][0:(i + 1)], trajectory[0][0:(i + 1)],even='first') \
                               + self.initial_condition[3]/ codata.c
            #trajectory[3][i] = integrate.simps(trajectory[6][0:(i + 1)], trajectory[0][0:(i + 1)]) \
            trajectory[3][i] = integrate.simps(trajectory[6][0:(i + 1)], trajectory[0][0:(i + 1)],even='first') \
                               +  self.initial_condition[5]/ codata.c

            # Az
        trajectory[9] = -(trajectory[7] * trajectory[4]) / trajectory[6]

        T=self.create_from_array(trajectory)
        return T

    def trajectory_from_magnetic_field_method_ODE(self, source):
        gamma = source.gamma()
        Beta_et = source.Beta_et()
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
        Zo=self.initial_condition[5]
        Z=np.linspace(Zo,-Zo,self.Nb_pts)
        trajectory[0] = Z / (Beta_et * codata.c)

        # trajectory[0] = np.linspace(B.z[0] / (Beta_et * codata.c),
        #                             B.z[-1] / (Beta_et * codata.c), self.Nb_pts)
        cst = -codata.e / (codata.m_e * gamma)
        initial_condition_for_ODE=self.initial_condition*(1.0/codata.c)
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
        #T.multiply_by((1.0/codata.c))
        return T

    def create_from_source(self, source):
        if (self.method == TRAJECTORY_METHOD_INTEGRATION or self.method == TRAJECTORY_METHOD_ODE):
            if (self.initial_condition == None):
                self.choise_initial_condition(source=source)
            if self.method == TRAJECTORY_METHOD_INTEGRATION:
                trajectory = self.trajectory_from_magnetic_field_method_INT(source=source)
            else:  # method=TRAJECTORY_METHOD_ODE
                trajectory = self.trajectory_from_magnetic_field_method_ODE(source=source)
        else:
            print(source.magnet_type())
            if source.magnet_type()==PLANE_UNDULATOR:
                trajectory = self.analytical_trajectory_plane_undulator(undulator=source)
                self.initial_condition = np.array([trajectory.v_x[0],trajectory.v_y[0],trajectory.v_z[0],
                                                   trajectory.x[0], trajectory.y[0], trajectory.z[0]])
                self.initial_condition *= codata.c
            else :
                if (self.initial_condition == None):
                    self.initial_condition = np.array([0.0, 0.0, source.Beta() * codata.c,
                                                       0.0, 0.0, source.Zo_analitic()])
                trajectory = self.analytical_trajectory_cst_magnf(bending_magnet=source)

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
        print '    method : %s' %self.get_method()
        print('    number of points : %d' %self.Nb_pts)
        print('    initial position (x,y,z) : %.3f '
              %self.initial_condition[3] )
        print('    initial velocity (x,y,z) : %.3f '
              % self.initial_condition[0])


if __name__ == "__main__" :
    undulator_test = Undulator(K=1.87, lambda_u=0.035, L=0.035 * 14)
    electron_beam_test = ElectronBeam(E=1.3e9, I=1.0)
    source_test=Source(magnetic_structure=undulator_test,electron_beam=electron_beam_test)


    print('Create trajectory with autamatic choise of initial condition and automatic magnetic field')

    trajectory_fact_ODE = TrajectoryFactory(Nb_pts=2000, method=TRAJECTORY_METHOD_ODE)
    trajectory1=trajectory_fact_ODE.create_from_source(source_test)
    print(' ')
    print('trajectory 1 create with ODE method')
    print('initial condition =')
    print(trajectory_fact_ODE.initial_condition * (1. / codata.c))
    trajectory1.plot_3D()

    trajectory_fact_ANALITIC = TrajectoryFactory(Nb_pts=2000, method=TRAJECTORY_METHOD_ANALYTIC)
    trajectory2=trajectory_fact_ANALITIC.create_from_source(source_test)
    print(' ')
    print('trajectory 2 create with ANALYTIC method')
    print('initial condition =')
    print(trajectory_fact_ANALITIC.initial_condition*(1./codata.c))
    trajectory2.plot_3D()


    print(' ')
    print('initial condition can be impose:')
    print('trajectory 1 modified with initial condition of trajectory2')
    trajectory_fact_ODE.initial_condition=trajectory_fact_ANALITIC.initial_condition
    #print(all(trajectory_fact_ODE.initial_condition==trajectory_fact_ANALITIC.initial_condition))
    trajectory1=trajectory_fact_ODE.create_from_source(source_test)
    trajectory1.plot_3D()

    print(' ')
    print('Now trajectory 1 and 2 have the same vector time, we make the difference : ')
    diff =trajectory1.difference_with(trajectory2)
    diff.plot()