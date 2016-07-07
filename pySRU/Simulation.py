import numpy as np
import matplotlib.pyplot as plt
from  scipy.constants import physical_constants
from mpl_toolkits.mplot3d import Axes3D
import scipy.constants as codata
import time
from scipy.interpolate import interp1d
from TrajectoryAnalitic import TrajectoryAnalitic
from TrajectoryArray import TrajectoryArray
from Radiation import Radiation
from MagneticField import MagneticField
from ParameterPlaneUndulator import ParameterPlaneUndulator as Undulator ,PLANE_UNDULATOR,BENDING_MAGNET

from TrajectoryFactory import TrajectoryFactory, TRAJECTORY_METHOD_ANALYTIC,TRAJECTORY_METHOD_ODE,\
                                        TRAJECTORY_METHOD_INTEGRATION
from RadiationFactoryAnalitic import RadiationFactoryAnalitic
from RadiationFactory import RadiationFactory,RADIATION_METHOD_NEAR_FIELD, \
                                RADIATION_METHOD_FARFIELD, RADIATION_METHOD_APPROX_FARFIELD


eV_to_J=1.602176487e-19

class Simulation(object):

    def __init__(self, parameter, trajectory_fact, magnetic_field, radiation_fact, trajectory, radiation):
        self.parameter=parameter
        self.magnetic_field=magnetic_field
        self.trajectory_fact=trajectory_fact
        self.radiation_fact=radiation_fact
        self.trajectory=trajectory
        self.radiation=radiation

# faire un autre constructeur pour la copie


    def copy(self):
        if self.magnetic_field==None :
            mag_f=None
        else :
            mag_f =self.magnetic_field.copy()
        return Simulation(parameter=self.parameter.copy(), trajectory_fact=self.trajectory_fact.copy(),
                                   magnetic_field=mag_f, radiation_fact=self.radiation_fact.copy(),
                                   trajectory=self.trajectory.copy(), radiation=self.radiation.copy())

    #change

    def change_distance(self,D):
        self.radiation.distance=D
        #update intensity
        self.radiation.intensity = self.radiation_fact.calculate_radiation_intensity(trajectory=self.trajectory,
                                                                                     parameter=self.parameter,
                                                                                     distance=self.radiation.distance,
                                                                                     X_arrays=self.radiation.X,
                                                                                     Y_arrays=self.radiation.Y)

    def change_radiation_method(self,method):
        self.radiation_fact.method=method
        #update intensity
        self.radiation.intensity= self.radiation_fact.calculate_radiation_intensity(trajectory=self.trajectory,
                                                                                    parameter=self.parameter,
                                                                                    distance=self.radiation.distance,
                                                                                    X_arrays=self.radiation.X,
                                                                                    Y_arrays=self.radiation.Y)

    def change_trajectory_method(self,method) :
        self.trajectory_fact.method = method
        self.trajectory =self.trajectory_fact.create_for_parameter(parameter=self.parameter,
                                                                         B=self.magnetic_field)
        self.radiation.intensity = self.radiation_fact.calculate_radiation_intensity(trajectory=self.trajectory,
                                                                                     parameter=self.parameter,
                                                                                     distance=self.radiation.distance,
                                                                                     X_arrays=self.radiation.X,
                                                                                     Y_arrays=self.radiation.Y)

    def change_initial_condition(self, initial_cond):
        self.trajectory_fact.initial_condition= initial_cond
        self.magnetic_field.change_Zo(Zo=initial_cond[5])
        self.trajectory = self.trajectory_fact.create_for_parameter(parameter=self.parameter,
                                                                    B=self.magnetic_field)

        self.radiation.intensity = self.radiation_fact.calculate_radiation_intensity(trajectory=self.trajectory,
                                                                                     parameter=self.parameter,
                                                                                     distance=self.radiation.distance,
                                                                                     X_arrays=self.radiation.X,
                                                                                     Y_arrays=self.radiation.Y)

    def change_omega(self, omega) :
        self.radiation_fact.omega=omega
        self.radiation.intensity = self.radiation_fact.calculate_radiation_intensity(trajectory=self.trajectory,
                                                                                     parameter=self.parameter,
                                                                                     distance=self.radiation.distance,
                                                                                     X_arrays=self.radiation.X,
                                                                                     Y_arrays=self.radiation.Y)

    def change_energy_eV(self, E) :
        omega = E * eV_to_J / codata.hbar
        self.change_omega(omega)

    # only the trajectory change, not the radiation
    # must be use in special case like method "time_radiation"
    def change_Nb_pts_trajectory_only(self, Nb_pts):
        self.trajectory_fact.Nb_pts = Nb_pts
        self.trajectory = self.trajectory_fact.create_for_plane_undulator(undulator=self.parameter,
                                                                          B=self.magnetic_filed)

    def change_Nb_pts_trajectory(self, Nb_pts):
        self.trajectory_fact.Nb_pts = Nb_pts
        self.trajectory = self.trajectory_fact.create_for_parameter(parameter =self.parameter,
                                                                          B=self.magnetic_field)
        self.radiation.intensity = self.radiation_fact.calculate_radiation_intensity(trajectory=self.trajectory,
                                                                                     parameter=self.parameter,
                                                                                     distance=self.radiation.distance,
                                                                                     X_arrays=self.radiation.X,
                                                                                     Y_arrays=self.radiation.Y)

    def change_Nb_pts_radiation(self, Nb_pts):
        self.radiation_fact.Nb_pts = Nb_pts
        self.radiation.X = np.linspace(self.radiation.X[0], self.radiation.X[-1], Nb_pts)
        self.radiation.Y = np.linspace(self.radiation.Y[0], self.radiation.Y[-1], Nb_pts)
        self.radiation.intensity = self.radiation_fact.calculate_radiation_intensity(trajectory=self.trajectory,
                                                                                     parameter=self.parameter,
                                                                                     distance=self.radiation.distance,
                                                                                     X_arrays=self.radiation.X,
                                                                                     Y_arrays=self.radiation.Y)

    # error

    def error_radiation_method_nb_pts_traj(self,method,nb_pts):
        sim2=self.copy()
        sim2.change_radiation_method(method)
        error=np.zeros_like(nb_pts)
        print(len(nb_pts))
        for i in range(len(nb_pts)) :
            print(i)
            self.change_Nb_pts_trajectory(nb_pts[i])
            sim2.change_Nb_pts_trajectory(nb_pts[i])
            error[i]=self.radiation.error_max(sim2.radiation)
        return error

    def error_radiation_method_distance(self,method,D):
        sim2=self.copy()
        sim2.change_radiation_method(method)
        error=np.zeros_like(D)
        print(len(D))
        for i in range(len(D)) :
            print(i)
            self.change_distance(D[i])
            sim2.change_distance(D[i])
            error[i]=self.radiation.error_max(sim2.radiation)
        return error

    def error_trajectory_method(self,method,nb_pts):
        sim2=self.copy()
        sim2.change_trajectory_method(method)
        error_rad=np.zeros_like(nb_pts)
        error_traj=np.zeros((10,len(nb_pts)))
        for i in range(len(nb_pts)) :
            self.change_Nb_pts_trajectory(nb_pts[i])
            sim2.change_Nb_pts_trajectory(nb_pts[i])
            error_traj[:,i]=self.trajectory.error_max(sim2.trajectory)
            error_rad[i]=self.radiation.error_max(sim2.radiation)
        traj_error_traj=self.trajectory_fact.create_from_array(error_traj)
        return error_rad, traj_error_traj

    # a changer #TODO
    def error_trajectory_cst_magnetic_field(self,nb_pts):
        Bo = (self.parameter.K / (93.4 * self.parameter.lambda_u))
        f = Bo * codata.e / (self.parameter.gamma() * codata.m_e)
        vz_0=self.trajectory_fact.initial_condition[2]/codata.c
        x_0=self.trajectory_fact.initial_condition[3]/codata.c
        error = np.zeros((10,len(nb_pts)))
        traj_ref=self.trajectory.copy()
        for i in range(len(nb_pts)):
            self.change_Nb_pts_trajectory_only(nb_pts[i])
            traj_ref_arrays= self.trajectory_fact.analytical_trajectory_cst_magnf(f=f,
                                                    t=self.trajectory.t,vz_0=vz_0,x_0=x_0)
            traj_ref.convert(traj_ref_arrays)
            error_max=self.trajectory.error_max(traj_ref)
            error[:,i]=error_max
            # print(error_max.shape)
            # for j in range(9) :
            #      error[j+1][i] = error_max[j+1]
        traj_err=self.trajectory_fact.create_from_array(error)
        return traj_err

    def relativ_error_radiation_method_distance(self,method,D):
        sim2=self.copy()
        sim2.change_radiation_method(method)
        error_relativ=np.array(len(D))
        for i in range(len(D)) :
            self.change_distance(D[i])
            sim2.change_distance(D[i])
            error_relativ[i]=self.radiation.relativ(sim2.radiation)
        return error_relativ

    #spectre

    def spectre(self,omega_array=None) :
        if omega_array == None :
            omega1=self.parameter.omega1()
            omega_array=np.arange(omega1*0.9,5.0*omega1*1.1,omega1*0.01)
        spectre=np.zeros_like(omega_array)
        print(len(spectre))
        # print('self.radiation.X.shape')
        # print(self.radiation.X.shape)
        # print('self.radiation.Y.shape')
        # print(self.radiation.Y.shape)
        for i in range(len(spectre)) :
            print(i)
            self.change_omega(omega_array[i])
            spectre[i]=self.radiation.integration()
        return spectre , omega_array

    def spectre2(self,omega_array=None) :
        if omega_array == None :
            omega1=self.parameter.omega1()
            omega_array=np.arange(omega1*0.9,5.0*omega1*1.1,omega1*0.01)
        spectre=np.zeros_like(omega_array)
        print(len(spectre))
        # print('self.radiation.X.shape')
        # print(self.radiation.X.shape)
        # print('self.radiation.Y.shape')
        # print(self.radiation.Y.shape)
        for i in range(len(spectre)) :
            print(i)
            self.change_omega(omega_array[i])
            spectre[i]=self.radiation.max()
        return spectre , omega_array

    def spectre3(self,omega_array=None) :
        if omega_array == None :
            omega1=self.parameter.omega1()
            omega_array=np.arange(omega1*0.9,5.0*omega1*1.1,omega1*0.01)
        spectre=np.zeros_like(omega_array)
        print(len(spectre))
        # print('self.radiation.X.shape')
        # print(self.radiation.X.shape)
        # print('self.radiation.Y.shape')
        # print(self.radiation.Y.shape)
        self.radiation.X=np.array([0.0])
        self.radiation.Y = np.array([0.0])
        for i in range(len(spectre)) :
            self.change_omega(omega_array[i])
            spectre[i]=self.radiation.max()
        return spectre , omega_array

    def spectre_max(self, omega_array=None):
        start_time = time.time()
        spectre,omega_array = self.spectre3(omega_array=omega_array)
        interval = time.time() - start_time
        # print("interval temps :")
        # print(interval)
        omega1=self.parameter.omega1()
        omega1_min=np.ceil(omega_array[0]/omega1)
        # print("harmonic min =")
        # print(omega1_min)
        omega1_max=np.floor(omega_array[-1]/omega1)
        # print("harmonic max =")
        # print(omega1_max)
        omega1_array = omega1 * np.ones_like(omega_array)
        plt.plot(omega_array, spectre)
        i=omega1_min
        while i <= omega1_max :
            plt.plot(i * omega1_array, spectre, color="red")
            i += 1
        print("omega1 =")
        print(omega1)
        omega_max = omega_array[np.argmax(spectre)]
        print("omega max =")
        print(omega_max)
        plt.show()
        return spectre , omega_array , omega_max

    # time

    def time_radiation(self,Nb_pts_trajectory,Nb_pts_radiation):
        N=len(Nb_pts_trajectory)
        M=len(Nb_pts_radiation)
        calc_time=np.zeros((N,M))
        for i in range(N) :
            print(i)
            self.change_Nb_pts_trajectory_only(Nb_pts_trajectory[i])
            for j in range(M) :
                start_time=time.time()
                self.change_Nb_pts_radiation(Nb_pts_radiation[j])
                calc_time[i,j]=time.time()-start_time
        return  calc_time






def create_simulation(parameter, trajectory_fact, radiation_fact=None, magnetic_field=None, X_max=None, Y_max=None,
                      distance=None):

    if radiation_fact==None :
        radiation_fact=RadiationFactory(method=RADIATION_METHOD_APPROX_FARFIELD,omega=parameter.omega1()
                                        ,Nb_pts=101)

    if (trajectory_fact.initial_condition == None and trajectory_fact.method != TRAJECTORY_METHOD_ANALYTIC) :

        trajectory_fact.initial_condition = np.array([0.0, 0.0, parameter.Beta() * codata.c,
                                                      0.0, 0.0, parameter.Zo_symetry()])
        # a changer
    elif(trajectory_fact.method==TRAJECTORY_METHOD_ANALYTIC):
        trajectory_fact.initial_condition = np.array([0.0, 0.0, parameter.Beta() * codata.c,
                                                      0.0, 0.0,  parameter.Zo_analitic()])

    Zo= trajectory_fact.initial_condition[5]
    Yo=trajectory_fact.initial_condition[4]

    if magnetic_field == None :
        if Zo== 0.0 :
            print("Warning : non symetriq magnet")
            Z=np.linspace(0.0,parameter.Zmax_no_symetry,trajectory_fact.Nb_pts)
        else :
            Z = np.linspace(Zo,-Zo, trajectory_fact.Nb_pts)
        if Yo!=0.0 :
            Y=np.linspace(-Yo,Yo,trajectory_fact.Nb_pts)
        else :
            Y=0.0
        harmonic_number=np.floor(radiation_fact.omega/parameter.omega1())

        ### A changer !!!!!!!!!!!
        if (harmonic_number==0) :
            harmonic_number=1

        magnetic_field = parameter.create_magnetic_field(Z=Z,Y=Y,X=0.0,
                                                harmonic_number=harmonic_number)

    else:
        if (magnetic_field.z != np.ndarray):
            if Zo == 0.0:
                Z = np.linspace(0.0, parameter.Zmax_no_symetry, trajectory_fact.Nb_pts)
            else:
                Z = np.linspace(Zo, -Zo, trajectory_fact.Nb_pts)
            magnetic_field.z=Z
        elif type(magnetic_field.By)==np.ndarray : #and trajectory_fact.method==TRAJECTORY_METHOD_ODE :
            magnetic_field.enlargement_vector_for_interpolation(nb_enlarg=np.floor(len(magnetic_field.z) * 0.1))
            magnetic_field = interp1d(magnetic_field.z, magnetic_field.By)


    trajectory = trajectory_fact.create_for_parameter(parameter=parameter, B=magnetic_field)
    Nb_pts=radiation_fact.Nb_pts
    if X_max ==None or Y_max== None :
        if (X_max != None) :
            X=np.linspace(0.0,X_max,Nb_pts)
            Y=np.linspace(0.0,X_max,Nb_pts)
        elif Y_max != None :
            X = np.linspace(0.0, Y_max, Nb_pts)
            Y = np.linspace(0.0, Y_max, Nb_pts)
        else :
            X=None
            Y=None
    else :
        X = np.linspace(0.0, X_max, Nb_pts)
        Y = np.linspace(0.0, Y_max, Nb_pts)
#
    #en attendant de faire mieux
    if type(trajectory)==TrajectoryAnalitic and type(radiation_fact == RadiationFactoryAnalitic):
        # radiation_fact=RadiationFactoryAnalitic(formula=radiation_fact.formula,omega=radiation_fact.omega,
        #                                          method=radiation_fact.method,Nb_pts=radiation_fact.Nb_pts )
        # radiation=radiation_fact.create_for_single_electron(trajectory=trajectory,
        #                                                    parameter=parameter,
        #                                                    distance=distance,
        #                                                    X=X, Y=Y)
        trajectory_array=trajectory.convert()

    else :
        trajectory_array=trajectory
    radiation = radiation_fact.create_for_single_electron(trajectory=trajectory_array,  parameter=parameter,
                                                          distance=distance,   X=X, Y=Y)

    return Simulation(parameter=parameter, trajectory_fact=trajectory_fact, magnetic_field=magnetic_field,
                               radiation_fact=radiation_fact, trajectory=trajectory_array, radiation=radiation)




if __name__ == "__main__" :
    und_test = Undulator(K=1.87, E=1.3e9, lambda_u=0.035, L=0.035 * 12, I=1.0)
    #initial_cond = np.array(
        #[0.0, 0.0, (und_test.Beta() * codata.c), 0.0, 0.0, -und_test.L / 2.0 - 5.0 * und_test.lambda_u])
    traj_test = TrajectoryFactory(Nb_pts=1001, method=TRAJECTORY_METHOD_ODE)
    #rad_test = RadiationFactory(method=RADIATION_METHOD_APPROX_FARFIELD, omega=und_test.omega1(), Nb_pts=101)
    distance=100
    # X = np.linspace(0.0, distance*1.01e-3, 101)
    # Y = np.linspace(0.0,distance*1.01e-3, 101)
    Xmax=distance*1e-3
    Ymax=distance*1e-3

    sim_test = create_simulation(parameter=und_test, trajectory_fact=traj_test,
                                 distance=distance, X_max=Xmax, Y_max=Ymax)




