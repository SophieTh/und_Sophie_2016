import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import scipy.constants as codata
import time
from scipy.interpolate import interp1d
from Trajectory import Trajectory
from Radiation import Radiation
from MagneticField import MagneticField
from UndulatorParameter import UndulatorParameters as Undulator
from TrajectoryFactory import TrajectoryFactory, TRAJECTORY_METHOD_ANALYTIC,TRAJECTORY_METHOD_ODE,\
                                        TRAJECTORY_METHOD_INTEGRATION
from RadiationFactory import RadiationFactory ,RADIATION_METHOD_NEAR_FIELD, \
                                RADIATION_METHOD_FARFIELD, RADIATION_METHOD_APPROX_FARFIELD



class UndulatorSimulation(object):

    def __init__(self,undulator,trajectory_fact,magnetic_field,radiation_fact,trajectory,radiation):
        self.undulator=undulator
        self.magnetic_filed=magnetic_field
        self.trajectory_fact=trajectory_fact
        self.radiation_fact=radiation_fact
        self.trajectory=trajectory
        self.radiation=radiation

# faire un autre constructeur pour la copie


    def copy(self):
        if self.magnetic_filed==None :
            mag_f=None
        else :
            mag_f =self.magnetic_filed.copy()
        return UndulatorSimulation(undulator=self.undulator.copy(),trajectory_fact=self.trajectory_fact.copy(),
                                   magnetic_field=mag_f,radiation_fact=self.radiation_fact.copy(),
                                   trajectory=self.trajectory.copy(),radiation=self.radiation.copy())

    def change_distance(self,D):
        self.radiation.distance=D
        #update intensity
        self.radiation.intensity = self.radiation_fact.calculate_radiation_intensity(trajectory=self.trajectory,
                                                                                     undulator=self.undulator,
                                                                                    distance=self.radiation.distance,
                                                                                    X_arrays=self.radiation.X,
                                                                                    Y_arrays=self.radiation.Y)


    def change_radiation_method(self,method):
        self.radiation_fact.method=method
        #update intensity
        self.radiation.intensity= self.radiation_fact.calculate_radiation_intensity(trajectory=self.trajectory,
                                                                                     undulator=self.undulator,
                                                                                    distance=self.radiation.distance,
                                                                                    X_arrays=self.radiation.X,
                                                                                    Y_arrays=self.radiation.Y)

    def change_trajectory_method(self,method) :
        self.trajectory_fact.method = method
        self.trajectory =self.trajectory_fact.create_for_plane_undulator(undulator=self.undulator,
                                                                         Z_By=self.magnetic_filed)

        self.radiation.intensity = self.radiation_fact.calculate_radiation_intensity(trajectory=self.trajectory,
                                                                                     undulator=self.undulator,
                                                                                 distance=self.radiation.distance,
                                                                                 X_arrays=self.radiation.X,
                                                                                 Y_arrays=self.radiation.Y)

    def change_omega(self, omega) :
        self.radiation_fact.omega=omega
        self.radiation.intensity = self.radiation_fact.calculate_radiation_intensity(trajectory=self.trajectory,
                                                    undulator=self.undulator, distance=self.radiation.distance,
                                                    X_arrays=self.radiation.X,Y_arrays=self.radiation.Y)

    def error_radiation_method(self,method,D):
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

    def relativ_error_radiation_method(self,method,D):
        sim2=self.copy()
        sim2.change_radiation_method(method)
        error_relativ=np.array(len(D))
        for i in range(len(D)) :
            self.change_distance(D[i])
            sim2.change_distance(D[i])
            error_relativ[i]=self.radiation.relativ(sim2.radiation)
        return error_relativ

    def spectre(self,omega_array=None) :
        if omega_array == None :
            omega1=self.undulator.omega1()
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
            omega1=self.undulator.omega1()
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
            omega1=self.undulator.omega1()
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
        omega1=self.undulator.omega1()
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

    # only the trajectory change, not the radiation
    # must be use in special case like method "time_radiation"
    def change_Nb_pts_trajectory_only(self,Nb_pts) :
        self.trajectory_fact.Nb_pts = Nb_pts
        self.trajectory=self.trajectory_fact.create_for_plane_undulator(undulator=self.undulator,
                                                                         B=self.magnetic_filed)


    def change_Nb_pts_trajectory(self,Nb_pts) :
        self.trajectory_fact.Nb_pts = Nb_pts
        self.trajectory=self.trajectory_fact.create_for_plane_undulator(undulator=self.undulator,
                                                                         B=self.magnetic_filed)
        self.radiation.intensity = self.radiation_fact.calculate_radiation_intensity(trajectory=self.trajectory,
                                                                                 undulator=self.undulator,
                                                                                 distance=self.radiation.distance,
                                                                                 X_arrays=self.radiation.X,
                                                                                 Y_arrays=self.radiation.Y)


    def change_Nb_pts_radiation(self,Nb_pts) :
        self.radiation_fact.Nb_pts=Nb_pts
        self.radiation.X=np.linspace(self.radiation.X[0],self.radiation.X[-1]*1.00001,Nb_pts)
        self.radiation.Y = np.linspace(self.radiation.Y[0], self.radiation.Y[-1] * 1.00001, Nb_pts)
        self.radiation.intensity = self.radiation_fact.calculate_radiation_intensity(trajectory=self.trajectory,
                                                                                     undulator=self.undulator,
                                                                                     distance=self.radiation.distance,
                                                                                     X_arrays=self.radiation.X,
                                                                                     Y_arrays=self.radiation.Y)

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







def create_simulation(undulator, trajectory_fact, radiation_fact=None, magnetic_field=None, X=None, Y=None,
                      distance=None):

    if (radiation_fact == None):
        if X== None or Y== None :
            Nb_pts=101
        else :
            Nb_pts=len(X)

        radiation_fact = RadiationFactory(method=RADIATION_METHOD_APPROX_FARFIELD,
                                          omega=undulator.omega1(),Nb_pts=Nb_pts)
    else:
        radiation_fact = radiation_fact

    if (trajectory_fact.initial_condition == None):
        Zo=-(undulator.L / 2.0 + 5.0 * undulator.lambda_u)
        trajectory_fact.initial_condition = np.array([0.0, 0.0,
                                np.sqrt(1.0 - (1.0 / (undulator.E / 0.511e6) ** 2)) * codata.c,
                                           0.0, 0.0, -Zo])
    Zo=-trajectory_fact.initial_condition[5]
    Yo=-trajectory_fact.initial_condition[4]

    if magnetic_field == None:
        Z = np.linspace(-Zo,
                        Zo, radiation_fact.Nb_pts)
        if Yo!=0.0 :
            Y=np.linspace(-Yo,Yo,radiation_fact.Nb_pts)
        else :
            Y=0.0
        harmonic_number=np.floor(radiation_fact.omega/undulator.omega1())
        ### A changer !!!!!!!!!!!
        if (harmonic_number==0) :
            harmonic_number=1
        magnetic_field = undulator.create_magnetic_field_plane_undulator(Z=Z,Y=Y,
                                            harmonic_number=harmonic_number)
    else:
        if type(magnetic_field.By)==np.ndarray :
            magnetic_field.enlargement_vector_for_interpolation(nb_enlarg=np.floor(len(magnetic_field.z) * 0.1))
            magnetic_field = interp1d(magnetic_field.z, magnetic_field.By)
        else :
            print('type de magnetic field inconnu')

    trajectory = trajectory_fact.create_for_plane_undulator(undulator=undulator, B=magnetic_field)

    radiation = radiation_fact.create_for_single_electron(trajectory=trajectory,
                                                          undulator=undulator,
                                                          distance=distance,
                                                          X=X, Y=Y)

    return UndulatorSimulation(undulator=undulator, trajectory_fact=trajectory_fact, magnetic_field=magnetic_field,
                               radiation_fact=radiation_fact, trajectory=trajectory, radiation=radiation)



if __name__ == "__main__" :
    und_test = Undulator(K=1.87, E=1.3e9, lambda_u=0.035, L=0.035 * 12, I=1.0)
    traj_test = TrajectoryFactory(Nb_pts=101, method=TRAJECTORY_METHOD_ODE)
    rad_test=RadiationFactory(method=RADIATION_METHOD_FARFIELD,omega=und_test.omega1())
    distance=100
    X = np.linspace(0.0, distance*1.01e-3, 101)
    Y = np.linspace(0.0,distance*1.01e-3, 101)
    sim_test = create_simulation(undulator=und_test, trajectory_fact=traj_test, distance=distance, X=X, Y=Y)

    # print('ok fin construction simulation')
    # sim_test.trajectory.draw()
    # sim_test.radiation.draw()
    # D_limite=und_test.D_max_plane_undulator(alpha=2)
    # print(D_limite)


    # print("calcul du spectre")
    # omega1=sim_test.undulator.omega1()
    # omega_array = np.arange(omega1 * 0.9, 5.0 * omega1 * 1.01, omega1*0.01)
    # sim_test.spectre_max(omega_array=omega_array)

    print("calcul du temps")
    Nb_period=und_test.L/und_test.lambda_u
    Nb_pts_trajectory=np.linspace(Nb_period*2+1 ,Nb_period*101,100)
    Nb_pts_radiation = np.linspace(10, 101,100)
    calc_time=sim_test.time_radiation(Nb_pts_trajectory, Nb_pts_radiation)
    X, Y = np.meshgrid(Nb_pts_trajectory,Nb_pts_radiation)
    fig = plt.figure()
    ax = Axes3D(fig)
    ax.plot_surface(X, Y,calc_time, rstride=1, cstride=1)
    ax.set_xlabel("trajectory pts")
    ax.set_ylabel('radiation pts')
    ax.set_zlabel("time")
    plt.show()



