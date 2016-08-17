import numpy as np
import matplotlib.pyplot as plt
from  scipy.constants import physical_constants
from mpl_toolkits.mplot3d import Axes3D
import scipy.constants as codata
import time
from scipy.interpolate import interp1d
from pySRU.Radiation import Radiation,RADIATION_LIST,RADIATION_GRID
from pySRU.MagneticField import MagneticField
from SourceUndulatorPlane import SourceUndulatorPlane
from SourceBendingmagnet import SourceBendingMagnet
from pySRU.MagneticStructureUndulatorPlane import MagneticStructureUndulatorPlane as Undulator
from pySRU.MagneticStructureBendingMagnet import MagneticStructureBendingMagnet as BM

from pySRU.TrajectoryFactory import TrajectoryFactory, TRAJECTORY_METHOD_ANALYTIC,TRAJECTORY_METHOD_ODE,\
                                       TRAJECTORY_METHOD_INTEGRATION
from pySRU.RadiationFactory import RadiationFactory,RADIATION_METHOD_NEAR_FIELD, \
                                RADIATION_METHOD_FARFIELD, RADIATION_METHOD_APPROX_FARFIELD


eV_to_J=1.602176487e-19

class Simulation(object):

    def __init__(self, source, trajectory_fact, radiation_fact, trajectory, radiation):
        self.source=source
        self.trajectory_fact=trajectory_fact
        self.radiation_fact=radiation_fact
        self.trajectory=trajectory
        self.radiation=radiation

# faire un autre constructeur pour la copie


    def copy(self):
        return Simulation(source=self.source.copy(), trajectory_fact=self.trajectory_fact.copy(),
                                   radiation_fact=self.radiation_fact.copy(),
                                   trajectory=self.trajectory.copy(), radiation=self.radiation.copy())

    #change

    def change_distance(self,D):
        self.radiation.distance=D
        #update intensity
        self.radiation.intensity = self.radiation_fact.calculate_radiation_intensity(trajectory=self.trajectory,
                                                                                     source=self.source,
                                                                                     distance=self.radiation.distance,
                                                                                     X_array=self.radiation.X,
                                                                                     Y_array=self.radiation.Y)

    def change_radiation_method(self,method):
        self.radiation_fact.method=method
        #update intensity
        self.radiation.intensity= self.radiation_fact.calculate_radiation_intensity(trajectory=self.trajectory,
                                                                                    source=self.source,
                                                                                    distance=self.radiation.distance,
                                                                                    X_array=self.radiation.X,
                                                                                    Y_array=self.radiation.Y)

    def change_trajectory_method(self,method) :
        self.trajectory_fact.method = method
        self.trajectory =self.trajectory_fact.create_from_source(source=self.source)
        self.radiation.intensity = self.radiation_fact.calculate_radiation_intensity(trajectory=self.trajectory,
                                                                                     source=self.source,
                                                                                     distance=self.radiation.distance,
                                                                                     X_array=self.radiation.X,
                                                                                     Y_array=self.radiation.Y)

    def change_initial_condition(self, initial_cond):
        self.trajectory_fact.initial_condition= initial_cond
        self.magnetic_field.change_Zo(Zo=initial_cond[5])
        self.trajectory = self.trajectory_fact.create_for_parameter(parameter=self.parameter,
                                                                    B=self.magnetic_field)

        self.radiation.intensity = self.radiation_fact.calculate_radiation_intensity(trajectory=self.trajectory,
                                                                                     source=self.parameter,
                                                                                     distance=self.radiation.distance,
                                                                                     X_arrays=self.radiation.X,
                                                                                     Y_arrays=self.radiation.Y)

    def change_omega(self, omega) :
        self.radiation_fact.omega=omega
        self.radiation.intensity = self.radiation_fact.calculate_radiation_intensity(trajectory=self.trajectory,
                                                                                     source=self.source,
                                                                                     distance=self.radiation.distance,
                                                                                     X_array=self.radiation.X,
                                                                                     Y_array=self.radiation.Y)

    def change_harmonic_number(self, harmonic_number) :
        self.radiation_fact.omega=harmonic_number*self.source.omega1()
        self.radiation.intensity = self.radiation_fact.calculate_radiation_intensity(trajectory=self.trajectory,
                                                                                     source=self.source,
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
        self.trajectory = self.trajectory_fact.create_from_source(source=self.source)

    def change_Nb_pts_trajectory(self, Nb_pts):
        self.trajectory_fact.Nb_pts = int(Nb_pts)
        self.trajectory = self.trajectory_fact.create_from_source(source=self.source)
        self.radiation.intensity = self.radiation_fact.calculate_radiation_intensity(trajectory=self.trajectory,
                                                                                     source=self.source,
                                                                                     distance=self.radiation.distance,
                                                                                     X_array=self.radiation.X,
                                                                                     Y_array=self.radiation.Y)

    def change_nb_period(self, Nb_period):
        self.source.magnetic_structure.length=Nb_period*self.source.magnetic_structure.period_length
        self.source.magnetic_field=self.source.magnetic_structure.create_magnetic_field()
        self.trajectory_fact.Nb_pts = int(self.source.choose_nb_pts_trajectory(2))
        self.trajectory = self.trajectory_fact.create_from_source(source=self.source)
        self.radiation.intensity = self.radiation_fact.calculate_radiation_intensity(trajectory=self.trajectory,
                                                                                     source=self.source,
                                                                                     distance=self.radiation.distance,
                                                                                     X_array=self.radiation.X,
                                                                                     Y_array=self.radiation.Y)

    def change_Nb_pts_radiation(self, Nb_pts):
        self.radiation_fact.Nb_pts = Nb_pts
        self.radiation.change_Nb_pts(Nb_pts)
        self.radiation.intensity = self.radiation_fact.calculate_radiation_intensity(trajectory=self.trajectory,
                                                                                     source=self.source,
                                                                                     distance=self.radiation.distance,
                                                                                     X_array=self.radiation.X,
                                                                                     Y_array=self.radiation.Y)


    def change_XY_radiation(self,X,Y):

        intensity = self.radiation_fact.calculate_radiation_intensity(trajectory=self.trajectory,
                                                                                     source=self.source,
                                                                                     distance=self.radiation.distance,
                                                                                     X_array=X,
                                                                                     Y_array=Y)
        self.radiation.X=X
        self.radiation.Y=Y
        self.radiation.intensity=intensity



    #TODO possible que pour Undulator, a revoir , a tester
    def calculate_until_wave_number(self,wave_number):
        harmonic_number=np.floor(self.radiation_fact.omega/self.source.harmonic_frequency(1))
        self.radiation_fact.omega=self.source.harmonic_frequency(harmonic_number)
        X=np.array([])
        Y = np.array([])
        t = np.linspace(0.0, 2.0* np.pi, self.radiation_fact.Nb_pts)
        for i in range(int(wave_number)+1) :
            Xi,Yi= self.source.describe_wave(distance=self.radiation.distance,harmonic_number=harmonic_number,
                                              wave_number=i,t=t)
            X=np.concatenate((X,Xi))
            Y=np.concatenate((Y,Yi))
        self.change_XY_radiation(X=X,Y=Y)

    def calculate_on_central_cone(self):
        theta_max=self.source.angle_deflection_central_cone()
        if self.radiation.distance==None :
            X_max=theta_max
            Y_max=theta_max
        else :
            X_max=self.radiation.distance*theta_max
            Y_max=self.radiation.distance*theta_max

        X=np.linspace(0.0,X_max,self.radiation_fact.Nb_pts)
        Y=np.linspace(0.0,Y_max,self.radiation_fact.Nb_pts)
        X,Y=np.meshgrid(X,Y)
        self.change_XY_radiation(X=X,Y=Y)



    def calculate_for_observation_angles(self,observation_angle):
        D=self.radiation.distance
        if D==None :
            D=1
        X = np.array([])
        Y = np.array([])
        t = np.linspace(0.0, 2.0 * np.pi, self.radiation_fact.Nb_pts)
        for theta in observation_angle:
            if theta==0.0 :
                Xi=np.array([0.0])
                Yi=np.array([0.0])
            else :
                Xi=np.cos(t)*D*theta
                Yi=np.sin(t)*D*theta
            X = np.concatenate((X, Xi))
            Y = np.concatenate((Y, Yi))
        #X,Y=np.meshgrid(X,Y)
        self.change_XY_radiation(X=X, Y=Y)


            # spectre




    def spectre(self, omega_array=None):
        if omega_array == None:
            omega1 = self.source. choose_photon_frequency()
            omega_array1 = np.arange(omega1 * 0.9, omega1 * 1.11, omega1 * 0.01)
            omega_array3 = np.arange(omega1 * 2.9, omega1 * 3.11, omega1 * 0.01)
            omega_array2=np.arange(omega1*1.2,omega1*2.85,omega1*0.2)
            print(len(omega_array1.shape))
            print(len(omega_array2.shape))
            print(len(omega_array3.shape))
            # omega_array=np.concatenate((omega_array2,omega_array3))
            # omega_array=np.concatenate((omega_array1,omega_array))
            omega_array=np.concatenate((omega_array1,omega_array2,omega_array3))
        spectre = np.zeros_like(omega_array)
        print(len(spectre))
        for i in range(len(spectre)):
            print(i)
            self.change_omega(omega_array[i])
            spectre[i] = self.radiation.integration()
        return spectre, omega_array

    def spectre_on_radiation_max(self, omega_array=None):
        if omega_array == None:
            omega1 = self.source.omega1()
            omega_array = np.arange(omega1 * 0.9, 5.0 * omega1 * 1.1, omega1 * 0.01)
        spectre = np.zeros_like(omega_array)
        print(len(spectre))
        for i in range(len(spectre)):
            print(i)
            self.change_omega(omega_array[i])
            spectre[i] = self.radiation.max()
        return spectre, omega_array

    def spectre_on_axis(self, omega_array=None):
        save_X=self.radiation.X
        X=np.array([0.0])
        save_Y=self.radiation.Y
        Y=np.array([0.0])
        self.change_XY_radiation(X=X,Y=Y)
        if omega_array == None:
            omega1 = self.source.harmonic_frequency(1)
            omega_array = np.arange(omega1 * 0.9, 5.0 * omega1 * 1.1, omega1 * 0.001)
        spectre = np.zeros_like(omega_array)
        for i in range(len(spectre)):
            self.change_omega(omega_array[i])
            spectre[i] = self.radiation.intensity.max()
        self.change_XY_radiation(X=save_X,Y=save_Y)
        return spectre, omega_array


    def plot_spectre(self,omega_array=None):
        import matplotlib.pyplot as plt
        spectre,omega_array=self.spectre(omega_array=omega_array)
        plt.plot(omega_array,spectre)
        plt.show()

    def plot_spectre_central_cone(self, omega_array=None):
        self.radiation_fact.distance=None
        self.calculate_on_central_cone()
        self.plot_spectre(omega_array=omega_array)

    def plot_spectre_on_axis(self, omega_array=None):
        import matplotlib.pyplot as plt
        spectre,omega_array=self.spectre_on_axis(omega_array=omega_array)
        omega1 = self.source.harmonic_frequency(1)
        harm_num_min = np.ceil(omega_array.min() / omega1)
        harm_num_max = np.floor(omega_array.max() / omega1)
        harm_number=np.arange(harm_num_min,harm_num_max+0.5,1)
        spectre_theoritical=np.zeros_like(harm_number)
        print(len(spectre_theoritical))
        for i in range(len(spectre_theoritical)):
            spectre_theoritical[i]=self.source.theorical_flux_on_axis(harm_number[i])
        print(spectre_theoritical)
        plt.plot(omega_array,spectre,label='calculad radiation on axis')
        plt.plot(harm_number*omega1,spectre_theoritical, "g*",label='theoritical value')
        plt.legend()
        plt.show()

    # plot

    def plot_magnetic_field_along_Z(self,X=None,Y=None):
        Z=self.trajectory.t*self.source.electron_speed()*codata.c
        if Y==None :
            Y=self.trajectory_fact.initial_condition[4]
        if X == None:
            X = self.trajectory_fact.initial_condition[3]
        self.source.magnetic_field.plot_z(X=X,Y=Y,Z=Z)

    def print_parameters(self):
        self.source.print_parameters()

        self.trajectory_fact.print_parameters()

        self.radiation_fact.print_parameters()

    def plot_everything(self):
        # self.trajectory.plot()
        # self.trajectory.plot_3D()
        self.radiation.plot()
        self.calculate_on_central_cone()
        self.radiation.plot()
        # self.plot_spectre()
        # self.plot_spectre_on_axis()


    # error

    def error_radiation_method_nb_period(self, method, nb_period):
        sim2=self.copy()
        sim2.change_radiation_method(method)
        error=np.zeros_like(nb_period)
        print(len(nb_period))
        for i in range(len(nb_period)) :
            print(i)
            self.change_nb_period(nb_period[i])
            print(self.trajectory_fact.Nb_pts)
            sim2.change_nb_period(nb_period[i])
            error[i]=self.radiation.error_max(sim2.radiation)
        return error

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
            print(i)
            self.change_Nb_pts_trajectory(nb_pts[i])
            sim2.change_Nb_pts_trajectory(nb_pts[i])
            error_traj[:,i]=self.trajectory.error_max(sim2.trajectory)
            error_rad[i]=self.radiation.error_max(sim2.radiation)
        traj_error_traj=self.trajectory_fact.create_from_array(error_traj)
        return error_rad, traj_error_traj

    def relativ_error_radiation_method_distance(self,method,D):
        sim2=self.copy()
        sim2.change_radiation_method(method)
        error_relativ=np.array(len(D))
        for i in range(len(D)) :
            self.change_distance(D[i])
            sim2.change_distance(D[i])
            error_relativ[i]=self.radiation.relativ(sim2.radiation)
        return error_relativ

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



    #TODO a tester attention possible que pour BENDING magnet
    def create_theoric_radiation(self):
        #TODO faire le cas None
        radiation=self.radiation.copy()
        shape1 = radiation.intensity.shape
        # X = radiation.X.flatten()
        # Y = radiation.X.flatten()
        # res=radiation.intensity.flatten()
        # cas grid :
        print(type(radiation))
        print(type(radiation.Y))
        Y=radiation.Y
        X=radiation.X
        for i in range(Y.shape[0]):
            for j in range(Y.shape[1]) :
                observation_angle=(Y[i,j]/radiation.distance)
                #TODO changer pour qu'on calcul avec des angles
                radiation.intensity[i,j]=self.source.radiation_theoric(
                    omega=self.radiation_fact.omega,observation_angle=observation_angle)
        return radiation



    def create_theoric_radiation2(self):
        radiation = self.radiation.copy()
        shape1 = radiation.intensity.shape
        # X = radiation.X.flatten()
        # Y = radiation.X.flatten()
        # res=radiation.intensity.flatten()
        # cas grid :
        print(type(radiation))
        print(type(radiation.Y))
        Y = radiation.Y
        X = radiation.X
        for i in range(Y.shape[0]):
            for j in range(Y.shape[1]):
                observation_angle = (Y[i, j] / radiation.distance)
                radiation.intensity[i, j] = self.source.radiation_theoric(
                    omega=self.radiation_fact.omega, observation_angle=observation_angle)
        return radiation

#

def create_simulation(magnetic_structure,electron_beam, magnetic_field=None, photon_energy=None,
                      traj_method=TRAJECTORY_METHOD_ANALYTIC,Nb_pts_trajectory=None,
                      rad_method=RADIATION_METHOD_APPROX_FARFIELD,  formula=1,
                      initial_condition=None, distance=None,XY_are_list=False,X=None,Y=None) :

    if type(magnetic_structure)==Undulator :
        source=SourceUndulatorPlane(undulator=magnetic_structure,
                                    electron_beam=electron_beam, magnetic_field=magnetic_field)
        print("source undulateur")
    elif type(magnetic_structure)==BM:
        source = SourceBendingMagnet(magnetic_structure=magnetic_structure,
                                      electron_beam=electron_beam, magnetic_field=magnetic_field)
        print("source BM")
    else :
        raise Exception('magnet type unknown')

    if photon_energy==None :
        omega=source.choose_photon_frequency()
    else :
        omega = photon_energy *eV_to_J / codata.hbar

    # print('omega=')
    # print(omega)

    if Nb_pts_trajectory==None :
        Nb_pts_trajectory = int(source.choose_nb_pts_trajectory(2))

    if distance==None and (rad_method==RADIATION_METHOD_NEAR_FIELD or rad_method==RADIATION_METHOD_FARFIELD) :
        distance=source.choose_distance_automatic(2)

    if X ==None or Y== None :
        if (X != None) :
            Y=X
        elif Y != None :
            X=Y
        else :
            if not (XY_are_list):
                theta_max=source.choose_angle_deflection_max()
                if distance==None :
                    X_max=theta_max
                    Y_max=theta_max
                else :
                    X_max = distance * theta_max
                    Y_max = distance * theta_max
                X = np.linspace(0.0, X_max, 101)
                Y = np.linspace(0.0, Y_max, 101)
            else :
                X=np.array([0.0])
                Y=np.array([0.0])

    if type(X) == float:
        X= np.linspace(0.0, X, 101)
    if type(Y) == float:
        Y = np.linspace(0.0, Y, 101)




    if X.shape != Y.shape :
        raise Exception('X and Y must have the same shape')
    Nb_pts_radiation=len(X.flatten())

    #print('step 1')
    traj_fact=TrajectoryFactory(Nb_pts=Nb_pts_trajectory,method=traj_method,initial_condition=initial_condition)
    if (traj_fact.initial_condition == None):
        # print('crearte initial cond automat')
        traj_fact.initial_condition = source.choose_initial_contidion_automatic()
        # print(traj_fact.initial_condition)


    #print('step 2')
    rad_fact=RadiationFactory(method=rad_method,omega=omega,Nb_pts=Nb_pts_radiation,formula=formula)

    #print('step 3')
    trajectory=traj_fact.create_from_source(source=source)
    #print('step 4')
    radiation = rad_fact.create_for_one_relativistic_electron(trajectory=trajectory, source=source, XY_are_list=XY_are_list,
                                                              distance=distance, X=X, Y=Y)

    #print('step 5')
    return Simulation(source=source, trajectory_fact=traj_fact,
                               radiation_fact=rad_fact, trajectory=trajectory, radiation=radiation)







# exemple

def Exemple_minimum():

    from  ElectronBeam import ElectronBeam

    undulator_test = Undulator(K=1.87, period_length=0.035, length=0.035 * 14)
    electron_beam_test = ElectronBeam(Electron_energy=1.3, I_current=1.0)

    simulation_test = create_simulation(magnetic_structure=undulator_test,electron_beam=electron_beam_test)

    simulation_test.print_parameters()

    simulation_test.trajectory.plot_3D()

    simulation_test.radiation.plot()

    simulation_test.change_distance(D=100)
    simulation_test.radiation.plot()

    simulation_test.calculate_on_central_cone()
    simulation_test.radiation.plot()


def Exemple_meshgrid():
    from ElectronBeam import ElectronBeam

    beam_test = ElectronBeam(Electron_energy=1.3, I_current=1.0)
    beam_ESRF = ElectronBeam(Electron_energy=6.0, I_current=0.2)
    und_test = Undulator(K=1.87, period_length=0.035, length=0.035 * 14)
    ESRF18 = Undulator(K=1.68, period_length=0.018, length=2.0)


    print('create simulation for a given screen')
    X=np.linspace(-0.02,0.02,100)
    Y=np.linspace(-0.02,0.02,100)
    simulation_test = create_simulation(magnetic_structure=ESRF18,electron_beam=beam_ESRF,
                                        traj_method=TRAJECTORY_METHOD_ANALYTIC,rad_method=RADIATION_METHOD_FARFIELD,
                                        distance= 100,X=X,Y=Y,photon_energy=7876.0)

    simulation_test.trajectory.plot_3D()
    simulation_test.radiation.plot()

    print('create simulation for a maximal X and Y given')
    simulation_test = create_simulation(magnetic_structure=ESRF18, electron_beam=beam_ESRF,
                                        traj_method=TRAJECTORY_METHOD_ANALYTIC, rad_method=RADIATION_METHOD_FARFIELD,
                                        distance=100,X=0.02,Y=0.02)

    simulation_test.radiation.plot()

    simulation_test.plot_spectre_on_axis()

    simulation_test.plot_spectre_central_cone()

    simulation_test.plot_spectre()






# TODO ne marche plus  qd on rentre a la main l'energy du photon?
def Exemple_list():
    from pySRU.ElectronBeam import ElectronBeam

    beam_ESRF = ElectronBeam(Electron_energy=6.0, I_current=0.2)
    ESRF18 = Undulator(K=1.68, period_length=0.018, length=2.0)

    #photon_energy = 7876.0,

    X = np.linspace(0.0,0.0002,1001)
    Y = np.linspace(0.0, 0.0002, 1001)
    simulation_test = create_simulation(magnetic_structure=ESRF18,electron_beam=beam_ESRF,
                                        X=X,Y=Y,XY_are_list=True)

    simulation_test.print_parameters()
    simulation_test.trajectory.plot_3D()

    simulation_test.radiation.plot()

    simulation_test.calculate_until_wave_number(2)


    simulation_test.radiation.plot()

    simulation_test.radiation.plot_wave(Nb_pts=simulation_test.radiation_fact.Nb_pts)

    observation_angle = np.linspace(0.0, simulation_test.source.angle_wave_number(1, 2), 51)
    simulation_test.calculate_for_observation_angles(observation_angle=observation_angle)
    simulation_test.radiation.plot()





if __name__ == "__main__" :


    Exemple_minimum()
    #Exemple_meshgrid()
    #Exemple_list()