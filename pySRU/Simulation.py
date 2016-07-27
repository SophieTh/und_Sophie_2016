import numpy as np
import matplotlib.pyplot as plt
from  scipy.constants import physical_constants
from mpl_toolkits.mplot3d import Axes3D
import scipy.constants as codata
import time
from scipy.interpolate import interp1d
from pySRU.RadiationList import RadiationList
from pySRU.RadiationGrid import RadiationGrid
from pySRU.MagneticField import MagneticField
from Source import Source
from pySRU.MagneticStructureUndulatorPlane import MagneticStructureUndulatorPlane as Undulator ,PLANE_UNDULATOR,BENDING_MAGNET

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
        self.trajectory_fact.Nb_pts = Nb_pts
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
        self.radiation.X=X
        self.radiation.Y=Y
        self.radiation.intensity = self.radiation_fact.calculate_radiation_intensity(trajectory=self.trajectory,
                                                                                     source=self.source,
                                                                                     distance=self.radiation.distance,
                                                                                     X_array=self.radiation.X,
                                                                                     Y_array=self.radiation.Y)



    #TODO changer et mettre possible que pour undulateur ?
    def calculate_on_wave(self, wave_number,t):
        if wave_number==0 :
            X=np.array([0.0])
            Y=X
        else :
            n=np.floor(self.source.harmonic_number(self.radiation_fact.omega))
            theta=self.source.theta(n=n, l=wave_number)
            R=self.radiation.distance*np.tan(theta)
            X=R*np.cos(t)
            Y=R*np.sin(t)
        return X,Y


    def calculate_until_wave_number(self,wave_number,harmonic_number=1.):
        self.radiation_fact.omega=self.source.omega1()*harmonic_number
        X=np.array([])
        Y = np.array([])
        t = np.linspace(0.0, 0.5 * np.pi, self.radiation_fact.Nb_pts)
        for i in range(int(wave_number)+1) :
            Xi,Yi= self.calculate_on_wave(i,t)
            X=np.concatenate((X,Xi))
            Y=np.concatenate((Y,Yi))
        self.change_XY_radiation(X=X,Y=Y)




            # spectre

# TODO spectre sur le cone , comparaison avec les formules)
    def spectre(self, omega_array=None):
        if omega_array == None:
            omega1 = self.source.omega1()
            omega_array = np.arange(omega1 * 0.9, 3.0 * omega1 * 1.05, omega1 * 0.01)
        spectre = np.zeros_like(omega_array)
        #print(len(spectre))
        for i in range(len(spectre)):
            #print(i)
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
            omega1 = self.source.omega1()
            omega_array = np.arange(omega1 * 0.9, 5.0 * omega1 * 1.1, omega1 * 0.01)
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
        xy_max=self.source.theta(n=1,l=1)/2.
        X = np.linspace(0.0,xy_max,self.radiation_fact.Nb_pts)
        Y = np.linspace(0.0, xy_max, self.radiation_fact.Nb_pts)
        X , Y = np.meshgrid(X,Y)
        if type(self.radiation)==RadiationList :
            self.radiation=RadiationGrid(intensity=0.0,X=0.0,Y=0.0)
        self.change_XY_radiation(X=X,Y=Y)
        self.plot_spectre(omega_array=omega_array)

    def plot_spectre_on_axis(self, omega_array=None):
        import matplotlib.pyplot as plt
        spectre,omega_array=self.spectre_on_axis(omega_array=omega_array)
        omega1 = self.source.omega1()
        harm_num_min = np.ceil(omega_array[0] / omega1)
        harm_num_max = np.floor(omega_array[-1] / omega1)
        harm_number=np.arange(harm_num_min,harm_num_max+0.5,1)
        spectre_theoritical=np.zeros_like(harm_number)
        #TODO changer doit y avoir un autre moyen ?
        for i in range(len(spectre_theoritical)):
            spectre_theoritical[i]=self.source.flux_on_axis_theoric(harm_number[i])
        plt.plot(omega_array,spectre,label='calculad radiation on axis')
        plt.plot(harm_number*omega1,spectre_theoritical, "g*",label='theoritical value')
        plt.legend()
        plt.show()


        # ???? TODO what is that ?

    def spectre_max(self, omega_array=None):
        start_time = time.time()
        spectre, omega_array = self.spectre3(omega_array=omega_array)
        interval = time.time() - start_time
        # print("interval temps :")
        # print(interval)
        omega1 = self.parameter.omega1()
        omega1_min = np.ceil(omega_array[0] / omega1)
        # print("harmonic min =")
        # print(omega1_min)
        omega1_max = np.floor(omega_array[-1] / omega1)
        # print("harmonic max =")
        # print(omega1_max)
        omega1_array = omega1 * np.ones_like(omega_array)
        plt.plot(omega_array, spectre)
        i = omega1_min
        while i <= omega1_max:
            plt.plot(i * omega1_array, spectre, color="red")
            i += 1
        print("omega1 =")
        print(omega1)
        omega_max = omega_array[np.argmax(spectre)]
        print("omega max =")
        print(omega_max)
        plt.show()
        return spectre, omega_array, omega_max

    # plot

    def plot_magnetic_field_along_Z(self,X=None,Y=None):
        Z=self.trajectory.t*self.source.Beta_et()*codata.c
        if Y==None :
            Y=self.trajectory_fact.initial_condition[4]
        if X == None:
            X = self.trajectory_fact.initial_condition[3]
        self.source.magnetic_field.plot_z(X=X,Y=Y,Z=Z)

    def print_parameters(self):
        self.source.print_parameters()

        self.trajectory_fact.print_parameters()

        self.radiation_fact.print_parameters()
        print('    harmonic number : %.3f' %(self.radiation_fact.omega/self.source.omega1()))


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
    # def error_trajectory_cst_magnetic_field(self,nb_pts):
    #     Bo = (self.parameter.K / (93.4 * self.parameter.lambda_u))
    #     f = Bo * codata.e / (self.parameter.gamma() * codata.m_e)
    #     vz_0=self.trajectory_fact.initial_condition[2]/codata.c
    #     x_0=self.trajectory_fact.initial_condition[3]/codata.c
    #     error = np.zeros((10,len(nb_pts)))
    #     traj_ref=self.trajectory.copy()
    #     for i in range(len(nb_pts)):
    #         self.change_Nb_pts_trajectory_only(nb_pts[i])
    #         traj_ref_arrays= self.trajectory_fact.analytical_trajectory_cst_magnf(f=f,
    #                                                 t=self.trajectory.t,vz_0=vz_0,x_0=x_0)
    #         traj_ref.convert(traj_ref_arrays)
    #         error_max=self.trajectory.error_max(traj_ref)
    #         error[:,i]=error_max
    #         # print(error_max.shape)
    #         # for j in range(9) :
    #         #      error[j+1][i] = error_max[j+1]
    #     traj_err=self.trajectory_fact.create_from_array(error)
    #     return traj_err

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


#

def create_simulation(magnetic_structure,electron_beam, magnetic_field=None, photon_energy=None,
                      traj_method=TRAJECTORY_METHOD_ODE,Nb_pts_trajectory=None,
                      rad_method=RADIATION_METHOD_APPROX_FARFIELD,  formule=1,
                      initial_condition=None, distance=None,XY_are_list=False,X=None,Y=None) :


    source=Source(magnetic_structure=magnetic_structure,electron_beam=electron_beam,magnetic_field=magnetic_field)

    if distance==None :
        distance= source.D_min(2) * 2.0

    if photon_energy==None :
        omega=source.omega1()
    else :
        omega = photon_energy * eV_to_J / codata.hbar


    if Nb_pts_trajectory==None :
        Nb_pts_trajectory = int(np.floor(source.n_min(2))) + 1

    if X ==None or Y== None :
        if (X != None) :
            Y=X
        elif Y != None :
            X=Y
        else :
            if not (XY_are_list):
                theta_max=source.theta_max()
                X_max = distance * theta_max
                Y_max = distance * theta_max
                X = np.linspace(0.0, X_max, 101)
                Y = np.linspace(0.0, Y_max, 101)
            else :
                X=np.array([0.0])
                Y=np.array([0.0])

    if type(X) == float or type(X) == int:
        X= np.linspace(0.0, X, 101)
    if type(Y) == float or type(Y) == int:
        Y = np.linspace(0.0, Y, 101)

    if X.shape != Y.shape :
        raise Exception('X and Y must have the same shape')
    Nb_pts_radiation=len(X.flatten())

    #print('step 1')
    traj_fact=TrajectoryFactory(Nb_pts=Nb_pts_trajectory,method=traj_method,initial_condition=initial_condition)
    #print('step 2')
    rad_fact=RadiationFactory(method=rad_method,omega=omega,Nb_pts=Nb_pts_radiation,formula=1)
    #print('step 3')
    if (traj_fact.initial_condition == None):
        traj_fact.initial_condition = traj_fact.choise_initial_condition(source)

    trajectory=traj_fact.create_from_source(source=source)
    #print('step 4')
    radiation = rad_fact.create_for_single_electron(trajectory=trajectory, source=source,XY_are_list=XY_are_list,
                                                          distance=distance, X=X, Y=Y)
    #print('step 5')
    return Simulation(source=source, trajectory_fact=traj_fact,
                               radiation_fact=rad_fact, trajectory=trajectory, radiation=radiation)







# exemple

def Exemple_minimum():

    from  ElectronBeam import ElectronBeam

    undulator_test = Undulator(K=1.87, lambda_u=0.035, L=0.035 * 14)
    electron_beam_test = ElectronBeam(E=1.3e9, I=1.0)

    simulation_test = create_simulation(magnetic_structure=undulator_test,electron_beam=electron_beam_test)

    simulation_test.print_parameters()

    simulation_test.trajectory.plot_3D()

    simulation_test.radiation.plot()


def Exemple_meshgrid():
    from MagneticStructureUndulatorPlane import MagneticStructureUndulatorPlane as Undulator
    from ElectronBeam import ElectronBeam

    undulator_test = Undulator(K=1.87, lambda_u=0.035, L=0.035 * 14)
    electron_beam_test = ElectronBeam(E=1.3e9, I=1.0)


    print('create simulation for a given sreen')
    X=np.linspace(-0.1,0.1,200)
    Y=np.linspace(-0.1,0.1,200)
    simulation_test = create_simulation(magnetic_structure=undulator_test,electron_beam=electron_beam_test,
                                        traj_method=TRAJECTORY_METHOD_ANALYTIC,rad_method=RADIATION_METHOD_FARFIELD,
                                        distance= 100,X=X,Y=Y)

    simulation_test.trajectory.plot_3D()
    simulation_test.radiation.plot()

    print('create simulation for a maximal X and Y given')
    simulation_test = create_simulation(magnetic_structure=undulator_test, electron_beam=electron_beam_test,
                                        traj_method=TRAJECTORY_METHOD_ANALYTIC, rad_method=RADIATION_METHOD_FARFIELD,
                                        distance=100, X=0.1, Y=0.1)

    simulation_test.radiation.plot()

    simulation_test.plot_spectre()

    simulation_test.plot_spectre_central_cone()

    simulation_test.plot_spectre_on_axis()



#
def Exemple_list():
    from MagneticStructureUndulatorPlane import MagneticStructureUndulatorPlane as Undulator
    from ElectronBeam import ElectronBeam

    beam_ESRF = ElectronBeam(E=6.0e9, I=0.2)
    ESRF18 = Undulator(K=1.68, lambda_u=0.018, L=2.0)


    X = np.linspace(0.0,0.05,1001)
    Y = np.linspace(0.0, 0.05, 1001)
    simulation_test = create_simulation(magnetic_structure=ESRF18,electron_beam=beam_ESRF, photon_energy=7876.0,
                                        X=X,Y=Y,XY_are_list=True)

    simulation_test.trajectory.plot_3D()

    simulation_test.radiation.plot()

    simulation_test.calculate_until_wave_number(2)

    simulation_test.radiation.plot()

    simulation_test.radiation.plot_wave(Nb_pts=simulation_test.radiation_fact.Nb_pts)







if __name__ == "__main__" :


    #Exemple_minimum()
    Exemple_meshgrid()
    #Exemple_list()