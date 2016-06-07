import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as codata
from Trajectory import Trajectory
from Radiation import Radiation
from MagneticField import MagneticField
from UndulatorParameter import UndulatorParameters as Undulator
from TrajectoryFactory import TrajectoryFactory, TRAJECTORY_METHOD_ANALYTIC,TRAJECTORY_METHOD_ODE,\
                                        TRAJECTORY_METHOD_INTEGRATION
from RadiationFactory import RadiationFactory ,RADIATION_METHOD_NEAR_FIELD, \
                                RADIATION_METHOD_FARFIELD, RADIATION_METHOD_APPROX_FARFIELD


class UndulatorSimulation(object):
    def __init__(self,undulator,trajectory_fact,radiation_fact=None,magnetic_field=None,
                 X=None, Y=None, distance=None):
        self.undulator=undulator
        self.magnetic_filed=magnetic_field
        self.trajectory_fact=trajectory_fact
        self.trajectory=self.trajectory_fact.create_for_plane_undulator(undulator=self.undulator,
                                                                              Z_By=self.magnetic_filed)
        if (radiation_fact == None):
            self.radiation_fact=RadiationFactory( method=RADIATION_METHOD_APPROX_FARFIELD,
                                                          omega=self.undulator.omega1())
        else:
            self.radiation_fact=radiation_fact

        self.radiation=self.radiation_fact.create_for_single_electron(trajectory=self.trajectory,
                                                                          undulator=self.undulator,
                                                                          distance=distance,
                                                                          X=X,Y=Y)


    def copy(self):
        return UndulatorSimulation(undulator=self.undulator.copy(),trajectory_fact=self.trajectory_fact.copy(),
                                   magnetic_field=self.magnetic_filed.copy(),radiation_fact=self.radiation_fact.copy())

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
        self.radiation.intensity = self.radiation_fact.calculate_radiation_intensity(trajectory=self.trajectory,
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
                                                                                     undulator=self.undulator,
                                                                                 distance=self.radiation.distance,
                                                                                 X_arrays=self.radiation.X,
                                                                                 Y_arrays=self.radiation.Y)

    def error_radiation_method(self,method,D):
        sim2=self.copy()
        sim2.change_radiation_method(method)
        error=np.array(len(D))
        for i in range(len(D)) :
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
        for i in range(len(spectre)) :
            print(i)
            self.change_omega(omega_array[i])
            spectre[i]=self.radiation.integration()
        return spectre , omega_array

    def spectre_max(self, omega_array=None):
        spectre,omega_array = self.spectre(omega_array=omega_array)
        omega1=self.undulator.omega1()
        omega1_min=np.ceil(omega_array[0]/omega1)
        print("harmonic min =")
        print(omega1_min)
        omega1_max=np.floor(omega_array[-1]/omega1)
        print("harmonic max =")
        print(omega1_max)
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

    def D_max_for_n_cst(self,alpha):
        lim=self.undulator.Beta()*codata.c*np.abs(self.trajectory.t[-1])+(self.undulator.K*self.undulator.K*self.undulator.lambda_u)/(8.0*self.undulator.gamma()*self.undulator.gamma()*2.0*np.pi)
        return lim *10**alpha


if __name__ == "__main__" :
    und_test = Undulator(K=1.87, E=1.3e9, lambda_u=0.035, L=0.035 * 12, I=1.0)
    traj_test = TrajectoryFactory(Nb_pts=201, method=TRAJECTORY_METHOD_ANALYTIC)

    sim_test = UndulatorSimulation(undulator=und_test, trajectory_fact=traj_test)


    # sim_test.trajectory.draw()
    #sim_test.radiation.draw()
    D_limite=sim_test.D_max_for_n_cst(alpha=2)
    print(D_limite)
    # print("calcul du spectre")
    # omega1=sim_test.undulator.omega1()
    # omega_array=np.arange(omega1*(1.0-1e-3),omega1*(1.0+1e-5),omega1*0.1e-4)
    # sim_test.spectre_max(omega_array=omega_array)

