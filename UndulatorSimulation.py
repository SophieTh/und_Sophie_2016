import numpy as np
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
    def __init__(self,undulator,trajectory_fact,radiation_fact=None,
                 X=None, Y=None, distance=None):
        self.undulator=undulator
        self.trajectory_fact=trajectory_fact
        self.trajectory=self.trajectory_fact.create_for_plane_undulator_ideal(undulator)
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
                                   radiation_fact=self.radiation_fact.copy())

    def change_distance(self,D):
        self.radiation.distance=D
        #update intensity
        self.radiation.intensity = self.radiation_fact.calculate_radiation_intensity(trajectory=self.trajectory,
                                                                                    distance=self.radiation.distance,
                                                                                    X_arrays=self.radiation.X,
                                                                                    Y_arrays=self.radiation.Y)
    def change_radiation_method(self,method):
        self.radiation.method=method
        #update intensity
        self.radiation.intensity = self.radiation_fact.calculate_radiation_intensity(trajectory=self.trajectory,
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




