import numpy as np
from abc import abstractmethod
import matplotlib.pyplot as plt



class Trajectory(object):
    def __init__(self, t, x, y, z, v_x, v_y, v_z, a_x, a_y, a_z):
        self.t = t
        self.x = x
        self.y = y
        self.z = z
        self.v_x = v_x
        self.v_y = v_y
        self.v_z = v_z
        self.a_x = a_x
        self.a_y = a_y
        self.a_z = a_z

    @abstractmethod
    def nb_points(self):
        raise Exception

    def convert(self):
        return

    @abstractmethod
    def multiply_by(self,cst):
        self.x *= cst
        self.y *= cst
        self.z *= cst
        self.v_x *= cst
        self.v_y *= cst
        self.v_z *= cst
        self.a_x *= cst
        self.a_y *= cst
        self.a_z *= cst

    @abstractmethod
    def error(self,trajec_test):
        return

    @abstractmethod
    def error_max(self,trajec_test):
        return

    @abstractmethod
    def error_rel(self,trajec_test):
        return

    @abstractmethod
    def copy(self):
        return Trajectory(t=self.t.copy(), x=self.x.copy(), y=self.y.copy(), z=self.z.copy(), v_x=self.v_x.copy(),
                          v_y=self.v_y.copy(), v_z=self.v_z.copy(), a_x=self.a_x.copy(),
                          a_y=self.a_y.copy(), a_z=self.a_z.copy())

    # draw all coordinate of the trajectory in function of the time
    @abstractmethod
    def plot_3D(self):
        pass

    @abstractmethod
    def plot(self):
        pass

    # draw all coordinate of the 2 trajectories in function of the time
    # It must have the same shape and the same fistr vector wich represent time.
    @abstractmethod
    def plot_2_trajectory(self,traj2):
        pass

if __name__ == "__main__" :
    t=np.linspace(0.0,20.0,1001)
    x = np.sin(t)
    y= 0.0*t
    z=np.cos(t)
    v_x=np.cos(t)
    v_y=0.0*t
    v_z=-np.sin(t)
    a_x=-np.sin(t)
    a_y=0.0*t
    a_z=-np.cos(t)
    Traj=Trajectory(t,x,y,z,v_x,v_y,v_z,a_x,a_y,a_z)
    Traj.plot()

    Traj2=Traj.copy()
    Traj2.multiply_by(3)

    Traj.plot_2_trajectory(Traj2)