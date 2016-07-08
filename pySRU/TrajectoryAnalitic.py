import numpy as np
import matplotlib.pyplot as plt
from pySRU.Trajectory import Trajectory
from pySRU.TrajectoryArray import TrajectoryArray


class TrajectoryAnalitic(Trajectory):
    def __init__(self, t, x, y, z, v_x, v_y, v_z, a_x, a_y, a_z):
        super(self.__class__, self).__init__(t,x,y,z,v_x,v_y,v_z,a_x,a_y,a_z)

    def nb_points(self):
        return len(self.t)

    def convert(self):
        trajectory_array=TrajectoryArray(self.t,
                                          self.x(self.t),self.y(self.t),self.z(self.t),
                                          self.v_x(self.t),self.v_y(self.t),self.v_z(self.t),
                                          self.a_x(self.t),self.a_y(self.t),self.a_z(self.t),
                                          )
        return trajectory_array

    def multiply_by(self,cst):
        self.x = (lambda t :cst*self.x(t) )
        self.y =  (lambda t :cst*self.y(t) )
        self.z = (lambda t :cst*self.z(t) )
        self.v_x = (lambda t :cst*self.v_x(t) )
        self.v_y = (lambda t :cst*self.v_y(t) )
        self.v_z = (lambda t :cst*self.v_z(t) )
        self.a_x = (lambda t :cst*self.a_x(t) )
        self.a_y = (lambda t :cst*self.a_y(t) )
        self.a_z = (lambda t :cst*self.a_z(t) )

    def error(self,trajec_test):
        if type(trajec_test)==TrajectoryAnalitic :
            error=self.error_analitic(trajec_test)
        else :
            error=self.error_array(trajec_test)
        return error

    def error_analitic(self,trajec_test):
        # faire erreur sur la difference plutot!
        if all((self.t- trajec_test.t) <= np.abs(self.t).max()*1e-6) :
            raise Exception("Problem : the two trajectory have not the same vector t ??")
        error = self.copy()
        error.x = (lambda t : np.abs(self.x(t) - trajec_test.x(t)))
        error.y = (lambda t : np.abs(self.y(t) - trajec_test.y(t)))
        error.z = (lambda t : np.abs(self.z(t) - trajec_test.z(t)))
        error.v_x = (lambda t : np.abs(self.v_x(t) - trajec_test.v_x(t)))
        error.v_y = (lambda t : np.abs(self.v_y(t) - trajec_test.v_y(t)))
        error.v_z = (lambda t : np.abs(self.v_z(t) - trajec_test.v_z(t)))
        error.a_x= (lambda t : np.abs(self.a_x(t) - trajec_test.a_x(t)))
        error.a_y = (lambda t : np.abs(self.a_y(t) - trajec_test.a_y(t)))
        error.a_z = (lambda t : np.abs(self.a_z(t) - trajec_test.a_z(t)))
        return error
        # return Trajectory(error[0],error[1],error[2],error[3],error[4],error[5],error[6]
        #                   ,error[7],error[8],error[9])

    def error_array(self, trajec_test):
        # faire erreur sur la difference plutot!
        if all((self.t - trajec_test.t) <= np.abs(self.t).max() * 1e-6):
            raise Exception("Problem : the two trajectory have not the same vector t ??")
        if type(trajec_test)==TrajectoryAnalitic :
            raise Exception("trajec_test must be a Trajectory Array")
        traj_array=self.convert_arrays()
        error = traj_array.error(trajec_test)
        return error


    def error_max(self,trajec_test):
        if all((self.t - trajec_test.t) <= np.abs(self.t).max() * 1e-6):
            raise Exception("Problem : the two trajectory have not the same vector t ??")
        if type(trajec_test)==TrajectoryAnalitic :
            raise Exception("trajec_test must be a Trajectory Array")
        traj_array=self.convert_arrays()
        error = traj_array.error(trajec_test)
        return error

    def error_rel(self,trajec_test):
        if all((self.t - trajec_test.t) <= np.abs(self.t).max() * 1e-6):
            raise Exception("Problem : the two trajectory have not the same vector t ??")
        if type(trajec_test)==TrajectoryAnalitic :
            raise Exception("trajec_test must be a Trajectory Array")
        traj_array=self.convert_arrays()
        error = traj_array.error(trajec_test)
        return error

    # draw all coordinate of the trajectory in function of the time
    def plot_3D(self):
        z=self.z(self.t)
        y = self.y(self.t)
        x = self.x(self.t)
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        ax.plot(z,x,y, label='trajectory')
        ax.legend()
        ax.set_xlabel("Z")
        ax.set_ylabel('X')
        ax.set_zlabel("Y")
        plt.show()


    def plot(self):
        traj_array=self.convert_arrays()
        traj_array.plot()

    # draw all coordinate of the 2 trajectories in function of the time
    # It must have the same shape and the same fistr vector wich represent time.
    def plot_2_trajectory(self,traj2):
        traj_array=self.convert_arrays()
        if type(traj2)==TrajectoryAnalitic:
            traj2_array=traj2.convert_arrays()
            traj_array.plot_2_trajectory(traj2_array)
        else :
            traj_array.plot_2_trajectory(traj2)

    def copy(self):
        return TrajectoryAnalitic(t=self.t.copy(), x=self.x, y=self.y, z=self.z, v_x=self.v_x,
                          v_y=self.v_y, v_z=self.v_z, a_x=self.a_x,
                          a_y=self.a_y, a_z=self.a_z)

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