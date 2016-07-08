import numpy as np
import matplotlib.pyplot as plt
from Trajectory import Trajectory




class TrajectoryArray(Trajectory):

    def __init__(self, t, x, y, z, v_x, v_y, v_z, a_x, a_y, a_z):
        super(self.__class__, self).__init__(t,x,y,z,v_x,v_y,v_z,a_x,a_y,a_z)


    def nb_points(self):
        N = len(self.t)
        if (len(self.x) == N and len(self.y) == N
            and len(self.z) == N and len(self.v_x) == N
            and len(self.v_y) == N and len(self.v_z) == N
            and len(self.a_x) == N and len(self.a_y) == N
            and len(self.a_z) == N):
            N = len(self.t)
        else: # raise Exception("different lenght in trajectory.")
            print("different lenght in trajectory.")
            N = 0
        return N


    def convert(self, T):
        if T.shape[0] == 10:
            self.t = T[0]
            self.x = T[1]
            self.y = T[2]
            self.z = T[3]
            self.v_x = T[4]
            self.v_y = T[5]
            self.v_z = T[6]
            self.a_x = T[7]
            self.a_y = T[8]
            self.a_z = T[9]


    def multiply_by(self, cst):
        self.x *= cst
        self.y *= cst
        self.z *= cst
        self.v_x *= cst
        self.v_y *= cst
        self.v_z *= cst
        self.a_x *= cst
        self.a_y *= cst
        self.a_z *= cst


    def error(self, trajectory_test):
        if type(trajectory_test)!=TrajectoryArray :
            if any((self.t != trajectory_test.t) ):
                print("the time vecto change for the analitical trajectory")
                trajectory_test.t=self.t
            trajec_test=trajectory_test.convert()
        elif any(np.abs(self.t - trajectory_test.t) > np.abs(self.t).max()*1e-6):
            raise Exception("Problem : the two trajectory have not the same vector t ??")
        else :
            trajec_test=trajectory_test
        error = np.zeros((10, self.nb_points()))
        error[0] = self.t.copy()
        error[1] = np.abs(self.x - trajec_test.x)
        error[2] = np.abs(self.y - trajec_test.y)
        error[3] = np.abs(self.z - trajec_test.z)
        error[4] = np.abs(self.v_x - trajec_test.v_x)
        error[5] = np.abs(self.v_y - trajec_test.v_y)
        error[6] = np.abs(self.v_z - trajec_test.v_z)
        error[7] = np.abs(self.a_x - trajec_test.a_x)
        error[8] = np.abs(self.a_y - trajec_test.a_y)
        error[9] = np.abs(self.a_z - trajec_test.a_z)
        return TrajectoryArray(error[0],error[1],error[2],error[3],error[4],error[5],error[6]
                           ,error[7],error[8],error[9])


    def error_max(self, trajectory_test):
        if type(trajectory_test) !=TrajectoryArray :
            if any((self.t != trajectory_test.t) ):
                print("the time vector change for the analitical trajectory")
                trajectory_test.t=self.t
            trajec_test=trajectory_test.convert()
        elif all((self.t - trajectory_test.t) <= np.abs(self.t).max() * 1e-6):
            trajec_test = trajectory_test

        else :
            raise Exception("Problem : the two trajectory have not the same vector t ??")
        error = np.zeros(10)
        error[0] = self.nb_points()
        error[1] = (np.abs(self.x - trajec_test.x)).max()
        error[2] = (np.abs(self.y - trajec_test.y)).max()
        error[3] = (np.abs(self.z - trajec_test.z)).max()
        error[4] = (np.abs(self.v_x - trajec_test.v_x)).max()
        error[5] = (np.abs(self.v_y - trajec_test.v_y)).max()
        error[6] = (np.abs(self.v_z - trajec_test.v_z)).max()
        error[7] = (np.abs(self.a_x - trajec_test.a_x)).max()
        error[8] = (np.abs(self.a_y - trajec_test.a_y)).max()
        error[9] = (np.abs(self.a_z - trajec_test.a_z)).max()
        return error


    def error_rel(self, trajectory_test):
        if type(trajectory_test)==TrajectoryAnalitic :
            if any((self.t != trajectory_test.t) ):
                print("the time vecto change for the analitical trajectory")
                trajectory_test.t=self.t
            trajec_test=trajectory_test.convert_array()
        elif all((self.t - trajectory_test.t) <= np.abs(self.t).max() * 1e-6):
            raise Exception("Problem : the two trajectory have not the same vector t ??")
        else :
            trajec_test=trajectory_test
        error = self.error(trajec_test)
        for i in range(self.nb_points()):
            if self.x[i] != 0.0:
                error[1][i] *= 1.0 / (np.abs(self.x[i]))
            # else :
            #     error[1][i]= -1

            if self.y[i] != 0.0:
                error[2][i] *= 1.0 / (np.abs(self.y[i]))
            # else :
            #     error[2][i]= -1

            if self.z[i] != 0.0:
                error[3][i] *= 1.0 / (np.abs(self.z[i]))
            # else :
            #     error[3][i]= -1

            if self.v_x[i] != 0.0:
                error[4][i] *= 1.0 / (np.abs(self.v_x[i]))
            # else:
            #     error[4][i] = -1

            if self.v_y[i] != 0.0:
                error[5][i] *= 1.0 / (np.abs(self.v_y[i]))
            # else:
            #     error[5][i] = -1

            if self.v_z[i] != 0.0:
                error[6][i] *= 1.0 / (np.abs(self.v_z[i]))
            # else:
            #     error[6][i] = -1

            if self.a_x[i] != 0.0:
                error[7][i] *= 1.0 / (np.abs(self.a_x[i]))
            # else:
            #     error[7][i] = -1

            if self.a_y[i] != 0.0:
                error[8][i] *= 1.0 / (np.abs(self.a_y[i]))
            # else:
            #     error[8][i] = -1

            if self.a_z[i] != 0.0:
                error[9][i] *= 1.0 / (np.abs(self.a_z[i]))
                # else:
                #     error[9][i] = -1

        return TrajectoryArray(error[0], error[1], error[2], error[3], error[4], error[5], error[6]
                          , error[7], error[8], error[9])

    def error_rel_max(self, trajectory_test):
        if type(trajectory_test) != TrajectoryArray:
            if any((self.t != trajectory_test.t)):
                print("the time vecto change for the analitical trajectory")
                trajectory_test.t = self.t
            trajec_test = trajectory_test.convert_array()
        elif any((self.t - trajectory_test.t) > np.abs(self.t).max() * 1e-6):
            raise Exception("Problem : the two trajectory have not the same vector t ??")
        else:
            trajec_test = trajectory_test
        error = self.error(trajec_test)
        xm=np.abs(self.x).max()
        if xm != 0.0:
            error.x *= 1.0 / xm

        ym = np.abs(self.y).max()
        if ym != 0.0 :
            error.y *= 1.0 /ym

        zm = np.abs(self.z).max()
        if zm != 0.0 :
            error.z *= 1.0 /zm

        vxm = np.abs(self.v_x).max()
        if vxm != 0.0:
            error.v_x *= 1.0 / vxm

        vym = np.abs(self.v_y).max()
        if vym != 0.0:
            error.v_y *= 1.0 / vym

        vzm = np.abs(self.v_z).max()
        if vzm != 0.0:
            error.v_z *= 1.0 / vzm

        axm = np.abs(self.a_x).max()
        if axm != 0.0:
            error.a_x *= 1.0 / axm

        aym = np.abs(self.a_y).max()
        if aym != 0.0:
            error.a_y *= 1.0 / aym

        azm = np.abs(self.a_z).max()
        if azm != 0.0:
            error.a_z *= 1.0 / azm

        return error


    # draw all coordinate of the trajectory in function of the time
    def plot_3D(self):
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        ax.plot(self.z, self.x, self.y, label='trajectory')
        ax.legend()
        ax.set_xlabel("Z")
        ax.set_ylabel('X')
        ax.set_zlabel("Y")
        plt.show()


    def plot(self):
        plt.plot(self.t, self.x)
        plt.title(" X = f(t) ")
        plt.xlabel('t')
        plt.ylabel('X')
        plt.show()

        plt.plot(self.t, self.y)
        plt.title(" Y = f(t) ")
        plt.xlabel('t')
        plt.ylabel('V')
        plt.show()

        plt.plot(self.t, self.z)
        plt.title(" Z  = f(t) ")
        plt.xlabel('t')
        plt.ylabel('Z')
        plt.show()

        print('average of Vz  =')
        Beta_et = np.sum(self.v_z) / len(self.v_z)
        print(Beta_et)

        Z = Beta_et * self.t
        plt.plot(self.t, self.z - Z)
        plt.title(" Z - Beta* t = f(t) ")
        plt.xlabel('t')
        plt.ylabel('Z - Beta* t')
        plt.show()

        plt.plot(self.t, self.v_x)
        plt.title(" Vx = f(t) ")
        plt.xlabel('t')
        plt.ylabel('Vx')
        plt.show()

        plt.plot(self.t, self.v_y)
        plt.title(" Vy = f(t) ")
        plt.xlabel('t')
        plt.ylabel('Vy')
        plt.show()

        plt.plot(self.t, self.v_z)
        plt.title(" Vz = f(t) ")
        plt.xlabel('t')
        plt.ylabel('Vz')
        plt.show()

        plt.plot(self.t, self.v_z - Beta_et)
        plt.title(" Vz -Beta* = f(t) ")
        plt.xlabel('t')
        plt.ylabel('Vz')
        plt.show()

        plt.plot(self.t, self.a_x)
        plt.title(" Ax = f(t) ")
        plt.xlabel('t')
        plt.ylabel('Ax')
        plt.show()

        plt.plot(self.t, self.a_y)
        plt.title(" Ay = f(t) ")
        plt.xlabel('t')
        plt.ylabel('Ay')
        plt.show()

        plt.plot(self.t, self.a_z)
        plt.title(" Az = f(t) ")
        plt.xlabel('t')
        plt.ylabel('Az')
        plt.show()


    # draw all coordinate of the 2 trajectories in function of the time
    # It must have the same shape and the same fistr vector wich represent time.
    def plot_2_trajectory(self, traj2):
        plt.plot(self.t, self.x)
        plt.plot(traj2.t, traj2.x)
        plt.title(" X = f(t) ")
        plt.xlabel('t')
        plt.ylabel('X')
        plt.show()

        plt.plot(self.t, self.y)
        plt.plot(traj2.t, traj2.y)
        plt.title(" Y = f(t) ")
        plt.xlabel('t')
        plt.ylabel('V')
        plt.show()

        plt.plot(self.t, self.z)
        plt.plot(traj2.t, traj2.z)
        plt.title(" Z  = f(t) ")
        plt.xlabel('t')
        plt.ylabel('Z')
        plt.show()

        Beta_et = np.sum(self.v_z) / len(self.v_z)
        plt.plot(self.t, self.z - Beta_et * self.t)
        plt.plot(traj2.t, traj2.z - Beta_et * self.t)
        plt.title(" Z -Beta* t = f(t) ")
        plt.xlabel('t')
        plt.ylabel('Z')
        plt.show()

        plt.plot(self.t, self.v_x)
        plt.plot(traj2.t, traj2.v_x)
        plt.title(" Vx = f(t) ")
        plt.xlabel('t')
        plt.ylabel('Vx')
        plt.show()

        plt.plot(self.t, self.v_y)
        plt.plot(traj2.t, traj2.v_y)
        plt.title(" Vy = f(t) ")
        plt.xlabel('t')
        plt.ylabel('Vy')
        plt.show()

        plt.plot(self.t, self.v_z)
        plt.plot(traj2.t, traj2.v_z)
        plt.title(" Vz = f(t) ")
        plt.xlabel('t')
        plt.ylabel('Vz')
        plt.show()

        plt.plot(self.t, self.a_x)
        plt.plot(traj2.t, traj2.a_x)
        plt.title(" Ax = f(t) ")
        plt.xlabel('t')
        plt.ylabel('Ax')
        plt.show()

        plt.plot(self.t, self.a_y)
        plt.plot(traj2.t, traj2.a_y)
        plt.title(" Ay = f(t) ")
        plt.xlabel('t')
        plt.ylabel('Ay')
        plt.show()

        plt.plot(self.t, self.a_z)
        plt.plot(traj2.t, traj2.a_z)
        plt.title(" Az = f(t) ")
        plt.xlabel('t')
        plt.ylabel('Az')
        plt.show()

    def copy(self):
        return TrajectoryArray(t=self.t.copy(), x=self.x.copy(), y=self.y.copy(), z=self.z.copy(), v_x=self.v_x.copy(),
                          v_y=self.v_y.copy(), v_z=self.v_z.copy(), a_x=self.a_x.copy(),
                          a_y=self.a_y.copy(), a_z=self.a_z.copy())