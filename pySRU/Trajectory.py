import numpy as np
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

    def nb_points(self):
        N=len(self.t)
        if (len(self.x)==N and len(self.y)==N
            and len(self.z)==N and len(self.v_x)==N
            and len(self.v_y) == N  and len(self.v_z)==N
            and len(self.a_x)==N and len(self.a_y)==N
            and len(self.a_z)==N ) :
            N = len(self.t)
        else :
            #raise Exception("different lenght in trajectory.")
            print("different lenght in trajectory.")
            N=0
        return N

    def convert(self,T):
        if T.shape[0]== 10 :
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


    def error(self,trajec_test):
        if all(self.t != trajec_test.t) :
            raise Exception("Problem : the two trajectory have not the same vector t.")
        error=np.zeros((10,self.nb_points()))
        error[0]=self.t.copy()
        error[1] = np.abs(self.x - trajec_test.x)
        error[2] = np.abs(self.y - trajec_test.y)
        error[3] = np.abs(self.z - trajec_test.z)
        error[4] = np.abs(self.v_x - trajec_test.v_x)
        error[5] = np.abs(self.v_y - trajec_test.v_y)
        error[6] = np.abs(self.v_z - trajec_test.v_z)
        error[7] = np.abs(self.a_x - trajec_test.a_x)
        error[8] = np.abs(self.a_y - trajec_test.a_y)
        error[9] = np.abs(self.a_z - trajec_test.a_z)
        return error
        # return Trajectory(error[0],error[1],error[2],error[3],error[4],error[5],error[6]
        #                   ,error[7],error[8],error[9])

    def error_max(self,trajec_test):
        if all(self.t != trajec_test.t) :
            raise Exception("Problem : the two trajectory have not the same vector t.")
        error=np.zeros((10))
        error[0]=self.nb_points()
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


    def copy(self):
        return Trajectory(t=self.t.copy(), x=self.x.copy(), y=self.y.copy(), z=self.z.copy(), v_x=self.v_x.copy(),
                          v_y=self.v_y.copy(), v_z=self.v_z.copy(), a_x=self.a_x.copy(),
                          a_y=self.a_y.copy(), a_z=self.a_z.copy())

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


    def draw(self):
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
    def draw_2_trajectory(self,traj2):
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
    Traj.draw()