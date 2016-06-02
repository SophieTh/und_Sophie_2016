import numpy as np
import matplotlib.pyplot as plt



class Trajectory(object):
    def __init__(self, t, x, y, z, v_x, v_y, v_z, a_x, a_y, a_z):
        self.t = t.copy()
        self.x = x.copy()
        self.y = y.copy()
        self.z = z.copy()
        self.v_x = v_x.copy()
        self.v_y = v_y.copy()
        self.v_z = v_z.copy()
        self.a_x = a_x.copy()
        self.a_y = a_y.copy()
        self.a_z = a_z.copy()

    def nb_points(self):
        N=len(self.t)
        if (len(self.x)==N and len(self.y)==N
            and len(self.z)==N and len(self.v_x)==N
            and len(self.v_y) == N  and len(self.v_z)==N
            and len(self.a_x)==N and len(self.a_y)==N
            and len(self.a_z)==N ) :
            N = len(self.t)
        else :
            N=0
        return N

    # def error(self,trajec_test):
    #     if self.t != trajec_test.t :
    #         raise Exception("Problem : the two trajectory have not the same vector t.")
    #     error=np.zeros((10,self.nb_points()))
    #     error[0]=

    def copy(self):
        return Trajectory(t=self.t.copy(), x=self.x.copy(), y=self.y.copy(), z=self.z.copy(), v_x=self.v_x.copy(),
                          v_y=self.v_y.copy(), v_z=self.v_z.copy(), a_x=self.a_x.copy(),
                          a_y=self.a_y.copy(), a_z=self.a_z.copy())

    # draw all coordinate of the trajectory in function of the time
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

