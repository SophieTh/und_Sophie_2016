import numpy as np
import scipy.integrate as integrate
from abc import abstractmethod

RADIATION_GRID=0
RADIATION_LIST=1


class Radiation(object):

    def __init__(self, intensity, X, Y, distance,radiation_type):
        self.radiation_type=radiation_type
        self.intensity=intensity
        self.X=X
        self.Y=Y
        self.distance=distance

    @abstractmethod
    def copy(self):
       return

    @abstractmethod
    def plot(self):
        pass
        # import matplotlib.pyplot as plt
        # from mpl_toolkits.mplot3d import Axes3D
        #
        # if self.X== None or self.Y== None:
        #     raise Exception(" X and Y must be grid or a list for plotting")
        #
        # fig = plt.figure()
        # if len(self.X.shape) ==2 :
        #     print('good')
        #     ax = Axes3D(fig)
        #     ax.plot_surface(self.X, self.Y, self.intensity, rstride=1, cstride=1)
        # else :
        #     ax = fig.gca(projection='3d')
        #     ax.plot(self.X, self.Y, self.intensity, label='radiation')
        #     ax.legend()
        # plt.show()


    def difference_with(self,radiation2):
        if self.radiation_type != radiation2.radiation_type :
            raise Exception('difference between two different radiation type not define')
        if not self.XY_are_like_in(radiation2):
            raise Exception('X and Y must be the same for each radiation')
        error=np.abs(self.intensity - radiation2.intensity)
        res=self.copy()
        res.intensity=error
        return res

    @abstractmethod
    def integration(self):
        return
        #
        # if (self.X == None or self.Y == None):
        #     print("pb ds Radiation . integration")
        #     res = 0
        # else:
        #     # print('(self.X).shape')
        #     # print((self.X).shape)
        #     # print('(self.Y).shape')
        #     # print((self.Y).shape)
        #     # print('(self.intensity).shape')
        #     # print((self.intensity).shape)
        #     # print('integrate.simps(self.intensity,self.X).shape')
        #     # print((integrate.simps(self.intensity,self.X)).shape)
        #     res = integrate.simps(integrate.simps(self.intensity, self.X), self.Y)
        # return res

    @abstractmethod
    def XY_are_like_in(self,rad2):
        return

    def max(self):
        res=(self.intensity).max()
        return res

    def error_max(self,radiation2):
        if self.intensity.shape != radiation2.intensity.shape :
            raise Exception('the two radiation must have the same shape')
        return (np.abs(self.intensity - radiation2.intensity)).max()

    def relativ_error(self,radiation2):
        if (self.max() ==0.0) :
           raise Exception("Problem : radiation max is null")
        res=self.error_max(radiation2)/self.max()
        return res

    @abstractmethod
    def change_Nb_pts(self,Nb_pts):
        return





