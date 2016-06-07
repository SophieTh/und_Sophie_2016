import numpy as np
import scipy.integrate as integrate
import matplotlib.pyplot as plt
from pylab import *
from mpl_toolkits.mplot3d import Axes3D

from TrajectoryFactory import TrajectoryFactory, TRAJECTORY_METHOD_ANALYTIC
from UndulatorParameter import UndulatorParameters as Undulator

class Radiation(object):
    def __init__(self, map,X,Y,distance=None):
        self.intensity=map
        self.X=X
        self.Y=Y
        self.distance=distance

    def copy(self):
        return Radiation( map=self.intensity.copy(),X=self.X.copy(),Y=self.Y.cpy(),distance=self.distance)

    # draw all coordinate of the trajectory in function of the time
    def draw(self):
        if self.X== None or self.Y== None:
            print("oups")
            if self.distance == None:
                distance = 100
            else :
                distance=self.distance
            self.X = np.arange(0.0, (distance * 1.01) * 1e-3, distance * 1e-5)
            self.Y = np.arange(0.0, (distance * 1.01) * 1e-3, distance * 1e-5)
        if len(self.X.shape) ==1 :
            X_grid,Y_grid = np.meshgrid(self.X, self.Y)
        else :
            X_grid=self.X.copy()
            Y_grid=self.Y.copy()
        fig = plt.figure()
        ax = Axes3D(fig)
        ax.plot_surface(X_grid, Y_grid, self.intensity, rstride=1, cstride=1)
        ax.set_xlabel("X")
        ax.set_ylabel('Y')
        ax.set_zlabel("flux")
        plt.show()

    def difference_with(self,radiation2):
      error=np.abs(self.map - radiation2.map)
      res=Radiation(error,self.X.copy().self.Y.copy())
      return res


    def max(self):
        return self.intensity.max()

    def error_max(self,radiation2):
        return (np.abs(self.intensity - radiation2.intensity)).max()

    def relativ_error(self,radiation2):
        if (self.max() ==0.0) :
           raise Exception("Problem : radiation max is null")
        res=self.error_max(radiation2)/self.max()
        return res

    def integration(self):
        if (self.X == None or self.Y == None ) :
            print("pb ds Radtion . integration")
            res = 0
        else :
            res=integrate.simps(integrate.simps(self.intensity,self.X),self.Y)
        return res


if __name__ == "__main__" :

    X=np.linspace(0.0,0.005,101)
    Y=np.linspace(0.0,0.005,101)
    X_grid,Y_grid=np.meshgrid(X,Y)
    map = X_grid ** 2 * Y_grid ** 2
    rad=Radiation(map=map,X=X,Y=Y)
    print(rad.intensity.shape)
    print(rad.integration())
    rad.draw()

