from pylab import *
import matplotlib.pyplot as plt
import scipy.constants as codata
from mpl_toolkits.mplot3d import Axes3D


class Radiation(object):
    def __init__(self, map,X,Y,distance,parameters):
        self.intensity=map
        self.X=X
        self.Y=Y
        self.distance=distance
        self.parameters=parameters


    # draw all coordinate of the trajectory in function of the time
    def draw(self):
        if len(self.X.shape) ==1 :
            X_grid,Y_grid = np.meshgrid(self.X, self.Y)
        else :
            X_grid=self.X.copy()
            Y_grid=self.Y.copy()
        fig = figure()
        ax = Axes3D(fig)
        ax.plot_surface(X_grid, Y_grid, self.intensity, rstride=1, cstride=1)
        ax.set_xlabel("X")
        ax.set_ylabel('Y')
        ax.set_zlabel("flux")
        show()

    def difference_with(self,radiation2):
      error=np.abs(self.map - radiation2.map)
      res=Radiation(error,self.X.copy().self.Y.copy(),self.parameters.copy())
      return res

    def change_distance(self,D):
        self.distance=D
        #update intensity
        self.intensity = self.parameters.calculate_radiation_intensity(distance=self.distance,
                                                                       X_arrays=self.X, Y_arrays=self.Y)


    def compare_with_ditance(self,radiation2,D):
        error_max=np.zeros(len(D))
        for i in range(len(D)):
            print(i)
            self.change_distance(D[i])
            radiation2.change_distance(D[i])
            error_max[i]=(np.abs(self.intensity - radiation2.intensity)).max()
        return error_max




