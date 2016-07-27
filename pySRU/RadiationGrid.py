import numpy as np
import scipy.integrate as integrate
from pySRU.Radiation import Radiation , RADIATION_LIST,RADIATION_GRID


class RadiationGrid(Radiation):

    def __init__(self, intensity, X, Y, distance):
        super(self.__class__, self).__init__(intensity=intensity,X=X,Y=Y,distance=distance,
                                             radiation_type=RADIATION_GRID)

    def copy(self):
        return RadiationGrid(intensity=self.intensity.copy(), X=self.X.copy(), Y=self.Y.copy(), distance=self.distance)

    def plot(self):
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D
        if self.X== None or self.Y== None:
            raise Exception(" X and Y must be grid or a list for plotting")
        fig = plt.figure()
        ax = Axes3D(fig)
        ax.plot_surface(self.X, self.Y, self.intensity, rstride=1, cstride=1)
        plt.show()

    def XY_are_like_in(self,rad2):
        if self.X.shape!= rad2.X.shape :
            return False
        elif self.Y.shape != rad2.Y.shape :
            return False
        else :
            egality=np.zeros(self.X.shape[0],dtype=bool)
            print(egality.shape)
            for i in range(len(egality)):
                egality[i]=all(self.X[i]==rad2.X[i])
            if not all(egality) :
                return False
            else :
                egality =np.zeros(self.Y.shape[0],dtype=bool)
                for i in range(len(egality)) :
                    egality[i] = all(self.Y[i] == rad2.Y[i])
                res=all(egality)
        return res


    def integration(self):
        if (self.X == None or self.Y == None):
            print("pb ds Radiation . integration")
            res = 0
        else:
            Nb_pts=len(self.X) # TODO pas sure que ca marche pour X != Y
            X = np.linspace(self.X.min(),self.X.max(),Nb_pts)
            Y = np.linspace(self.Y.min(), self.Y.max(),Nb_pts)
            res1=integrate.trapz(self.intensity, X)
            res = integrate.trapz(res1, Y)
        return res

    def change_Nb_pts(self, Nb_pts):
        X= np.linspace(self.X[0].min(), self.X[0].max(), Nb_pts)
        Y= np.linspace(self.Y[0].min(), self.Y[0].max(), Nb_pts)
        self.X,self.Y=np.meshgrid(X,Y)


if __name__ == "__main__" :

    X=np.linspace(0.0,0.005,101)
    Y=np.linspace(0.0,0.005,101)
    distance=100
    X_grid,Y_grid=np.meshgrid(X,Y)
    intensity = (X_grid*1e3) ** 2 * (Y_grid*1e3) ** 2
    rad=RadiationGrid(intensity=intensity, X=X_grid, Y=Y_grid,distance=100)
    print('maximum intensity radiation')
    print(rad.max())
    print(' integration of the intensity on the grid X,Y')
    print(rad.integration())
    print(' radiation intensity plot')
    rad.plot()

    print(' ')
    print('create a second radiation on the same grid')
    rad2=rad.copy()
    rad2.intensity = (X_grid*1e3) ** 2 + (Y_grid*1e3) ** 2
    print('maximum intensity radiation')
    print(rad2.max())
    print(' integration of the intensity on the grid X,Y')
    print(rad2.integration())
    print(' radiation intensity plot')
    rad2.plot()

    print(' ')
    print('create a third radiation which is the different between the 2 radiation before')
    diff=rad.difference_with(rad2)
    print('maximum intensity radiation')
    print(diff.max())
    print(' integration of the intensity on the grid X,Y')
    print(diff.integration())
    print(' radiation intensity plot')
    diff.plot()