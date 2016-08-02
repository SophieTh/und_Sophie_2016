import numpy as np
import scipy.integrate as integrate
from pySRU.Radiation import Radiation,RADIATION_GRID,RADIATION_LIST


class RadiationList(Radiation):

    def __init__(self, intensity, X, Y, distance):
        super(self.__class__, self).__init__(intensity=intensity,X=X,Y=Y,distance=distance,
                                             radiation_type=RADIATION_LIST)

    def copy(self):
        return RadiationList(intensity=self.intensity.copy(), X=self.X.copy(), Y=self.Y.copy(), distance=self.distance)

    def plot(self):
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D

        if self.X == None or self.Y == None:
            raise Exception(" X and Y must be array for plotting")
        if self.X.shape != self.Y.shape:
            raise Exception(" X and Y must have the same shape")
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        ax.plot(self.X, self.Y, self.intensity)
        plt.show()


    def plot_wave(self,Nb_pts):
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D

        if self.X == None or self.Y == None:
            raise Exception(" X and Y must be grid or a list for plotting")
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        X=np.array([self.X[0]])
        Y = np.array([self.Y[0]])
        intensity=np.array([self.intensity[0]])
        ax.plot(X, Y, intensity, label='wave number 0')
        wave_number=1
        while (wave_number*Nb_pts < len(self.X) ):
            X=self.X[(wave_number-1)*Nb_pts+1:wave_number*Nb_pts+1]
            Y=self.Y[(wave_number-1)*Nb_pts+1:wave_number*Nb_pts+1]
            intensity=self.intensity[(wave_number-1)*Nb_pts+1:wave_number*Nb_pts+1]
            ax.plot(X, Y, intensity , label='wave number %d'%wave_number)
            wave_number += 1
        ax.set_xlabel("X")
        ax.set_ylabel('Y')
        ax.set_zlabel("itensity")
        ax.legend()
        plt.show()


    def XY_are_like_in(self,rad2):
        if any(self.X != rad2.X) :
            return False
        else :
            if  any(self.Y != rad2.Y) :
                return False
            else :
                return True


    def integration(self):
        if (self.X == None or self.Y == None):
            raise Exception ("X and Y must be define ")
        if len(self.X) == 1:
            res = self.intensity[0]
        else:
            XY = np.zeros_like(self.X)
            for i in range(1,len(self.X)) :
                XY[i]=np.sqrt((self.X[i]-self.X[0])*2+(self.Y[i]-self.Y[0])**2)
            res=np.trapz(self.intensity,XY)
        return res

    def change_Nb_pts(self, Nb_pts):
        self.X= np.linspace(self.X.min(), self.X.max(), Nb_pts)



if __name__ == "__main__":
    X =np.linspace(0.0, 0.005, 101)
    Y =np.linspace(0.0025, 0.0025, 101)
    distance = 100

    intensity = (X * 1e3) ** 2 * (Y* 1e3) ** 2
    rad = RadiationList(intensity=intensity, X=X, Y=Y, distance=distance)
    print('maximum intensity radiation')
    print(rad.max())
    print(' integration of the intensity on X,Y')
    print(rad.integration())
    print(' radiation intensity plot')
    rad.plot()

    print(' ')
    print('create a second radiation on the same X,Y')
    rad2 = rad.copy()
    rad2.intensity = (X * 1e3) ** 2 + (Y * 1e3) ** 2
    print('maximum intensity radiation')
    print(rad2.max())
    print(' integration of the intensity on X,Y')
    print(rad2.integration())
    print(' radiation intensity plot')
    rad2.plot()

    print(' ')
    print('create a third radiation which is the different between the 2 radiation before')
    diff = rad.difference_with(rad2)
    print('maximum intensity radiation')
    print(diff.max())
    print(' integration of the intensity on X,Y')
    print(diff.integration())
    print(' radiation intensity plot')
    diff.plot()