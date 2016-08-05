import numpy as np
import scipy.integrate as integrate
from abc import abstractmethod

RADIATION_GRID=0
RADIATION_LIST=1


class Radiation(object):

    def __init__(self, intensity, X, Y, distance,radiation_type=None):
        self.radiation_type=radiation_type
        self.intensity=intensity
        self.X=X
        self.Y=Y
        self.distance=distance

    def copy(self):
        return Radiation(intensity=self.intensity.copy(), X=self.X.copy(), Y=self.Y.copy(), distance=self.distance)



    def plot(self):
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D

        if self.X == None or self.Y == None:
            raise Exception(" X and Y must be array for plotting")
        if self.X.shape != self.Y.shape:
            raise Exception(" X and Y must have the same shape")
        fig = plt.figure()
        if len(self.X.shape) ==2 :
            ax = Axes3D(fig)
            ax.plot_surface(self.X, self.Y, self.intensity, rstride=1, cstride=1)
        else :
            ax = fig.gca(projection='3d')
            ax.plot(self.X, self.Y, self.intensity, label='brigthness')
            ax.legend()
        plt.xlabel('X')
        plt.ylabel('Y')
        plt.show()


    def difference_with(self,radiation2):
        if self.radiation_type != radiation2.radiation_type :
            raise Exception('difference between two different radiation type not define')
        if not self.XY_are_like_in(radiation2):
            raise Exception('X and Y must be the same for each radiation')
        error=np.abs(self.intensity - radiation2.intensity)
        res=self.copy()
        res.intensity=error
        return res

    def relativ_difference_with(self,radiation2):
        diff=self.difference_with(radiation2=radiation2)
        rad_max_ref=self.intensity.max()
        diff.intensity *= 1./rad_max_ref
        return diff

    def integration(self):
        if (self.X == None or self.Y == None):
            raise Exception(" X and Y must be array for integration")
        if self.X.shape != self.Y.shape:
            raise Exception(" X and Y must have the same shape")
        if len(self.X.shape)==2 :
            X = np.linspace(self.X[0].min(), self.X[0].max(), len(self.X[0]))
            Y = np.linspace(self.Y[:, 0].min(), self.Y[:, 0].max(),len(self.Y[:,0]))
            res1=integrate.trapz(self.intensity, X)
            res = integrate.trapz(res1, Y)
            #res = integrate.simps(integrate.simps(self.intensity, self.X), self.Y)
        else : # X and Y are 1d array
            if len(self.X) == 1:
                res = self.intensity[0]
            else: # choix arbitraire
                XY = np.zeros_like(self.X)
                for i in range(1, len(self.X)):
                    XY[i] = np.sqrt((self.X[i] - self.X[0]) * 2 + (self.Y[i] - self.Y[0]) ** 2)
                res = np.trapz(self.intensity, XY)
        return res


    def XY_are_like_in(self,rad2):
        if self.X.shape!= rad2.X.shape :
            return False
        elif self.Y.shape != rad2.Y.shape :
            return False
        else :
            if len(self.X.shape)==2 and len(self.X.shape)==2 :
                egality=np.zeros(self.X.shape[0],dtype=bool)
                for i in range(len(egality)):
                    egality[i]=all(self.X[i]==rad2.X[i])
                if not all(egality) :
                    return False
            else :
                 if not all(self.X==rad2.X):
                    return False
            if len(self.Y.shape)==2 and len(self.Y.shape)==2 :
                egality=np.zeros(self.X.shape[0],dtype=bool)
                for i in range(len(egality)):
                    egality[i]=all(self.Y[i]==rad2.Y[i])
                if not all(egality) :
                    return False
            else :
                 if not all(self.Y==rad2.Y):
                    return False
        return  True

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


    def change_Nb_pts(self,Nb_pts):
        if len(self.X.shape)==2 or len(self.Y.shape)==2 :
            X = np.linspace(self.X[0].min(), self.X[0].max(), Nb_pts)
            Y = np.linspace(self.Y[:, 0].min(), self.Y[:, 0].max(), Nb_pts)
            self.X, self.Y = np.meshgrid(X, Y)
        else :
            self.X = np.linspace(self.X.min(), self.X.max(), Nb_pts)
            self.Y= np.linspace(self.Y.min(), self.Y.max(), Nb_pts)


# def plot_wave(self, Nb_pts):
#     import matplotlib.pyplot as plt
#     from mpl_toolkits.mplot3d import Axes3D
#
#     if self.X == None or self.Y == None:
#         raise Exception(" X and Y must be grid or a list for plotting")
#     fig = plt.figure()
#     ax = fig.gca(projection='3d')
#     X = np.array([self.X[0]])
#     Y = np.array([self.Y[0]])
#     intensity = np.array([self.intensity[0]])
#     ax.plot(X, Y, intensity, '^', label='wave number 0')
#     wave_number = 1
#     while (wave_number * Nb_pts < len(self.X)):
#         X = self.X[(wave_number - 1) * Nb_pts + 1:wave_number * Nb_pts + 1]
#         Y = self.Y[(wave_number - 1) * Nb_pts + 1:wave_number * Nb_pts + 1]
#         intensity = self.intensity[(wave_number - 1) * Nb_pts + 1:wave_number * Nb_pts + 1]
#         ax.plot(X, Y, intensity, label='wave number %d' % wave_number)
#         wave_number += 1
#     ax.set_xlabel("X")
#     ax.set_ylabel('Y')
#     ax.set_zlabel("itensity")
#     ax.legend()
#     plt.show()


def Exemple_Grid():
    X = np.linspace(0.0, 0.005, 101)
    Y = np.linspace(0.0, 0.005, 101)
    X_grid, Y_grid = np.meshgrid(X, Y)
    intensity = (X_grid * 1e3) ** 2 * (Y_grid * 1e3) ** 2
    rad = Radiation(intensity=intensity, X=X_grid, Y=Y_grid, distance=100)
    print('maximum intensity radiation')
    print(rad.max())
    print(' integration of the intensity on the grid X,Y')
    print(rad.integration())
    print(' radiation intensity plot')
    rad.plot()

    print(' ')
    print('create a second radiation on the same grid')
    rad2 = rad.copy()
    rad2.intensity = (X_grid * 1e3) ** 2 + (Y_grid * 1e3) ** 2
    print('maximum intensity radiation')
    print(rad2.max())
    print(' integration of the intensity on the grid X,Y')
    print(rad2.integration())
    print(' radiation intensity plot')
    rad2.plot()

    print(' ')
    print('create a third radiation which is the different between the 2 radiation before')
    diff = rad.difference_with(rad2)
    print('maximum intensity radiation')
    print(diff.max())
    print(' integration of the intensity on the grid X,Y')
    print(diff.integration())
    print(' radiation intensity plot')
    diff.plot()



def Exemple_List():
    X = np.linspace(0.0, 0.005, 101)
    Y = np.linspace(0.0025, 0.0025, 101)
    distance = 100

    intensity = (X * 1e3) ** 2 * (Y * 1e3) ** 2
    rad = Radiation(intensity=intensity, X=X, Y=Y, distance=distance)
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

if __name__ == "__main__" :

    Exemple_Grid()
    Exemple_List()
