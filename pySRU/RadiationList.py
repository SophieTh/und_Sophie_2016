import numpy as np
import scipy.integrate as integrate
from pySRU.Radiation import Radiation,RADIATION_GRID,RADIATION_LIST


class RadiationList(Radiation):

    def __init__(self, intensity, X, Y, distance):
        super(self.__class__, self).__init__(intensity=intensity,X=X,Y=Y,distance=distance,
                                             radiation_type=RADIATION_LIST)


    def plot(self):
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D

        if self.X == None or self.Y == None:
            raise Exception(" X and Y must be grid or a list for plotting")
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        ax.plot(self.X, self.Y, self.intensity, label='wave number')
        ax.legend()
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


    #TODO a tester ! et changer !
    def integration(self):
        if (self.X == None or self.Y == None):
            raise Exception ("X and Y must be define ")
        if len(self.X.shape)==1:
            if len(self.X) == 1:
                res = self.X[0]
            else:
                XY = np.zeros_like(self.X)
                for i in range(1,len(self.X)) :
                    XY[i]=np.linalg.norm([self.X[i]-self.X[i-1],self.Y[i]-self.Y[i-1]])
                res=np.trapz(self.intensity,XY)
        else : # on suppose qu'on ne va pas plus loin que 2 et que tout est a la bonne shape #TODO
            integral=np.zeros_like(self.X)
            for i in range(self.X.shape[0]):
                X=self.X[i]
                Y=self.Y[i]
                if len(X)==1 :
                    integral[i]==X[0]
                else :
                    XY = np.zeros_like(X)
                    for j in range(1, len(X)):
                        XY[j] = np.linalg.norm([X[j] - X[j - 1], Y[j] - Y[j - 1]])
                integral[i] = np.trapz(self.intensity[i], XY)
                #TODO criticable ....
            res=np.sum(integral)
        return res



    def change_Nb_pts(self, Nb_pts):
        self.X= np.linspace(self.X.min(), self.X.max(), Nb_pts)


#TODO
if __name__ == "__main__":
    pass