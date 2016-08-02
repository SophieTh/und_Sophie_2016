import numpy as np
import matplotlib.pyplot as plt





class MagneticField(object):
    def __init__(self,Bx,By,Bz):
        self.Bx = Bx
        self.By = By
        self.Bz = Bz

    def copy(self): #TODO a changer ?
        # if type(self.x)==np.ndarray:
        #     x = self.x.copy()
        # else :
        #     x=self.x
        # if type(self.y)==np.ndarray:
        #     y = self.y.copy()
        # else :
        #     y=self.y
        # if type(self.z)==np.ndarray:
        #     z = self.z.copy()
        # else :
        #    z=self.z
        if type(self.Bx)==np.ndarray:
            Bx = self.Bx.copy()
        else :
            Bx=self.Bx
        if type(self.By)==np.ndarray:
            By = self.By.copy()
        else :
            By=self.By
        if type(self.Bz)==np.ndarray:
            Bz = self.Bz.copy()
        else :
            Bz=self.Bz
        return MagneticField(Bx=Bx, By=By, Bz=Bz)

    #faire une fct aui fait l'interpolation  et elargissement
    #

    # enelarge 2 vector for the interpolation of the magnetic field use in trajectory_undulator_from_magnetic_field2
    def enlargement_vector_for_interpolation(self,Z,Y ,nb_enlarg=0):
        dz = Z[1] -Z[0]
        enlarg_1 = np.linspace(Z[0] - nb_enlarg * dz, Z[0] - dz, nb_enlarg)
        enlarg_2 = np.linspace(Z[- 1] + dz, Z[- 1] + nb_enlarg * dz, nb_enlarg)
        Z = np.concatenate((enlarg_1, Z))
        Z = np.concatenate((Z, enlarg_2))
        if (type(Y )==np.ndarray) :
            dy = Y [1] - Y [0]
            enlarg_1 = np.linspace(Y [0] - nb_enlarg * dz, Y [0] - dy, nb_enlarg)
            enlarg_2 = np.linspace(Y [- 1] + dy, Y [- 1] + nb_enlarg * dy, nb_enlarg)
            Y  = np.concatenate((enlarg_1, Y ))
            Y = np.concatenate((Y , enlarg_2))
        if (type(self.By)==np.ndarray) :
            enlarg_3 = np.zeros(nb_enlarg)
            self.By = np.concatenate((enlarg_3, self.By))
            self.By = np.concatenate((self.By, enlarg_3))

# rajouter un code qui qffiche la legende ??
    def plot_z(self,X,Y,Z):
        import matplotlib.pyplot as plt

        plt.plot(Z, self.Bx(z=Z,y=Y,x=X))
        plt.title(" Bx = f(z) ")
        plt.ylabel('Bx')
        plt.xlabel('Z')
        plt.show()

        plt.plot(Z, self.By(z=Z,y=Y,x=X))
        plt.title(" By = f(z) ")
        plt.ylabel('By')
        plt.xlabel('Z')
        plt.show()

        plt.plot(Z, self.Bz(z=Z,y=Y,x=X))
        plt.title(" Bz = f(z) ")
        plt.ylabel('Bz')
        plt.xlabel('Z')
        plt.show()

# # TODO A CHANGER COMME PLOT_z
        #TODO quiver plot ?
#     def plot_y(self):
#         plt.plot(self.y, self.Bx(self.z, self.y))
#         plt.title(" Bx = f(y) ")
#         plt.xlabel('Bx')
#         plt.ylabel('y')
#         plt.show()
#
#         plt.plot(self.y, self.By(self.z, self.y))
#         plt.title(" By = f(y) ")
#         plt.xlabel('By')
#         plt.ylabel('y')
#         plt.show()
#
#         plt.plot(self.y, self.Bz(self.z, self.y))
#         plt.title(" Bz = f(y) ")
#         plt.xlabel('Bz')
#         plt.ylabel('y')
#         plt.show()