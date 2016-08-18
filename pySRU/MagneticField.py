import numpy as np
import matplotlib.pyplot as plt





class MagneticField(object):
    def __init__(self,Bx,By,Bz):
        self.Bx = Bx
        self.By = By
        self.Bz = Bz

    def copy(self):
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

    #TODO a modifier
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
        return Z,Y




    def plot_z(self,X,Y,Z,title="Magnetic field",field_component=None):
        import matplotlib.pyplot as plt

        if ( (field_component == None) or (field_component == 0) or (field_component == 'x') ):
            plt.plot(Z, self.Bx(z=Z,y=Y,x=X))

            plt.ylabel('Bx [T]')
            plt.xlabel('Z [m]')
            plt.title(title+" Bx = f(z) ")
            plt.show()

        if ( (field_component == None) or (field_component == 1) or (field_component == 'y') ):
            plt.plot(Z, self.By(z=Z,y=Y,x=X))
            plt.ylabel('By [T]')
            plt.xlabel('Z [m]')
            plt.title(title+" By = f(z) ")
            plt.show()

        if ( (field_component == None) or (field_component == 2) or (field_component == 'z') ):
            plt.plot(Z, self.Bz(z=Z,y=Y,x=X))
            plt.ylabel('Bz [T]')
            plt.xlabel('Z [m]')
            plt.title(title+" Bz = f(z) ")
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


if __name__ == "__main__" :
    Bx=lambda x,y,z : np.cos(y*z)
    By = lambda x,y,z : np.sin(x)*np.cos(z)
    Bz = lambda x,y,z : np.sin(y*z)

    B=MagneticField(Bx=Bx,By=By,Bz=Bz)

    Z=np.linspace(0.0,2.*np.pi,10000)
    Y=np.linspace(0.0,2.*np.pi,10000)
    X=np.linspace(0.0,2.*np.pi,10000)

    print('plot X=pi/2 , Y=1 , Z in [0,2 pi]')
    B.plot_z(X=np.pi*0.5,Y=1.0,Z=Z)
    print('plot X= pi/2 , Y and Z in [0,2 pi]')
    B.plot_z(X=np.pi*0.5,Y=Y,Z=Z)
    print('plot X , Y and Z in [0,2 pi]')
    B.plot_z(X=X,Y=Y,Z=Z)