import numpy as np






class MagneticField(object):
    def __init__(self, x,y,z,Bx,By,Bz):
        self.x = x
        self.y = y
        self.z = z
        self.Bx = Bx
        self.By = By
        self.Bz = Bz

    def copy(self):
        return MagneticField(x=self.x.copy(), y=self.y.copy(), z=self.z.copy(),
                             Bx=self.Bx.copy(), By=self.By.copy(), Bz=self.Bz.copy())

    #faire une fct aui fait l'interpolation  et elargissement
    #

    # enelarge 2 vector for the interpolation of the magnetic field use in trajectory_undulator_from_magnetic_field2
    def enlargement_vector_for_interpolation(self, nb_enlarg=0):
        dz = self.z[1] - self.z[0]
        enlarg_1 = np.linspace(self.z[0] - nb_enlarg * dz, self.z[0] - dz, nb_enlarg)
        enlarg_2 = np.linspace(self.z[- 1] + dz, self.z[- 1] + nb_enlarg * dz, nb_enlarg)
        self.z = np.concatenate((enlarg_1, self.z))
        self.z = np.concatenate((self.z, enlarg_2))
        if (type(self.y)==np.ndarray) :
            dy = self.y[1] - self.y[0]
            enlarg_1 = np.linspace(self.y[0] - nb_enlarg * dz, self.y[0] - dy, nb_enlarg)
            enlarg_2 = np.linspace(self.y[- 1] + dy, self.y[- 1] + nb_enlarg * dy, nb_enlarg)
            self.y = np.concatenate((enlarg_1, self.y))
            self.y = np.concatenate((self.y, enlarg_2))
        if (type(self.By)==np.ndarray) :
            enlarg_3 = np.zeros(nb_enlarg)
            self.By = np.concatenate((enlarg_3, self.By))
            self.By = np.concatenate((self.By, enlarg_3))

