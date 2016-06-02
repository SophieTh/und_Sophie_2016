import numpy as np






class MagneticField(object):
    def __init__(self, x,y,z,Bx,By,Bz):
        self.x = x
        self.y = y
        self.z = z
        self.Bx = Bx
        self.By = By
        self.Bz = Bz

    #faire une fct aui fait l'interpolation  et elargissement
    #

    # enelarge 2 vector for the interpolation of the magnetic field use in trajectory_undulator_from_magnetic_field2
    def enlargement_vector_for_interpolation(self, nb_enlarg=0):
        Z=self.z
        By=self.By
        dz = Z[1] - Z[0]
        enlarg_1 = np.linspace(Z[0] - nb_enlarg * dz, Z[0] - dz, nb_enlarg)
        enlarg_2 = np.linspace(Z[len(Z) - 1] + dz, Z[len(Z) - 1] + nb_enlarg * dz, nb_enlarg)
        Z = np.concatenate((enlarg_1, Z))
        Z = np.concatenate((Z, enlarg_2))
        enlarg_3 = np.zeros(nb_enlarg)
        By = np.concatenate((enlarg_3, By))
        By = np.concatenate((By, enlarg_3))
        self.z=Z
        self.By=By

