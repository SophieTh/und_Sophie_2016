import numpy as np
import scipy.constants as codata
import scipy.integrate as integrate
from RadiationFactory import RadiationFactory , RADIATION_METHOD_APPROX_FARFIELD,RADIATION_METHOD_FARFIELD,\
                                            RADIATION_METHOD_NEAR_FIELD,RADIATION_METHOD_AUTOMATIC


class RadiationFactoryAnalitic(RadiationFactory):
    def __init__(self,method,formula,omega,Nb_pts):
        super(self.__class__, self).__init__(method=method,formula=formula,omega=omega,Nb_pts=Nb_pts)

# approximation of the energy radiated by an electron in a PLANE undulator
    # warning : trajectory is the trajectory difine like the function "undulator trajectory" before :
    def energy_radiated_approximation_and_farfield(self,trajectory, x, y, distance):


        # N = trajectory.shape[1]
        N = trajectory.nb_points()
        if distance == None:
            # in radian :
            n_chap = (lambda t :np.array([x, y, 1.0 - 0.5 * (x ** 2 + y ** 2)]))
        #in meters :
        else :
            R = np.sqrt(x ** 2 + y ** 2 + distance ** 2)
            n_chap = np.array([x, y, distance]) / R

        E = np.zeros((3,), dtype=np.complex)

        A1=(lambda t: (n_chap[0]*trajectory.a_x(t) +n_chap[1]*trajectory.a_y(t) +n_chap[2]*trajectory.a_z(t) ))
        A2=(lambda  t :(n_chap[0]*(n_chap[0]-trajectory.v_x(t))+n_chap[1]*(n_chap[1]-trajectory.v_y(t))
            +n_chap[2]*(n_chap[2]-trajectory.v_z(t))))
        Alpha2 = (lambda t :np.exp(
            0. + 1j * self.omega * (t - n_chap[0] * trajectory.x(t)
                                    -n_chap[1] * trajectory.y(t)- n_chap[2] * trajectory.z(t))))
        Alpha1 = (lambda t : (1.0 / (1.0 - n_chap[0] * trajectory.v_x(t)
                         -n_chap[1] * trajectory.v_y(t)- n_chap[2] * trajectory.v_z(t))) ** 2 )

        integrand0=(lambda t: (A1(t)*(n_chap[0]-trajectory.v_x(t))-A2(t)*trajectory.a_x(t)) * Alpha2(t)*Alpha1(t))
        integrand1 = (lambda t : (A1(t) * (n_chap[1] - trajectory.v_y(t)) - A2(t) * trajectory.a_y(t))
                                 * Alpha2(t)*Alpha1(t))
        integrand2 = (lambda t: (A1(t)*(n_chap[2]-trajectory.v_z(t))-A2(t)*trajectory.a_z(t)) * Alpha2(t)*Alpha1(t))
        # integrand=integrand0,integrand1,integrand2
        # for k in range(3):
        #     # E[k] = np.trapz(integrand[k], self.trajectory.t)
        #     E[k] = integrate.quad(integrand[k], trajectory.t[0],trajectory.t[-1])
        E[0] , err= integrate.quad(integrand0, trajectory.t[0],trajectory.t[-1],limit=2*N)
        E[1],err = integrate.quad(integrand1, trajectory.t[0], trajectory.t[-1],limit=2*N)
        E[2],err = integrate.quad(integrand2, trajectory.t[0], trajectory.t[-1],limit=2*N)
        return (np.abs(E[0]) ** 2 + np.abs(E[1]) ** 2 + np.abs(E[2]) ** 2)

    def energy_radiated_approximation_and_farfield2(self, trajectory, x, y, distance):
        # N = trajectory.shape[1]
        N = trajectory.nb_points()
        if distance == None:
            # in radian :
            n_chap = np.array([x, y, 1.0 - 0.5 * (x ** 2 + y ** 2)])
        # in meters :
        else:
            R = np.sqrt(x ** 2 + y ** 2 + distance ** 2)
            n_chap = np.array([x, y, distance]) / R

        E = np.zeros((3,), dtype=np.complex)
        integrand = np.zeros((3, N), dtype=np.complex)
        A1 = (n_chap[1] * trajectory.v_z - n_chap[2] * trajectory.v_y)
        A2 = (-n_chap[0] * trajectory.v_z + n_chap[2] * trajectory.v_x)
        A3 = (n_chap[0] * trajectory.v_y - n_chap[1] * trajectory.v_x)
        Alpha2 = np.exp(
            0. + 1j * self.omega * (trajectory.t - n_chap[0] * trajectory.x
                                    - n_chap[1] * trajectory.y - n_chap[2] * trajectory.z))


        integrand[0] += ( n_chap[1]*A3 - n_chap[2]*A2) * Alpha2
        integrand[1] += (- n_chap[0]*A3 + n_chap[2]*A1) * Alpha2
        integrand[2] += ( n_chap[0]*A2 - n_chap[1]*A1) * Alpha2

        for k in range(3):
            # E[k] = np.trapz(integrand[k], self.trajectory.t)
            E[k] = integrate.simps(integrand[k], trajectory.t)
        E *= self.omega * 1j

        Alpha3 = (np.exp(1j*self.omega*(trajectory.t[-1] - n_chap[2] * trajectory.z[-1]))-
                  np.exp(1j*self.omega*(trajectory.t[0] - n_chap[2] * trajectory.z[0])))
        terme_bord = np.full((3,), 0. + 1j * 0., dtype=np.complex)
        terme_bord[0] += n_chap[0] * n_chap[2] * trajectory.v_z[0]
        terme_bord[1] += n_chap[1] * n_chap[2] * trajectory.v_z[0]
        terme_bord[2] += trajectory.v_z[0] * (n_chap[2]**2 - 1.0)
        terme_bord *= Alpha3 / (1.0 - n_chap[2] * trajectory.v_z[0])
        E += terme_bord

        return (np.abs(E[0]) ** 2 + np.abs(E[1]) ** 2 + np.abs(E[2]) ** 2)

#
# ### marche pas !!!
#         # approximation of the energy radiated by an electron in a PLANE undulator
#         # warning : trajectory is the trajectory difine like the function "undulator trajectory" before :
#     def energy_radiated_approximation_and_farfield3(self, trajectory, x, y, distance):
#         # N = trajectory.shape[1]
#         N = trajectory.nb_points()
#         if distance == None:
#                 # in radian :
#             n_chap = np.array([x, y, 1.0 - 0.5 * (x ** 2 + y ** 2)])
#         # in meters :
#         else:
#             R = np.sqrt(x ** 2 + y ** 2 + distance ** 2)
#             n_chap = np.array([x, y, distance]) / R
#
#         E = np.zeros((3,), dtype=np.complex)
#         integrand = np.zeros((3, N), dtype=np.complex)
#         A1 = (n_chap[0] * trajectory.a_x + n_chap[1] * trajectory.a_y + n_chap[2] * trajectory.a_z)
#         Alpha2 = np.exp(
#                 0. + 1j * self.omega * (trajectory.t - n_chap[0] * trajectory.x
#                                         - n_chap[1] * trajectory.y - n_chap[2] * trajectory.z))
#         Alpha1 = (1.0 / (1.0 - n_chap[0] * trajectory.v_x
#                              - n_chap[1] * trajectory.v_y - n_chap[2] * trajectory.v_z))
#
#         integrand[0] += (A1 * (n_chap[0] - trajectory.v_x) * Alpha1 - trajectory.a_x) *Alpha1* Alpha2
#         integrand[1] += (A1 * (n_chap[1] - trajectory.v_y) * Alpha1 - trajectory.a_y) *Alpha1* Alpha2
#         integrand[2] += (A1 * (n_chap[2] - trajectory.v_z) * Alpha1 - trajectory.a_z) *Alpha1* Alpha2
#         for k in range(3):
#                 # E[k] = np.trapz(integrand[k], self.trajectory.t)
#             E[k] = integrate.simps(integrand[k], trajectory.t)
#         return (np.abs(E[0]) ** 2 + np.abs(E[1]) ** 2 + np.abs(E[2]) ** 2)
#
#     # exact equation for the energy radiated
#     # warning !!!!!!! far field approximation
#     def energy_radiated_farfield(self,trajectory,gamma, x, y,distance):
#         N = trajectory.nb_points()
#         n_chap = np.array([x - trajectory.x * codata.c, y - trajectory.y * codata.c, distance - trajectory.z * codata.c])
#         R = np.zeros(n_chap.shape[1])
#         for i in range(n_chap.shape[1]):
#             R[i] = np.linalg.norm(n_chap[:, i])
#             n_chap[:, i] /= R[i]
#
#         E = np.zeros((3,), dtype=np.complex)
#         integrand = np.zeros((3, N), dtype=np.complex)
#         A1 = (n_chap[0] * trajectory.a_x + n_chap[1] * trajectory.a_y + n_chap[2] * trajectory.a_z)
#         A2 = (n_chap[0] * (n_chap[0] - trajectory.v_x) + n_chap[1] * (n_chap[1] - trajectory.v_y)
#               + n_chap[2] * (n_chap[2] - trajectory.v_z))
#         Alpha2 = np.exp(
#             0. + 1j * self.omega * (trajectory.t - n_chap[0] * trajectory.x
#                                     - n_chap[1] * trajectory.y - n_chap[2] * trajectory.z))
#         Alpha1 = (1.0 / (1.0 - n_chap[0] * trajectory.v_x
#                          - n_chap[1] * trajectory.v_y - n_chap[2] * trajectory.v_z)) ** 2
#         integrand[0] += ((A1 * (n_chap[0] - trajectory.v_x) - A2 * trajectory.a_x)
#                          ) * Alpha2 * Alpha1
#         integrand[1] += ((A1 * (n_chap[1] - trajectory.v_y) - A2 * trajectory.a_y)
#                          ) * Alpha2 * Alpha1
#         integrand[2] += ((A1 * (n_chap[2] - trajectory.v_z) - A2 * trajectory.a_z)
#                          ) * Alpha2 * Alpha1
#         for k in range(3):
#             # E[k] = np.trapz(integrand[k], trajectory[0])
#             E[k] = integrate.simps(integrand[k], trajectory.t)
#         return (np.abs(E[0]) ** 2 + np.abs(E[1]) ** 2 + np.abs(E[2]) ** 2)
#
#     # pb ne marche pas
#     def energy_radiated_farfield2(self, trajectory, gamma, x, y, distance):
#         N = trajectory.nb_points()
#         n_chap = np.array([x - trajectory.x * codata.c, y - trajectory.y * codata.c, distance - trajectory.z * codata.c])
#         R = np.zeros(n_chap.shape[1])
#         for i in range(n_chap.shape[1]):
#             R[i] = np.linalg.norm(n_chap[:, i])
#             n_chap[:, i] /= R[i]
#
#         E = np.zeros((3,), dtype=np.complex)
#         integrand = np.zeros((3, N), dtype=np.complex)
#
#         A1 = (n_chap[1] * trajectory.v_z - n_chap[2] * trajectory.v_y)
#         A2 = (-n_chap[0] * trajectory.v_z + n_chap[2] * trajectory.v_x)
#         A3 = (n_chap[0] * trajectory.v_y - n_chap[1] * trajectory.v_x)
#
#         Alpha1 = np.exp(
#             0. + 1j * self.omega * (trajectory.t +R))
#         integrand[0] += (n_chap[1] * A3 - n_chap[2] * A2    ) * Alpha1
#         integrand[1] += (-n_chap[0] * A3 + n_chap[2] * A1    ) * Alpha1
#         integrand[2] += (n_chap[0] * A2 - n_chap[1] * A1   ) * Alpha1
#         for k in range(3):
#             # E[k] = np.trapz(integrand[k], self.trajectory.t)
#             E[k] = integrate.simps(integrand[k], trajectory.t)
#         E *= self.omega * 1j
#
#         # Alpha3 = (np.exp(1j * self.omega * (trajectory.t[-1] - n_chap[2][-1] * trajectory.z[-1])) -
#         #           np.exp(1j * self.omega * (trajectory.t[0] - n_chap[2][0] * trajectory.z[0])))
#         # terme_bord = np.full((3,), 0. + 1j * 0., dtype=np.complex)
#         # terme_bord[0] += (n_chap[0][0] * n_chap[2][0] * trajectory.v_z[0]
#         #                     )*Alpha3 / (1.0 - n_chap[2][0] * trajectory.v_z[0])
#         # terme_bord[1] += (n_chap[1][0] * n_chap[2][0] * trajectory.v_z[0]
#         #             )*Alpha3 / (1.0 - n_chap[2][0] * trajectory.v_z[0])
#         # terme_bord[2] += (trajectory.v_z[0] * (n_chap[2][0] ** 2 - 1.0)
#         #                 )*Alpha3 / (1.0 - n_chap[2][0] * trajectory.v_z[0])
#         # E += terme_bord
#
#         return (np.abs(E[0]) ** 2 + np.abs(E[1]) ** 2 + np.abs(E[2]) ** 2)
#
#     # energy radiated without the the far filed approxiamation
#     def energy_radiated_near_field(self,gamma,trajectory, x, y,distance):
#         N = trajectory.nb_points()
#         n_chap = np.array([x - trajectory.x * codata.c, y - trajectory.y * codata.c, distance - trajectory.z * codata.c])
#         R = np.zeros(n_chap.shape[1])
#         for i in range(n_chap.shape[1]):
#             R[i] = np.linalg.norm(n_chap[:, i])
#             n_chap[:, i] /= R[i]
#
#         E = np.zeros((3,), dtype=np.complex)
#         integrand = np.zeros((3, N), dtype=np.complex)
#         A1 = (n_chap[0] * trajectory.a_x + n_chap[1] * trajectory.a_y + n_chap[2] * trajectory.a_z)
#         A2 = (n_chap[0] * (n_chap[0] - trajectory.v_x) + n_chap[1] * (n_chap[1] - trajectory.v_y)
#               + n_chap[2] * (n_chap[2] - trajectory.v_z))
#         Alpha2 = np.exp(
#             0. + 1j * self.omega * (trajectory.t - n_chap[0] * trajectory.x
#                                     -n_chap[1] * trajectory.y- n_chap[2] * trajectory.z))
#         Alpha1 = (1.0 / (1.0 - n_chap[0] * trajectory.v_x
#                          -n_chap[1] * trajectory.v_y- n_chap[2] * trajectory.v_z)) ** 2
#         Alpha3 = codata.c / (gamma ** 2 * R)
#         integrand[0] += ((A1 * (n_chap[0] - trajectory.v_x) - A2 * trajectory.a_x)
#                          + Alpha3 * (n_chap[0] - trajectory.v_x)
#                          ) * Alpha2 * Alpha1
#         integrand[1] += ((A1 * (n_chap[1] - trajectory.v_y) - A2 * trajectory.a_y)
#                          + Alpha3 * (n_chap[1] - trajectory.v_y)
#                          ) * Alpha2 * Alpha1
#         integrand[2] += ((A1 * (n_chap[2] - trajectory.v_z) - A2 * trajectory.a_z)
#                          + Alpha3 * (n_chap[2] - trajectory.v_z)
#                          ) * Alpha2 * Alpha1
#         for k in range(3):
#             # E[k] = np.trapz(integrand[k], trajectory[0])
#             E[k] = integrate.simps(integrand[k], trajectory.t)
#         return (np.abs(E[0]) ** 2 + np.abs(E[1]) ** 2 + np.abs(E[2]) ** 2)
#
#     # energy radiated without the the far filed approxiamation
#     # faux pour le mmt,
#     def energy_radiated_near_field3(self, gamma, trajectory, x, y, distance):
#         N = trajectory.nb_points()
#         n_chap = np.array(
#             [x - trajectory.x * codata.c, y - trajectory.y * codata.c, distance - trajectory.z * codata.c])
#         R = np.zeros(n_chap.shape[1])
#         for i in range(n_chap.shape[1]):
#             R[i] = np.linalg.norm(n_chap[:, i])
#             n_chap[:, i] /= R[i]
#
#         E = np.zeros((3,), dtype=np.complex)
#         integrand = np.zeros((3, N), dtype=np.complex)
#         A1 = (n_chap[0] * trajectory.a_x + n_chap[1] * trajectory.a_y + n_chap[2] * trajectory.a_z)
#         A2 = (n_chap[0] * (n_chap[0] - trajectory.v_x) + n_chap[1] * (n_chap[1] - trajectory.v_y)
#               + n_chap[2] * (n_chap[2] - trajectory.v_z))
#         Alpha2 = np.exp(
#             0. + 1j * self.omega * (trajectory.t - n_chap[0] * trajectory.x
#                                     - n_chap[1] * trajectory.y - n_chap[2] * trajectory.z))
#         Alpha1 = (1.0 / (1.0 - n_chap[0] * trajectory.v_x
#                          - n_chap[1] * trajectory.v_y - n_chap[2] * trajectory.v_z)) ** 2
#         Alpha3 = codata.c / (gamma ** 2 * R)
#         integrand[0] += ((A1 * (n_chap[0] - trajectory.v_x) - A2 * trajectory.a_x)
#                          + Alpha3 * (n_chap[0] - trajectory.v_x)
#                          ) * Alpha2 * Alpha1
#         integrand[1] += ((A1 * (n_chap[1] - trajectory.v_y) - A2 * trajectory.a_y)
#                          + Alpha3 * (n_chap[1] - trajectory.v_y)
#                          ) * Alpha2 * Alpha1
#         integrand[2] += ((A1 * (n_chap[2] - trajectory.v_z) - A2 * trajectory.a_z)
#                          + Alpha3 * (n_chap[2] - trajectory.v_z)
#                          ) * Alpha2 * Alpha1
#         for k in range(3):
#             # E[k] = np.trapz(integrand[k], trajectory[0])
#             E[k] = np.trapz(integrand[k], trajectory.t)
#         return (np.abs(E[0]) ** 2 + np.abs(E[1]) ** 2 + np.abs(E[2]) ** 2)
#
#
#     # attention ici c'est faux
#     # multiplier par 4 ???
#     def energy_radiated_near_field2(self, trajectory, gamma,x, y, distance):
#         # N = trajectory.shape[1]
#         N = trajectory.nb_points()
#         n_chap = np.array(
#             [x - trajectory.x * codata.c, y - trajectory.y * codata.c, distance - trajectory.z * codata.c])
#         R = np.zeros(n_chap.shape[1])
#         for i in range(n_chap.shape[1]):
#             R[i] = np.linalg.norm(n_chap[:, i])
#             n_chap[:, i] /= R[i]
#
#         E = np.zeros((3,), dtype=np.complex)
#         integrand = np.zeros((3, N), dtype=np.complex)
#
#         A1 = (n_chap[1] * trajectory.v_z - n_chap[2] * trajectory.v_y)
#         A2 = (-n_chap[0] * trajectory.v_z + n_chap[2] * trajectory.v_x)
#         A3 = (n_chap[0] * trajectory.v_y - n_chap[1] * trajectory.v_x)
#         Alpha1 = np.exp(
#             0. + 1j * self.omega * (trajectory.t - n_chap[0] * trajectory.x
#                                     - n_chap[1] * trajectory.y - n_chap[2] * trajectory.z))
#         Alpha2 = codata.c / (gamma ** 2 * R)
#         integrand[0] += (n_chap[1] * A3 - n_chap[2] * A2
#                          + Alpha2 * (n_chap[0] - trajectory.v_x)
#                          ) * Alpha1
#         integrand[1] += (-n_chap[0] * A3 + n_chap[2] * A1
#                          + Alpha2 * (n_chap[1] - trajectory.v_y)
#                          ) * Alpha1
#         integrand[2] += (n_chap[0] * A2 - n_chap[1] * A1
#                          + Alpha2 * (n_chap[2] - trajectory.v_z)
#                          ) * Alpha1
#         for k in range(3):
#             # E[k] = np.trapz(integrand[k], self.trajectory.t)
#             E[k] = integrate.simps(integrand[k], trajectory.t)
#         E *= self.omega * 1j
#
#         # terme_bord = np.full((3), 0. + 1j * 0., dtype=np.complex)
#         # unit_z = np.array([0.0, 0.0, 1.0])
#         # terme_bord += ((n_chap[2][-1] * n_chap[:, -1] - unit_z) *
#         #                np.exp(1j * self.omega * (trajectory.t[-1] - n_chap[2][-1] * trajectory.z[-1])) / (
#         #                    1.0 - n_chap[2][-1] * trajectory.v_z[-1]
#         #                ))
#         # terme_bord -= ((n_chap[2][0] * n_chap[:, 0] - unit_z) *
#         #                np.exp(1j * self.omega * (trajectory.t[-1] - n_chap[2][0] * trajectory.z[0])) / (
#         #                    1.0 - n_chap[2][0] * trajectory.v_z[0]
#         #                ))
#         # terme_bord *= trajectory.v_z[0]
#         # E += terme_bord
#
#         return (np.abs(E[0]) ** 2 + np.abs(E[1]) ** 2 + np.abs(E[2]) ** 2)
