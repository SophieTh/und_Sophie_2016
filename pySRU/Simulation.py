import numpy as np
import matplotlib.pyplot as plt
from  scipy.constants import physical_constants
from mpl_toolkits.mplot3d import Axes3D
import scipy.constants as codata
import time
from scipy.interpolate import interp1d
from pySRU.Radiation import Radiation
from pySRU.MagneticField import MagneticField
from pySRU.SourceUndulatorPlane import SourceUndulatorPlane
from pySRU.SourceBendingmagnet import SourceBendingMagnet
from pySRU.MagneticStructureUndulatorPlane import MagneticStructureUndulatorPlane as Undulator
from pySRU.MagneticStructureBendingMagnet import MagneticStructureBendingMagnet as BM

from pySRU.TrajectoryFactory import TrajectoryFactory, TRAJECTORY_METHOD_ANALYTIC,TRAJECTORY_METHOD_ODE

from pySRU.RadiationFactory import RadiationFactory,RADIATION_METHOD_NEAR_FIELD, \
                                 RADIATION_METHOD_APPROX_FARFIELD


eV_to_J = codata.e

class Simulation(object):

    def __init__(self, source, trajectory_fact, radiation_fact, trajectory, radiation):
        self.source=source
        self.trajectory_fact=trajectory_fact
        self.radiation_fact=radiation_fact
        self.trajectory=trajectory
        self.radiation=radiation

# faire un autre constructeur pour la copie


    def copy(self):
        return Simulation(source=self.source.copy(), trajectory_fact=self.trajectory_fact.copy(),
                                   radiation_fact=self.radiation_fact.copy(),
                                   trajectory=self.trajectory.copy(), radiation=self.radiation.copy())


    def update_radiation(self):
        self.radiation.intensity = self.radiation_fact.calculate_radiation_intensity(trajectory=self.trajectory,
                                                                                     source=self.source,
                                                                                     distance=self.radiation.distance,
                                                                                     X_array=self.radiation.X,
                                                                                     Y_array=self.radiation.Y)
    #change

    def change_distance(self,D,update_radiation=1):
        self.radiation.distance=D
        if update_radiation:
            self.update_radiation()

    def change_radiation_method(self,method,update_radiation=1):
        self.radiation_fact.method=method
        if update_radiation:
            self.update_radiation()

    def change_trajectory_method(self,method,update_radiation=1) :
        self.trajectory_fact.method = method
        self.trajectory =self.trajectory_fact.create_from_source(source=self.source)
        if update_radiation:
            self.update_radiation()

    def change_initial_condition(self, initial_cond,update_radiation=1):
        #TODO change this
        self.trajectory_fact.initial_condition= initial_cond
        self.magnetic_field.change_Zo(Zo=initial_cond[5])
        self.trajectory = self.trajectory_fact.create_for_parameter(parameter=self.parameter,B=self.magnetic_field)
        if update_radiation:
            self.update_radiation()

    def change_omega(self,omega,update_radiation=1) :
        self.radiation_fact.omega=omega
        if update_radiation:
            self.update_radiation()

    def change_harmonic_number(self,harmonic_number,update_radiation=1) :
        self.radiation_fact.omega=harmonic_number*self.source.harmonic_frequency()
        if update_radiation:
            self.update_radiation()

    def change_energy_eV(self,E,update_radiation=1) :
        omega = E * eV_to_J / codata.hbar
        self.change_omega(omega)
        if update_radiation:
            self.update_radiation()

    def change_Nb_pts_trajectory(self,Nb_pts,update_radiation=1):
        self.trajectory_fact.Nb_pts = Nb_pts
        self.trajectory = self.trajectory_fact.create_from_source(source=self.source)
        if update_radiation:
            self.update_radiation()

    def change_nb_period(self,Nb_period,update_radiation=1):
        self.source.magnetic_structure.length=Nb_period*self.source.magnetic_structure.period_length
        self.source.magnetic_field=self.source.magnetic_structure.create_magnetic_field()
        self.trajectory_fact.Nb_pts = int(self.source.choose_nb_pts_trajectory(2))
        self.trajectory = self.trajectory_fact.create_from_source(source=self.source)
        if update_radiation:
            self.update_radiation()

    def change_Nb_pts_radiation(self,Nb_pts,update_radiation=1):
        self.radiation_fact.Nb_pts = Nb_pts
        self.radiation.change_Nb_pts(Nb_pts)
        if update_radiation:
            self.update_radiation()

    def change_XY_radiation(self,X,Y,update_radiation=1):

        self.radiation.X=X
        self.radiation.Y=Y
        if update_radiation:
            self.update_radiation()

    #TODO possible que pour Undulator, a revoir , a tester
    def calculate_on_ring_number(self, ring_number_max,update_radiation=1):
        harmonic_number=np.floor(self.radiation_fact.omega/self.source.harmonic_frequency(1))
        self.radiation_fact.omega=self.source.harmonic_frequency(harmonic_number)
        X=np.array([])
        Y = np.array([])
        t = np.linspace(0.0, 2.0* np.pi, self.radiation_fact.Nb_pts)
        for i in range(int(ring_number_max)+1) :
            Xi,Yi= self.source.describe_ring(distance=self.radiation.distance,harmonic_number=harmonic_number,
                                              ring_number=i,t=t)
            X=np.concatenate((X,Xi))
            Y=np.concatenate((Y,Yi))
        self.change_XY_radiation(X=X,Y=Y,update_radiation=update_radiation)


    def calculate_on_central_cone(self,harmonic_number=1,npoints_x=None,npoints_y=None,update_radiation=1):
        theta_max=self.source.angle_deflection_central_cone(harmonic_number)
        if self.radiation.distance==None :
            X_max=theta_max
            Y_max=theta_max
        else :
            X_max=self.radiation.distance*theta_max
            Y_max=self.radiation.distance*theta_max

        if npoints_x is None:
            npoints_x = self.radiation_fact.Nb_pts
        if npoints_y is None:
            npoints_y = self.radiation_fact.Nb_pts

        X=np.linspace(0.0,X_max,npoints_x)
        Y=np.linspace(0.0,Y_max,npoints_y)
        X,Y=np.meshgrid(X,Y)
        self.change_harmonic_number(harmonic_number,update_radiation=0)
        self.change_XY_radiation(X=X,Y=Y,update_radiation=update_radiation)

    def calculate_until_ring_number(self,ring_number,harmonic_number=1, XY_are_list=True,update_radiation=1,
                                    npoints=None):
        # harm_num=self.radiation_fact.omega/self.source.harmonic_frequency(1)
        if npoints != None:
            self.radiation_fact.Nb_pts = npoints

        self.change_harmonic_number(harmonic_number,update_radiation=0)
        observation_angle = np.linspace(0.0, self.source.angle_ring_number(harmonic_number=harmonic_number,
                                        ring_number=ring_number),self.radiation_fact.Nb_pts)
        self.calculate_for_observation_angles(observation_angle,XY_are_list=XY_are_list,update_radiation=update_radiation)



    def calculate_for_observation_angles(self,observation_angle,XY_are_list=True,update_radiation=1):
        D=self.radiation.distance
        if D==None :
            D=1
        if XY_are_list :
            X = np.array([])
            Y = np.array([])
            t = np.linspace(0.0, 2.0 * np.pi, self.radiation_fact.Nb_pts)
            for theta in observation_angle:
                if theta==0.0 :
                    Xi=np.array([0.0])
                    Yi=np.array([0.0])
                else :
                    Xi=np.cos(t)*D*theta
                    Yi=np.sin(t)*D*theta
                X = np.concatenate((X, Xi))
                Y = np.concatenate((Y, Yi))
        else :
            X=D*observation_angle
            Y=D*observation_angle
            X,Y=np.meshgrid(X,Y)
        self.change_XY_radiation(X=X,Y=Y,update_radiation=update_radiation)


    #######
    # spectrum
    #######

    def spectrum(self, omega_array=None, is_quadrant=0):
        if omega_array is None:
            omega1 = self.source. choose_photon_frequency()
            omega_array1 = np.arange(omega1 * 0.8, omega1 * 1.11, omega1 * 0.02)
            omega_array3 = np.arange(omega1 * 1.8, omega1 * 2.11, omega1 * 0.02)
            omega_array5 = np.arange(omega1 * 2.8, omega1 * 3.11, omega1 * 0.02)

            omega_array2=np.arange(omega1*1.12,omega1*1.79,omega1*0.4)
            omega_array4=np.arange(omega1*2.12,omega1*2.79,omega1*0.4)

            omega_array = np.concatenate((omega_array1,omega_array2,omega_array3,omega_array4,omega_array5))

        print("Calculating spectrum for %d energy points..."%(omega_array.size))
        spectrum = np.zeros_like(omega_array)

        for i in range(len(spectrum)):
            print("Point %d of %d..."%(1+i,spectrum.size))
            self.change_omega(omega_array[i])
            spectrum[i] = self.radiation.integration(is_quadrant=is_quadrant)

        return spectrum, omega_array

    def spectrum_on_radiation_max(self, omega_array=None):
        if omega_array == None:
            omega1 = self.source.omega1()
            omega_array = np.arange(omega1 * 0.9, 5.0 * omega1 * 1.1, omega1 * 0.01)
        spectrum = np.zeros_like(omega_array)
        print(len(spectrum))
        for i in range(len(spectrum)):
            print(i)
            self.change_omega(omega_array[i])
            spectrum[i] = self.radiation.max()
        return spectrum, omega_array

    def spectrum_on_axis(self, omega_array=None):
        save_X=self.radiation.X
        X=np.array([0.0])
        save_Y=self.radiation.Y
        Y=np.array([0.0])
        self.change_XY_radiation(X=X,Y=Y)
        if omega_array == None:
            omega1 = self.source.harmonic_frequency(1)
            omega_array = np.arange(omega1 * 0.9, 5.0 * omega1 * 1.1, omega1 * 0.001)
        spectrum = np.zeros_like(omega_array)
        for i in range(len(spectrum)):
            self.change_omega(omega_array[i])
            spectrum[i] = self.radiation.intensity.max()
        self.change_XY_radiation(X=save_X,Y=save_Y)
        return spectrum, omega_array


    def calculate_spectrum_on_slit(self,abscissas_array=None,use_eV=0,is_quadrant=0,do_plot=1):

        if use_eV:
            conversion_factor = codata.hbar / eV_to_J
            xlabel = "Photon energy [eV]"
        else:
            conversion_factor = 1.0
            xlabel = "Photon frequency [Hz]"

        if abscissas_array is None:
            omega_array = None
        else:
            omega_array = abscissas_array / conversion_factor

        spectrum,omega_array=self.spectrum(omega_array=omega_array,is_quadrant=is_quadrant)

        if do_plot:
            import matplotlib.pyplot as plt
            plt.plot(conversion_factor*omega_array, spectrum, label='calculad radiation on slit')
            plt.legend()
            plt.xlabel(xlabel)
            plt.ylabel("Flux [phot/s/0.1%bw]")
            plt.show()

        print("Maximum flux is %g"%spectrum.max())
        return conversion_factor*omega_array,spectrum



    #TODO possible que pour undultor et c'est le quart du cone central !
    def calculate_spectrum_central_cone(self,harmonic_number=1,theoretical_value=False,abscissas_array=None,use_eV=0,do_plot=1):


        if use_eV:
            conversion_factor = codata.hbar / eV_to_J
            xlabel = "Photon energy [eV]"
        else:
            conversion_factor = 1.0
            xlabel = "Photon frequency [Hz]"

        if abscissas_array is None:
            omega_array = None
        else:
            omega_array = abscissas_array / conversion_factor

        #self.radiation_fact.distance=None
        theta_max = self.source.angle_deflection_central_cone()
        if self.radiation.distance == None:
            X_max = theta_max
            Y_max = theta_max
        else:
            X_max = self.radiation.distance * theta_max
            Y_max = self.radiation.distance * theta_max

        X = np.linspace(0.0, X_max, 21)
        Y = np.linspace(0.0, Y_max, 21)
        XY = np.meshgrid(X,Y)

        self.change_XY_radiation(XY[0],XY[1],update_radiation=0)
        # if do_plot:
        #     self.radiation.plot()


        # main calculation in one quadrant
        if omega_array is None: # pre-calculate default
            omega1 = self.source.harmonic_frequency(1)
            omega_array = np.linspace(omega1*harmonic_number*0.95,omega1*harmonic_number*1.05,30)

        spectrum,omega_array=self.spectrum(omega_array=omega_array,is_quadrant=1)

        if do_plot:
            import matplotlib.pyplot as plt

            if theoretical_value :
                spectrum_theoretical=self.source.theoretical_flux_integrated_central_cone(harmonic_number)
                plt.plot(conversion_factor*harmonic_number*omega1, spectrum_theoretical, "g*", label='theoretical value')
                print("Theoretical flux",spectrum_theoretical,"Max spectrum: ",spectrum.max())
            plt.plot(conversion_factor*omega_array, spectrum, label='calculad radiation on central cone for harmonic %d'%harmonic_number)
            plt.legend()
            plt.xlabel(xlabel)
            plt.ylabel("Flux [phot/s/0.1%bw]")
            plt.show()

        return conversion_factor*omega_array,spectrum




    def calculate_spectrum_on_axis(self, omega_array=None, title="On axis flux", use_eV=0, do_plot=1):

        if use_eV:
            conversion_factor = codata.hbar / eV_to_J
            xlabel = "Photon energy [eV]"
        else:
            conversion_factor = 1.0
            xlabel = "Photon frequency [Hz]"

        spectrum,omega_array=self.spectrum_on_axis(omega_array=omega_array)

        if do_plot:

            # compute theoretical spectrum (to compare with)
            omega1 = self.source.harmonic_frequency(1)
            harm_num_min = np.ceil(omega_array.min() / omega1)
            harm_num_max = np.floor(omega_array.max() / omega1)
            harm_number=np.arange(harm_num_min,harm_num_max+0.5,1)

            spectrum_theoretical_scan = self.source.theoretical_flux_on_axis_omega_scan(harm_number,omega_array)
            # not needed
            spectrum_theoretical=np.zeros_like(harm_number)
            for i in range(len(spectrum_theoretical)):
                spectrum_theoretical[i]=self.source.theoretical_flux_on_axis(harm_number[i])

            # plot
            import matplotlib.pyplot as plt
            plt.plot(conversion_factor*omega_array,spectrum,label='calculated radiation on axis')
            plt.plot(conversion_factor*omega_array,spectrum_theoretical_scan,label='theoretical radiation on axis')
            plt.plot(conversion_factor*harm_number*omega1,spectrum_theoretical, "g*",label='theoretical peak value')
            plt.legend()
            plt.xlabel(xlabel)
            plt.ylabel("Flux (phot/s/0.1%bw/mrad2")
            plt.title(title)
            plt.show()

        return conversion_factor*omega_array,spectrum
    # plot

    def plot_magnetic_field_along_Z(self,X=None,Y=None):
        Z=self.trajectory.t*self.source.electron_speed()*codata.c
        if Y==None :
            Y=self.trajectory_fact.initial_condition[4]
        if X == None:
            X = self.trajectory_fact.initial_condition[3]
        self.source.magnetic_field.plot_z(X=X,Y=Y,Z=Z)

    def print_parameters(self):
        self.source.print_parameters()

        self.trajectory_fact.print_parameters()

        self.radiation_fact.print_parameters()

        print(" Others")
        print("     Distance: ",self.radiation.distance)
        print("     Central Cone width (urad) for  n=1: %6.3f, n=3: %6.3f, n=5: %6.3f "%
              ( 1e6 * self.source.angle_deflection_central_cone(harmonic_number=1),
                1e6 * self.source.angle_deflection_central_cone(harmonic_number=3),
                1e6 * self.source.angle_deflection_central_cone(harmonic_number=5) ))
        print("     First ring position (urad) for n=1: %6.3f, n=3: %6.3f, n=5: %6.3f "%
              ( 1e6 * self.source.angle_ring_number(harmonic_number=1,ring_number=1),
                1e6 * self.source.angle_ring_number(harmonic_number=3,ring_number=1),
                1e6 * self.source.angle_ring_number(harmonic_number=5,ring_number=1) ))

    def plot_everything(self):
        self.trajectory.plot()
        self.trajectory.plot_3D()
        self.radiation.plot()
        self.calculate_on_central_cone()
        self.radiation.plot()
        self.plot_spectrum()
        self.calculate_spectrum_on_axis()



    #TODO a tester attention possible que pour BENDING magnet
    def create_theoretical_radiation(self):
        #TODO faire le cas None
        radiation=self.radiation.copy()
        shape1 = radiation.intensity.shape
        # X = radiation.X.flatten()
        # Y = radiation.X.flatten()
        # res=radiation.intensity.flatten()
        # cas grid :
        print(type(radiation))
        print(type(radiation.Y))
        Y=radiation.Y
        X=radiation.X
        for i in range(Y.shape[0]):
            for j in range(Y.shape[1]) :
                observation_angle=(Y[i,j]/radiation.distance)
                #TODO changer pour qu'on calcul avec des angles
                radiation.intensity[i,j]=self.source.radiation_theoretical(
                    omega=self.radiation_fact.omega,observation_angle=observation_angle)
        return radiation

    def create_theoretical_radiation2(self):
        radiation = self.radiation.copy()
        shape1 = radiation.intensity.shape
        # X = radiation.X.flatten()
        # Y = radiation.X.flatten()
        # res=radiation.intensity.flatten()
        # cas grid :
        print(type(radiation))
        print(type(radiation.Y))
        Y = radiation.Y
        X = radiation.X
        for i in range(Y.shape[0]):
            for j in range(Y.shape[1]):
                observation_angle = (Y[i, j] / radiation.distance)
                radiation.intensity[i, j] = self.source.radiation_theoretical(
                    omega=self.radiation_fact.omega, observation_angle=observation_angle)
        return radiation

#

def create_simulation(magnetic_structure,electron_beam, magnetic_field=None, photon_energy=None,
                      traj_method=TRAJECTORY_METHOD_ANALYTIC,Nb_pts_trajectory=None,
                      rad_method=RADIATION_METHOD_APPROX_FARFIELD, Nb_pts_radiation=101,
                      initial_condition=None, distance=None,XY_are_list=False,X=None,Y=None) :

    if type(magnetic_structure)==Undulator :
        source=SourceUndulatorPlane(undulator=magnetic_structure,
                                    electron_beam=electron_beam, magnetic_field=magnetic_field)
        print("Calculating undulator source...")
    elif type(magnetic_structure)==BM:
        source = SourceBendingMagnet(magnetic_structure=magnetic_structure,
                                      electron_beam=electron_beam, magnetic_field=magnetic_field)
        print("Calculating bending magnet source...")
    else :
        raise Exception('magnet type unknown')

    if photon_energy==None :
        omega=source.choose_photon_frequency(harmonic_number=1)
    else :
        omega = photon_energy * eV_to_J / codata.hbar

    if Nb_pts_trajectory==None :
        Nb_pts_trajectory = int(source.choose_nb_pts_trajectory(2))

    if distance==None and (rad_method==RADIATION_METHOD_NEAR_FIELD) :
        distance=source.choose_distance_automatic(2)

    if X is None or Y is None :
        if (X != None) :
            Y=X
        elif Y != None :
            X=Y
        else :
            theta_max=source.choose_angle_deflection_max()
            if distance==None :
                X_max=theta_max
                Y_max=theta_max
            else :
                X_max = distance * theta_max
                Y_max = distance * theta_max
            X = np.linspace(0.0, X_max, Nb_pts_radiation)
            Y = np.linspace(0.0, Y_max, Nb_pts_radiation)


    if type(X) == float:
        X= np.linspace(0.0, X, Nb_pts_radiation)
    if type(Y) == float:
        Y = np.linspace(0.0, Y, Nb_pts_radiation)

    if X.shape != Y.shape :
        raise Exception('X and Y must have the same shape')
    Nb_pts_radiation=len(X.flatten())

    #print('step 1')
    traj_fact=TrajectoryFactory(Nb_pts=Nb_pts_trajectory,method=traj_method,initial_condition=initial_condition)
    if (traj_fact.initial_condition == None):
        traj_fact.initial_condition = source.choose_initial_contidion_automatic()


    #print('step 2')
    rad_fact=RadiationFactory(method=rad_method,omega=omega,Nb_pts=Nb_pts_radiation)

    #print('step 3')
    trajectory=traj_fact.create_from_source(source=source)
    #print('step 4')
    radiation = rad_fact.create_for_one_relativistic_electron(trajectory=trajectory, source=source, XY_are_list=XY_are_list,
                                                              distance=distance, X=X, Y=Y)

    #print('step 5')
    return Simulation(source=source, trajectory_fact=traj_fact,
                               radiation_fact=rad_fact, trajectory=trajectory, radiation=radiation)


# exemples

def Example_minimum():

    from  pySRU.ElectronBeam import ElectronBeam

    print("======================================================================")
    print("======      Undulator from X-ray data booklet                  =======")
    print("====== fig 2.5 in  http://xdb.lbl.gov/Section2/Sec_2-1.html    =======")
    print("======================================================================")

    # note that the flux in the reference fig 2.6 is a factor 10 smaller than the calculated here.
    # This factor comes from the units:
    #     here: phot / s  / A / 0.1%bw / (mrad)^2
    #     ref : phot / s  / A /   1%bw / (0.1 mrad)^2

    undulator_test = Undulator(K=1.87, period_length=0.035, length=0.035 * 14)
    electron_beam_test = ElectronBeam(Electron_energy=1.3, I_current=1.0)

    simulation_test = create_simulation(magnetic_structure=undulator_test,electron_beam=electron_beam_test,
                        magnetic_field=None, photon_energy=None,
                        traj_method=TRAJECTORY_METHOD_ANALYTIC,Nb_pts_trajectory=None,
                        rad_method=RADIATION_METHOD_APPROX_FARFIELD, Nb_pts_radiation=101,
                        initial_condition=None, distance=None,XY_are_list=False,X=None,Y=None)


    simulation_test.print_parameters()

    simulation_test.trajectory.plot_3D(title="Electron Trajectory")

    simulation_test.radiation.plot(title="Flux in far field vs angle")

    #  second calculation, change distance and also change the grid to accept the central cone
    simulation_test.change_distance(D=100)
    simulation_test.calculate_on_central_cone()

    simulation_test.radiation.plot(title="New distance: %3.1f m"%simulation_test.radiation.distance)



def Example_meshgrid():
    from pySRU.ElectronBeam import ElectronBeam

    print("======================================================================")
    print("======      Undulator U18 from ESRF with K=1.68                =======")
    print("======================================================================")

    beam_ESRF = ElectronBeam(Electron_energy=6.0, I_current=0.2)
    ESRF18 = Undulator(K=1.68, period_length=0.018, length=2.0)

    #
    # radiation in a defined meah
    #
    print('create radiation a given screen (40mm x 40mm @ 100 m )')
    X = np.linspace(-0.02,0.02,101)
    Y = np.linspace(-0.02,0.02,101)
    distance = 100.0
    simulation_test = create_simulation(magnetic_structure=ESRF18,electron_beam=beam_ESRF,
                        magnetic_field=None, photon_energy=None,
                        traj_method=TRAJECTORY_METHOD_ANALYTIC,Nb_pts_trajectory=None,
                        rad_method=RADIATION_METHOD_APPROX_FARFIELD, Nb_pts_radiation=101,
                        initial_condition=None, distance=distance,XY_are_list=False,X=X,Y=Y)

    simulation_test.print_parameters()
    simulation_test.trajectory.plot_2D()
    simulation_test.radiation.plot(title=" radiation in a defined screen (40mm x 40mm @ 100 m )")

    #
    # up to a maximum X and Y
    #
    print('create simulation for a given maximum X and Y ')
    simulation_test = create_simulation(magnetic_structure=ESRF18, electron_beam=beam_ESRF,
                                        traj_method=TRAJECTORY_METHOD_ANALYTIC,
                                        rad_method=RADIATION_METHOD_APPROX_FARFIELD,
                                        distance=100,X=0.01,Y=0.01)

    simulation_test.radiation.plot(title='simulation for a maximum X=0.01 and Y=0.01')


def Example_meshgrid_on_central_cone_and_rings():
    from pySRU.ElectronBeam import ElectronBeam

    print("======================================================================")
    print("======      Undulator U18 from ESRF with K=1.68                =======")
    print("======================================================================")

    beam_ESRF = ElectronBeam(Electron_energy=6.0, I_current=0.2)
    ESRF18 = Undulator(K=1.68, period_length=0.018, length=2.0)

    #
    # radiation in a defined mesh with only one point to save time
    #

    distance = None
    simulation_test = create_simulation(magnetic_structure=ESRF18,electron_beam=beam_ESRF,
                        magnetic_field=None, photon_energy=None,
                        traj_method=TRAJECTORY_METHOD_ANALYTIC,Nb_pts_trajectory=None,
                        rad_method=RADIATION_METHOD_APPROX_FARFIELD, Nb_pts_radiation=101,
                        initial_condition=None, distance=distance,XY_are_list=False,X=np.array([0]),Y=np.array([0]))


    #
    # central cone
    #
    harmonic_number = 1
    print('create radiation in a screen including central cone for harmonic %d'%harmonic_number)
    simulation_test.calculate_on_central_cone(harmonic_number=harmonic_number,npoints_x=50,npoints_y=50)
    simulation_test.print_parameters()
    simulation_test.radiation.plot(title=("radiation in a screen including central cone, harmonic %d,at D="+repr(distance))%harmonic_number)

    #
    # up to a given ring
    #
    harmonic_number = 1
    ring_number = 1
    print('create radiation a given screen including ring %d for harmonic %d'%(ring_number,harmonic_number))
    simulation_test.calculate_until_ring_number(harmonic_number=harmonic_number, ring_number=ring_number, XY_are_list=False, npoints=51)
    simulation_test.print_parameters()
    simulation_test.radiation.plot(title=" radiation in a screen containing harmonic %d up to ring %d ring"%
                                         (harmonic_number,ring_number))





def Example_spectrum_on_axis():
    from pySRU.ElectronBeam import ElectronBeam

    print("======================================================================")
    print("======      Undulator U18 from ESRF with K=1.68                =======")
    print("======================================================================")

    beam_ESRF = ElectronBeam(Electron_energy=6.0, I_current=0.2)
    ESRF18 = Undulator(K=1.68, period_length=0.018, length=2.0)

    # note that distance=None is important to get results in angles and therefore flux in ph/s/0.1%bw/mrad2
    simulation_test = create_simulation(magnetic_structure=ESRF18,electron_beam=beam_ESRF,
                        magnetic_field=None, photon_energy=None,
                        traj_method=TRAJECTORY_METHOD_ANALYTIC,Nb_pts_trajectory=None,
                        rad_method=RADIATION_METHOD_APPROX_FARFIELD, Nb_pts_radiation=101,
                        initial_condition=None, distance=None,XY_are_list=False,X=np.array([0]),Y=np.array([0]))


    x,y = simulation_test.calculate_spectrum_on_axis(use_eV=1,do_plot=1)

    # dump file
    filename = "spectrum.spec"
    f = open(filename,'w')
    f.write("#F %s\n\n#S 1 undulator flux using pySRU\n"%filename)
    f.write("#N 2\n#L Photon energy [eV]  Flux on axis [phot/s/0.1%bw/mrad2]\n")
    for i in range(x.size):
        f.write("%f  %g\n"%(x[i],y[i]))
    f.close()
    print("File written to disk: %s"%filename)

def Example_spectrum_on_slit():
    from pySRU.ElectronBeam import ElectronBeam

    print("======================================================================")
    print("======      Undulator U18 from ESRF with K=1.68                =======")
    print("======================================================================")

    beam_ESRF = ElectronBeam(Electron_energy=6.0, I_current=0.2)
    ESRF18 = Undulator(K=1.68, period_length=0.018, length=2.0)

    print('create radiation a given screen (.025 x .025 mrad) accepting central cone')
    is_quadrant = 1
    distance = None # using angles

    if is_quadrant:
        X = np.linspace(0,0.0125e-3,26)
        Y = np.linspace(0,0.0125e-3,26)
    else:
        X = np.linspace(-0.0125e-3,0.0125e-3,51)
        Y = np.linspace(-0.0125e-3,0.0125e-3,51)

    if distance != None:
        X *= distance
        Y *= distance

    simulation_test = create_simulation(magnetic_structure=ESRF18,electron_beam=beam_ESRF,
                        magnetic_field=None, photon_energy=None,
                        traj_method=TRAJECTORY_METHOD_ANALYTIC,Nb_pts_trajectory=None,
                        rad_method=RADIATION_METHOD_APPROX_FARFIELD, Nb_pts_radiation=101,
                        initial_condition=None, distance=distance,XY_are_list=False,X=X,Y=Y)

    simulation_test.print_parameters()
    simulation_test.radiation.plot(title=("radiation in a screen for first harmonic"))
    print("Integrated flux at resonance: %g photons/s/0.1bw"%(simulation_test.radiation.integration(is_quadrant=is_quadrant)))

    x,y = simulation_test.calculate_spectrum_on_slit(abscissas_array=None,use_eV=1,is_quadrant=is_quadrant)

    # dump file
    filename = "spectrum.spec"
    f = open(filename,'w')
    f.write("#F %s\n\n#S 1 undulator flux using pySRU\n"%filename)
    f.write("#N 2\n#L Photon energy [eV]  Flux on slit [phot/s/0.1%bw]\n")
    for i in range(x.size):
        f.write("%f  %g\n"%(x[i],y[i]))
    f.close()
    print("File written to disk: %s"%filename)


def Example_spectrum_on_central_cone():
    from pySRU.ElectronBeam import ElectronBeam

    print("======================================================================")
    print("======      Undulator U18 from ESRF with K=1.68                =======")
    print("======================================================================")

    beam_ESRF = ElectronBeam(Electron_energy=6.0, I_current=0.2)
    ESRF18 = Undulator(K=1.68, period_length=0.018, length=2.0)


    simulation_test = create_simulation(magnetic_structure=ESRF18,electron_beam=beam_ESRF,
                        magnetic_field=None, photon_energy=None,
                        traj_method=TRAJECTORY_METHOD_ANALYTIC,Nb_pts_trajectory=None,
                        rad_method=RADIATION_METHOD_APPROX_FARFIELD, Nb_pts_radiation=101,
                        initial_condition=None, distance=None,XY_are_list=False,X=np.array([0]),Y=np.array([0]))

    harmonic_number = 1
    x,y = simulation_test.calculate_spectrum_central_cone(harmonic_number=harmonic_number,theoretical_value=True,abscissas_array=None,use_eV=1)

    # dump file
    filename = "spectrum.spec"
    f = open(filename,'w')
    f.write("#F %s\n\n#S 1 undulator flux using pySRU\n"%filename)
    f.write("#N 2\n#L Photon energy [eV]  Flux on central cone for harmonic %d [phot/s/0.1bw]\n"%harmonic_number)
    for i,xi in enumerate(x):
        f.write("%f  %g\n"%(xi,y[i]))
    f.close()
    print("File written to disk: %s"%filename)

# TODO ne marche plus  qd on rentre a la main l'energy du photon?

def Example_list():
    from pySRU.ElectronBeam import ElectronBeam

    beam_ESRF = ElectronBeam(Electron_energy=6.0, I_current=0.2)
    ESRF18 = Undulator(K=1.68, period_length=0.018, length=2.0)

    X = np.linspace(0.0, 0.0002,1001)
    Y = np.linspace(0.0, 0.0002, 1001)
    simulation_test = create_simulation(magnetic_structure=ESRF18,electron_beam=beam_ESRF,
                                        X=X,Y=Y,XY_are_list=True)

    simulation_test.print_parameters()

    simulation_test.trajectory.plot_3D()

    simulation_test.radiation.plot()

    simulation_test.calculate_on_ring_number(ring_number_max=2)

    simulation_test.radiation.plot()

    simulation_test.radiation.plot_ring(Nb_pts=simulation_test.radiation_fact.Nb_pts)

    observation_angle = np.linspace(0.0, simulation_test.source.angle_ring_number(1, 2), 51)

    simulation_test.calculate_for_observation_angles(XY_are_list=True,observation_angle=observation_angle)

    simulation_test.radiation.plot()



if __name__ == "__main__" :
    pass
    # Example_minimum()
    # Example_meshgrid()
    # Example_meshgrid_on_central_cone_and_rings()
    # Example_spectrum_on_axis()
    # Example_spectrum_on_slit()
    Example_spectrum_on_central_cone()
    # Example_list()

