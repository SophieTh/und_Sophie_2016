import numpy as np
import scipy.constants as codata
from scipy.special import jn,yn,jv,yv
from pySRU.MagneticField import MagneticField
from pySRU.ElectronBeam import ElectronBeam
from pySRU.MagneticStructureUndulatorPlane import MagneticStructureUndulatorPlane as Undulator
from pySRU.Source import Source,PLANE_UNDULATOR,BENDING_MAGNET



class SourceUndulatorPlane(Source):
    def __init__(self, electron_beam, undulator, magnetic_field=None):
        super(self.__class__, self).__init__(electron_beam=electron_beam,
                                             magnetic_field=magnetic_field,
                                             magnetic_structure=undulator)



    def copy(self):
        return SourceUndulatorPlane(electron_beam=self.electron_beam.copy(),
                                    magnetic_field=self.magnetic_field.copy(),
                                    undulator=self.magnetic_structure.copy())

    def period_number(self):
        return self.magnetic_structure.period_number()

    def magnetic_field_strength(self):
        return self.magnetic_structure.magnetic_field_strength()

    def angle_wave_number(self,harmonic_number,wave_number):
        racine=(wave_number/harmonic_number)*(1.+(self.magnetic_structure.K**2)/2.)
        theta = np.sqrt(racine)/self.Lorentz_factor()
        return theta

    # constante de construction

    def choose_distance_automatic(self, alpha):
        lim = self.magnetic_structure.length
        return lim * 10 ** alpha

    def choose_nb_pts_trajectory(self, alpha):
        racine = (np.pi * self.magnetic_structure.period_number() * 10 ** (alpha)) / (45.)
        n_mini = (2.0 * np.pi * self.period_number()) * np.exp(0.25 * np.log(racine))
        return n_mini*2.0

    def choose_initial_contidion_automatic(self):
        Zo=-self.magnetic_structure.length/2.-5.*self.magnetic_structure.period_length
        ic=np.array([0.0,0.0,self.electron_speed()*codata.c,0.0,0.0,Zo])
        return ic

    def harmonic_frequency(self,harmonic_number=1):
        gamma = self.Lorentz_factor()
        first_harm = ((2.0 * gamma ** 2) / (1.0 + (self.magnetic_structure.K ** 2) / 2.0)) * (
                (2.0 * np.pi * codata.c) / self.magnetic_structure.period_length)
        return harmonic_number*first_harm

    def choose_photon_frequency(self):
        first_harm=self.harmonic_frequency(harmonic_number=1)
        return first_harm

    #se demander quel angles choisir automatiquement
    def choose_angle_deflection_max(self):
        gamma=self.Lorentz_factor()
        return 1.5*self.magnetic_structure.K/gamma

    def angle_deflection_max(self):
        gamma=self.Lorentz_factor()
        return self.magnetic_structure.K/gamma

    def angle_deflection_central_cone(self):
        wave1=self.angle_wave_number(harmonic_number=1,wave_number=1)
        return wave1/4.

    def average_z_speed_in_undulator(self):
        Beta_et = 1.0 - (1.0 / (2.0 * self.Lorentz_factor() ** 2)) * (
            1.0 + (self.magnetic_structure.K** 2) / 2.0)
        return Beta_et

    def analytical_times_vector(self,Nb_pts):
        to = self.magnetic_structure.length / (2.*self.average_z_speed_in_undulator()*codata.c)
        time=np.linspace(-to,to,Nb_pts)
        return time

    def construct_times_vector(self, initial_contition,Nb_pts):
        to = initial_contition[5]/ (self.average_z_speed_in_undulator() * codata.c)
        if to < 0.0 :
            time = np.linspace(to, -to, Nb_pts)
        elif to==0.0 :
            t1=self.magnetic_structure.length/2.0 + 5.*self.magnetic_structure.period_length
            time=np.linspace(to,t1)
        else :
            delta_t=to-self.magnetic_structure.length/(2.*self.average_z_speed_in_undulator() * codata.c)
            if delta_t >0.0 :
                time= np.linspace(to,to+2.*delta_t,Nb_pts)
            else :
                time=np.linspace(to,2.*to,Nb_pts)
        return time

    def Fn(self,n):
        K=self.magnetic_structure.K
        cst1=((n*K)/(1.+(K**2)/2.))**2
        cst2 = (n * K**2) / (4. + (2.*K ** 2))
        Fn=cst1*(jn(0.5*(n-1),cst2)-jn(0.5*(n+1),cst2))**2
        return Fn

    #electron energy en Gev
    def critical_frequency(self,electron_energy):
        critical_energy=self.magnetic_field_strength()*(electron_energy/2.9)**2*5.59e3
        return critical_energy/codata.hbar
    #TODO TALMAN

    #TODO a changer , faire des exemples
    def describe_wave(self,distance,harmonic_number,wave_number,t=None):
        if wave_number==0 :
            X=np.array([0.0])
            Y=X
        else :
            if t==None :
                t=np.linspace(0.0,0.5*np.pi,101)
            n=harmonic_number
            theta=self.angle_wave_number(harmonic_number=n, wave_number=wave_number)
            R=distance*np.tan(theta)
            X=R*np.cos(t)
            Y=R*np.sin(t)
        return X,Y

    # in photon /sec /1% /mrad*mrad
    def theorical_flux_on_axis(self,n):
        if n%2==1 :
            cst=1.744e14*((self.period_number()*self.Electron_energy())**2)*self.I_current()
            result=cst*self.Fn(n)
        else :
            result=0.0
        return  result

def Exemple1_undulator(undulator):

    undulator.print_parameters()

    print( " speed average in Z direction (m/s)/c")
    print(undulator.average_z_speed_in_undulator())
    print(" first frequency")
    print(undulator.harmonic_frequency(1))
    print(" angle of the first wave for the first number (radian) ")
    print(undulator.angle_wave_number(1,1))

    print(' magnetic field intensity (T)')
    print(undulator.magnetic_field_strength())

    print(' Magnetic field intensity in (0,0,0)')
    Bx = undulator.magnetic_field.Bx(z=0.0, y=0.0, x=0.0)
    By = undulator.magnetic_field.By(z=0.0, y=0.0, x=0.0)
    Bz = undulator.magnetic_field.Bz(z=0.0, y=0.0, x=0.0)
    B = np.array([Bx, By, Bz])
    print(B)

    print('Theorical flux on axis')
    print(undulator.Fn(1))




if __name__ == "__main__" :
    electron_beam_test = ElectronBeam(Electron_energy=1.3, I_current=1.0)
    undulator_test=Undulator( K=1.87, period_length=0.035, length=0.035 * 14)
    source=SourceUndulatorPlane(electron_beam=electron_beam_test, undulator=undulator_test)

    Exemple1_undulator(source)