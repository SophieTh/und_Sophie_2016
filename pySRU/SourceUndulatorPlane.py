import numpy as np
import scipy.constants as codata
from scipy.special import jn,yn,jv,yv
from pySRU.MagneticField import MagneticField
from pySRU.ElectronBeam import ElectronBeam
from pySRU.MagneticStructureUndulatorPlane import MagneticStructureUndulatorPlane
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

    ###############
    # choix automatic lors de la construction d'une simulation
    ################

    #TODO useful ?
    def choose_distance_automatic(self, alpha=2):
        lim = self.magnetic_structure.length
        return lim * 10 ** alpha

    def choose_nb_pts_trajectory(self, alpha=2):
        racine = (np.pi * self.magnetic_structure.period_number() * 10 ** (alpha)) / (45.)
        n_mini = (2.0 * np.pi * self.period_number()) * np.exp(0.25 * np.log(racine))
        return n_mini*2.0

    def choose_initial_contidion_automatic(self):
        Zo=-self.magnetic_structure.length/2.-5.*self.magnetic_structure.period_length
        ic=np.array([0.0,0.0,self.electron_speed()*codata.c,0.0,0.0,Zo])
        return ic

    def choose_photon_frequency(self,harmonic_number=1):
        first_harm=self.harmonic_frequency(harmonic_number=harmonic_number)
        return first_harm

    #se demander quel angles choisir automatiquement
    def choose_angle_deflection_max(self):
        gamma=self.Lorentz_factor()
        return self.magnetic_structure.K/gamma


    ###########
    #property of the trajectory
    ##################

    def average_z_speed_in_undulator(self):
        Beta_et = 1.0 -( 1.+(self.magnetic_structure.K ** 2) / 2.0) / (2.0 * self.Lorentz_factor() ** 2)
        return Beta_et

    def analytical_times_vector(self, Nb_pts):
        to = self.magnetic_structure.length / (2. * self.average_z_speed_in_undulator() * codata.c)
        time = np.linspace(-to, to, Nb_pts)
        return time

    def construct_times_vector(self, initial_contition, Nb_pts):
        to = initial_contition[5] / (self.average_z_speed_in_undulator() * codata.c)
        if to < 0.0:
            time = np.linspace(to, -to, Nb_pts)
        elif to == 0.0:
            t1 = self.magnetic_structure.length / 2.0 + 5. * self.magnetic_structure.period_length
            time = np.linspace(to, t1)
        else:
            delta_t = to - self.magnetic_structure.length / (2. * self.average_z_speed_in_undulator() * codata.c)
            if delta_t > 0.0:
                time = np.linspace(to, to + 2. * delta_t, Nb_pts)
            else:
                time = np.linspace(to, 2. * to, Nb_pts)
        return time


    ###########
    #property of the radiation
    ##################

    def angle_deflection_max(self):
        gamma=self.Lorentz_factor()
        return self.magnetic_structure.K/gamma

    # def angle_deflection_central_cone(self):
    #     # ring1=self.angle_ring_number(harmonic_number=1,ring_number=1)
    #     # return ring1/4.
    #     N=self.magnetic_structure.period_number()
    #     gamma=self.Lorentz_factor()
    #     return 1.5/(gamma*np.sqrt(N))

    def angle_deflection_central_cone(self,harmonic_number=1):
        # this is the angle of central cone sigma
        N=self.magnetic_structure.period_number()
        gamma=self.Lorentz_factor()
        return (1 / gamma) * np.sqrt( (1 + 0.5 * self.magnetic_structure.K**2) / N / harmonic_number)


    def angle_ring_number(self,harmonic_number,ring_number):
        racine=(ring_number/harmonic_number)*(1.+(self.magnetic_structure.K**2)/2.)
        theta = np.sqrt(racine)/self.Lorentz_factor()
        return theta

    def harmonic_frequency(self,harmonic_number=1):
        gamma = self.Lorentz_factor()
        first_harm = ((2.0 * gamma ** 2) / (1.0 + (self.magnetic_structure.K ** 2) / 2.0)) * (
                (2.0 * np.pi * codata.c) / self.magnetic_structure.period_length)
        return harmonic_number*first_harm



    #######3
    #theoretical result
    ##############

    def Fn(self,n):
        if n%2==1 :
            K=self.magnetic_structure.K
            cst1=((n*K)/(1.+(K**2)/2.))**2
            cst2 = (n * K**2) / (4. + (2.*K ** 2))
            Fn=cst1*(jn(0.5*(n-1),cst2)-jn(0.5*(n+1),cst2))**2
        else :
            Fn=0.0
        return Fn

    def Qn(self,n):
        if n==0 :
            raise Exception(' the harmonic number can not be 0')
        K=self.magnetic_structure.K
        res=(1.+0.5*K**2)*self.Fn(n)/n
        return res

    #electron energy en Gev
    def critical_frequency(self,electron_energy):
        critical_energy=self.magnetic_field_strength()*(electron_energy/2.9)**2*5.59e3
        return critical_energy/codata.hbar
    #TODO TALMAN

    #TODO a changer , faire des exemples
    def describe_ring(self,distance,harmonic_number,ring_number,t=None):
        if ring_number==0 :
            X=np.array([0.0])
            Y=X
        else :
            if t is None :
                t=np.linspace(0.0,0.5*np.pi,101)
            n=harmonic_number
            theta=self.angle_ring_number(harmonic_number=n, ring_number=ring_number)
            if distance==None :
                R=theta
            else :
                R=distance*np.tan(theta)
            X=R*np.cos(t)
            Y=R*np.sin(t)
        return X,Y

    # in photon /sec /1% /mrad*mrad
    def theoretical_flux_on_axis(self,n):
        if n%2==1 :
            # see X-ray data booklet pag 2.7
            cst=1.744e14*((self.period_number()*self.Electron_energy())**2)*self.I_current()
            result=cst*self.Fn(n)
        else :
            result=0.0
        return  result

    def theoretical_flux_on_axis_omega_scan(self,n,omega_array):
        if isinstance(n,(int,float)):
            n = [n]

        result = np.zeros_like(omega_array)
        for i in range(len(n)):
            cte = self.theoretical_flux_on_axis(n[i])
            arg = np.pi * self.period_number() * (omega_array / self.harmonic_frequency(1) - n[i])
            sinc = np.sin(arg) / arg
            result += cte*sinc**2
        return  result

    def theoretical_flux_integrated_central_cone(self,n):
        # see X-ray data booklet pag 2.8
        N=self.magnetic_structure.period_number()
        I=self.I_current()
        res=1.431e14*N*I*self.Qn(n)
        return res





def Example1_undulator(undulator):

    undulator.print_parameters()

    print( " velocity average in Z direction (c units)")
    print(undulator.average_z_speed_in_undulator())
    print(" first harmonic frequency (Hz)")
    print(undulator.harmonic_frequency(1))
    print(" angle (rad) of the first ring for the first harmonic ")
    print(undulator.angle_ring_number(1,1))

    print(' magnetic field intensity (T)')
    print(undulator.magnetic_field_strength())

    print(' Magnetic field intensity in (0,0,0)')
    Bx = undulator.magnetic_field.Bx(z=0.0, y=0.0, x=0.0)
    By = undulator.magnetic_field.By(z=0.0, y=0.0, x=0.0)
    Bz = undulator.magnetic_field.Bz(z=0.0, y=0.0, x=0.0)
    B = np.array([Bx, By, Bz])
    print(B)

    print('On axis peak intensity Fn(K)')
    print(undulator.Fn(1))




if __name__ == "__main__" :
    electron_beam_test = ElectronBeam(Electron_energy=1.3, I_current=1.0)
    undulator_test = MagneticStructureUndulatorPlane( K=1.87, period_length=0.035, length=0.035 * 14)
    source=SourceUndulatorPlane(electron_beam = electron_beam_test, undulator=undulator_test)

    Example1_undulator(source)