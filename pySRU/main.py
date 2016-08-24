import numpy as np
import scipy.constants as codata
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import scipy.integrate as integrate
import scipy.interpolate as interpolate
import time
from pySRU.MagneticStructureUndulatorPlane import MagneticStructureUndulatorPlane as Undulator
from pySRU.MagneticStructureBendingMagnet import MagneticStructureBendingMagnet  as BM
from pySRU.ElectronBeam import ElectronBeam
from pySRU.SourceUndulatorPlane import SourceUndulatorPlane
from pySRU.Simulation import Simulation ,create_simulation
from pySRU.TrajectoryFactory import TRAJECTORY_METHOD_ANALYTIC,TRAJECTORY_METHOD_ODE,\
                                        TRAJECTORY_METHOD_INTEGRATION
from pySRU.RadiationFactory import RADIATION_METHOD_NEAR_FIELD, RADIATION_METHOD_APPROX,\
                                RADIATION_METHOD_FARFIELD, RADIATION_METHOD_APPROX_FARFIELD


eV_to_J=1.602176487e-19
######################################################

#BM_test=MagneticStructureBendingMagnet(Bo=0.8, div=5e-3, R=25.0)

E_1=7876.0

beam_test=ElectronBeam(Electron_energy=1.3, I_current=1.0)
beam_ESRF=ElectronBeam(Electron_energy=6.0, I_current=0.2)
und_test=Undulator(  K = 1.87,  period_length= 0.035, length=0.035 * 14)
ESRF18=Undulator( K = 1.68, period_length = 0.018, length=2.0)
ESRFBM=BM(Bo=0.8,L=0.1249994791673177)

#omega=E_1*eV_to_J/codata.hbar

# # recuperation des donnees de B en array en fonction de z
# reference=np.load("../tests/x_ray_booklet_field.npz")
# Z=reference['ct']
# Z -= (Z[-1])/2.0
# By=reference['B_y']


def compare_2_traj_method(magnetic_structure, electron_beam, rad_method, traj_method_ref, traj_method_test, distance=None):
#

    #TODO trouver un moyen de le faire directement sur le central cone...
    sim_test = create_simulation(magnetic_structure=magnetic_structure, electron_beam=electron_beam,formula=1,
                             rad_method=rad_method, traj_method=traj_method_ref,distance=distance)
    sim_test.calculate_on_central_cone()


    ref = sim_test.copy()

    rad_max = (sim_test.radiation.max())
    print('rad max reference')
    print(rad_max)

    sim_test.change_trajectory_method(traj_method_test)

    traj_err=ref.trajectory.relativ_difference_with(sim_test.trajectory)

    traj_err.plot()

    sim_test.radiation.plot()
    sim_test.radiation.plot()

    rad_err=ref.radiation.relativ_difference_with(sim_test.radiation)
    print('rad_max simulation test')
    print(rad_err.max())
    rad_err.plot()

    omega_array,spectre_ref=sim_test.spectre()
    spectre_test=sim_test.spectre(omega_array=omega_array)[1]
    plt.plot(omega_array,spectre_ref)
    plt.plot(omega_array,spectre_test)
    plt.show()


def compare_interpolation(magnetic_structure, electron_beam,traj_method=TRAJECTORY_METHOD_ODE, distance=None):

    sim = create_simulation(magnetic_structure=magnetic_structure, electron_beam=electron_beam,
                            traj_method=traj_method,distance=100)

    z_ref = sim.trajectory.z * codata.c

    z0 = sim.trajectory.z[0] * codata.c
    z1 = sim.trajectory.z[-1] * codata.c
    z_test = np.linspace(z0, z1, 100)
    By_array = sim.source.magnetic_field.By(x=0.0, y=0.0, z=z_test)
    plt.plot(z_test, By_array)
    plt.show()

    By_linear = interpolate.interp1d(z_test, By_array, kind='linear')
    By_spline = interpolate.interp1d(z_test, By_array, kind='cubic')

    plt.plot(z_ref, sim.source.magnetic_field.By(x=0.0, y=0.0, z=z_ref))
    plt.plot(z_ref, By_linear(z_ref))
    plt.plot(z_ref, By_spline(z_ref))
    plt.show()


#

#compare_2_traj_method(magnetic_structure=ESRF18,electron_beam=beam_ESRF,rad_method=RADIATION_METHOD_APPROX_FARFIELD,
                    #traj_method_ref=TRAJECTORY_METHOD_ANALYTIC,traj_method_test=TRAJECTORY_METHOD_ODE,distance=100)

# sim=create_simulation(magnetic_structure=und_test,electron_beam=beam_test,Nb_pts_trajectory=10000,traj_method=TRAJECTORY_METHOD_ODE)
#
# beta_et=sim.source.average_z_speed_in_undulator()
# print('beta et')
# print (beta_et)
# beta2=beta_et + sim.source.magnetic_structure.K**2/(4.*sim.source.Lorentz_factor()**2)
# print('beta2')
# print (beta2)
# beta3=beta_et - sim.source.magnetic_structure.K**2/(4.*sim.source.Lorentz_factor()**2)
# print('beta3')
# print (beta3)
# #sim.change_trajectory_method(TRAJECTORY_METHOD_ODE)
# t=sim.trajectory.t
# z=sim.trajectory.z
# vz=sim.trajectory.v_z
# zo=und_test.length/2
# Mv=0
# Mz=0
# compt=0
# for i in range(len(z)) :
#     if z[i] > -zo and z[i] < zo :
#         Mz += z[i]/t[i]
#         Mv += vz[i]
#         compt += 1
#
# Mv /= compt
# Mz /= compt
# print('Mv')
# print(Mv)
# print('Mz')
# print(Mz)
#
# rel=(und_test.period_length*und_test.K**2)/(8*np.pi*sim.source.Lorentz_factor()**2*beta_et)
# print(rel)
# rel2=(und_test.K**2)/(sim.source.Lorentz_factor()*beta_et)
# print(rel2)
# #sim.trajectory.plot()





def observation_err_formule12_analytic(und,beam):
    sim_test = create_simulation(magnetic_structure=und, electron_beam=beam,
                                 traj_method=TRAJECTORY_METHOD_ANALYTIC,formula=1,
                                 Nb_pts_trajectory=1000,distance=100)

    sim_test.change_harmonic_number(5)


    print(sim_test.source.theorical_flux_on_axis(photon_frequency=sim_test.radiation_fact.omega))
    sim_test.radiation.plot()
    sim_test2=create_simulation(magnetic_structure=und, electron_beam=beam,
                                 traj_method=TRAJECTORY_METHOD_ANALYTIC,formula=2,
                                 Nb_pts_trajectory=1000,distance=100)
    sim_test2.change_omega(sim_test.radiation_fact.omega)
    print(sim_test2.radiation.max())
    sim_test2.radiation.plot()
    diff=sim_test.radiation.relativ_difference_with(sim_test2.radiation)
    print(diff.max())
    diff.plot()


observation_err_formule12_analytic(und_test,beam_test)

def erreur_entre_ODE_ANALYTIC_sans_mm_condition_initial(und,beam,rad):
    sim_test_analy = create_simulation(magnetic_structure=und, electron_beam=beam,
                                 traj_method=TRAJECTORY_METHOD_ANALYTIC, formula=1,rad_method=rad,
                                 Nb_pts_trajectory=20000,distance=100, X=0.002,Y=0.002)

    #sim_test_analy.radiation.plot(title='ESRF18 radiation intensity for analitycal trajectory')
    print('inital donc anali')
    print(sim_test_analy.trajectory_fact.initial_condition)
    sim_test_ODE=create_simulation(magnetic_structure=und, electron_beam=beam,
                                 traj_method=TRAJECTORY_METHOD_ODE, formula=1,rad_method=rad,
                                 Nb_pts_trajectory=20000,distance=100,X=0.002,Y=0.002)

    #sim_test_analy.radiation.plot(title='ESRF18 radiation intensity for ODE trajectory')


    t_milieu=sim_test_analy.trajectory.t
    delta_t=t_milieu[1]-t_milieu[0]
    t_debut=(-und.length*0.5-und.period_length*5.)/(sim_test_ODE.source.average_z_speed_in_undulator()*codata.c)
    t_fin = (und.length * 0.5 + und.period_length * 5.) / (
    sim_test_ODE.source.average_z_speed_in_undulator() * codata.c)
    t_avt = np.arange(t_debut,t_milieu[0],delta_t)
    indice_t_milieu=len(t_avt)
    t_apr = np.arange(t_milieu[-1]+delta_t, t_fin+delta_t, delta_t)
    new_t=np.concatenate((t_avt,t_milieu,t_apr))
    sim_test_ODE.trajectory = sim_test_ODE.trajectory_fact.create_from_source(source=sim_test_ODE.source,
                                                                              t=new_t)
    # Z_srw=np.linspace(-0.45,0.45,15000)
    # plt.plot(Z_srw, sim_test_ODE.source.magnetic_field.By(z=Z_srw, x=0.0, y=0.0),label='SRW magnetic field')
    # Z_anali = np.linspace(-und.length/2., und.length/2., 10000)
    # plt.plot(Z_anali,sim_test_ODE.source.magnetic_field.By(z=Z_anali,x=0.0,y=0.0),label='analitycal magnetic field')
    # plt.legend()
    # plt.xlabel('z')
    # plt.ylabel('By')
    # plt.title('y component of the magnetiq field')
    # plt.show()

    # sim_test_ODE.trajectory.multiply_by(codata.c)
    # sim_test_analy.trajectory.multiply_by(codata.c)

    print(all(sim_test_ODE.trajectory.t==new_t))
    # plt.plot(new_t,sim_test_ODE.trajectory.x,label='trajectory from ODE')
    # plt.plot(sim_test_analy.trajectory.t, sim_test_analy.trajectory.x,label='analytical trajectory')
    # plt.legend()
    # plt.xlabel('t')
    # plt.ylabel('X')
    # plt.title('X component of the trajectory')
    # plt.show()


    # plt.plot(new_t, sim_test_ODE.trajectory.z, label='trajectory from ODE')
    # plt.plot(sim_test_analy.trajectory.t, sim_test_analy.trajectory.z, label='analytical trajectory')
    # plt.legend()
    # plt.xlabel('t')
    # plt.ylabel('Z')
    # plt.title('Z component of the trajectory')
    # plt.show()


    # beta_et=sim_test_analy.source.average_z_speed_in_undulator()*codata.c
    # plt.plot(new_t, sim_test_ODE.trajectory.v_z, label='trajectory from ODE')
    # plt.plot(sim_test_analy.trajectory.t,
    #          sim_test_analy.trajectory.v_z, label='analytical trajectory')
    # plt.legend()
    # plt.xlabel('t')
    # plt.ylabel('Z-B*t')
    # plt.title('z component of the velocity')
    # plt.show()

    #sim_test_analy.radiation.plot()
    # print ('ODE')
    # sim_test_ODE.radiation.plot()
    #
    diff=sim_test_analy.radiation.relativ_difference_with(sim_test_ODE.radiation)
    diff.plot(title="relative difference between radiation from ODE and analitycal trajectory")

    # ## recherche des condition initial de l'ondulateur :
    # Zo= -und.length/2.
    # z=sim_test_ODE.trajectory.z*codata.c
    # i=0
    # while i < (len(z)) and z[i]<= Zo :
    #     i +=1
    #
    i=indice_t_milieu
    print('t milieu')
    print(new_t[i])
    print(t_milieu[0])
    initial_cond=np.array([sim_test_ODE.trajectory.v_x[i],
                               sim_test_ODE.trajectory.v_y[i],
                               sim_test_ODE.trajectory.v_z[i],
                               sim_test_ODE.trajectory.x[i],
                               sim_test_ODE.trajectory.y[i],
                               sim_test_ODE.trajectory.z[i]])*codata.c
    print('inital donc construite')
    print(initial_cond)



    print(sim_test_analy.trajectory_fact.initial_condition-initial_cond)

    sim_test_ODE.trajectory_fact.initial_condition = initial_cond
    sim_test_ODE.trajectory=sim_test_ODE.trajectory_fact.create_from_source(source=sim_test_ODE.source,
                                                                            t=sim_test_analy.trajectory.t)
    rad_old=sim_test_ODE.radiation.copy()
    sim_test_ODE.radiation.intensity=sim_test_ODE.radiation_fact.calculate_radiation_intensity(source=sim_test_ODE.source,
                                trajectory=sim_test_ODE.trajectory,X_array=rad_old.X,Y_array=rad_old.Y,
                                                                                     distance=rad_old.distance)


    print('inital cond ODE')
    print(sim_test_ODE.trajectory_fact.initial_condition)
    b_a=sim_test_analy.trajectory.average_vz()
    print(' vz moyen analitiq')
    print(sim_test_analy.trajectory.average_vz())
    print(' vz moyen ODE')
    b_o=sim_test_ODE.trajectory.average_vz()
    print(sim_test_ODE.trajectory.average_vz()/sim_test_ODE.source.electron_speed())
    print(' err rel des deux oyenne par rapport a la trajectoire analytiq')
    print((b_a-b_o)/b_a)
    print(' err rel des deux oyenne par rapport beta *')
    print((b_a - b_o)/sim_test_ODE.source.electron_speed())

    print('err de t :')
    print(abs(sim_test_analy.trajectory.t-sim_test_ODE.trajectory.t).max()/ np.abs(sim_test_ODE.trajectory.t).max())


    # print ('beta')
    # print(sim_test_ODE.source.electron_speed())
    # print ('beta et')
    # print(sim_test_ODE.source.average_z_speed_in_undulator())
    # t=sim_test_analy.trajectory.t
    # sim_test_analy.trajectory.plot_3D()
    # sim_test_ODE.trajectory.plot_3D()
    # sim_test_analy.trajectory.plot()
    # sim_test_ODE.trajectory.plot()
    # plt.plot(t,sim_test_ODE.trajectory.z-b_o*t)
    # plt.plot(t, sim_test_analy.trajectory.z - b_a * t)
    # plt.show()
    # plt.plot(t,sim_test_ODE.trajectory.v_z-b_o)
    # plt.plot(t, sim_test_analy.trajectory.v_z - b_a )
    # plt.show()
    # sim_test_analy.trajectory.plot_2_trajectory(sim_test_ODE.trajectory)
    # diff=sim_test_ODE.trajectory.difference_with(sim_test_analy.trajectory)
    # diff1=sim_test_ODE.trajectory.relativ_difference_with(sim_test_analy.trajectory)
    # diff1.plot()
    # second_terme=sim_test_analy.trajectory

    # second_terme.z = second_terme.z-beta_et*t
    # second_terme.v_z = second_terme.v_z - beta_et
    # diff.plot_2_trajectory(second_terme)


    # print('rad ODE max')
    # print(sim_test_ODE.radiation.max())
    # sim_test_ODE.radiation.plot()
    diff2=sim_test_ODE.radiation.relativ_difference_with(sim_test_analy.radiation)
    diff2.plot()

    # spectre_analy,omega_array=sim_test_analy.spectre()
    # spectre_ODE,omega_array=sim_test_ODE.spectre(omega_array=omega_array)
    # plt.plot(omega_array,spectre_analy)
    # plt.plot(omega_array, spectre_ODE)
    # plt.show()
    #
    # diff_spectre=np.abs(spectre_analy-spectre_ODE)
    # plt.plot(omega_array,diff_spectre)
    # plt.show()

#cas B
#erreur_entre_ODE_ANALYTIC_sans_mm_condition_initial(und=ESRF18,beam=beam_ESRF,rad=RADIATION_METHOD_NEAR_FIELD)

def erreur_entre_ODE_ANALYTIC_sans_mm_condition_initial_suivant_nb_pts(und,beam,nb_pts):


    sim_test_ODE = create_simulation(magnetic_structure=und, electron_beam=beam,
                                     traj_method=TRAJECTORY_METHOD_ODE, formula=2,
                                     Nb_pts_trajectory=30000)
    ## recherche des condition initial de l'ondulateur :
    Zo = -und.length / 2.
    z = sim_test_ODE.trajectory.z * codata.c
    i = 0
    print('ok')
    while i < (len(z)) and z[i] <= Zo:
        i += 1

    initial_cond = np.array([sim_test_ODE.trajectory.v_x[i],
                             sim_test_ODE.trajectory.v_y[i],
                             sim_test_ODE.trajectory.v_z[i],
                             sim_test_ODE.trajectory.x[i],
                             sim_test_ODE.trajectory.y[i],
                             sim_test_ODE.trajectory.z[i]]) * codata.c

    print('initial condition')
    print(initial_cond)
    sim_test_ODE.trajectory_fact.initial_condition = initial_cond

    beta=sim_test_ODE.source.electron_speed()
    average_error=np.zeros(len(nb_pts))
    rad_max_ODE = np.zeros(len(nb_pts))
    rad_max_analy = np.zeros(len(nb_pts))

    error_traj=np.zeros((10,len(nb_pts)))
    for i in range(len(nb_pts)) :
        sim_test_analy = create_simulation(magnetic_structure=und, electron_beam=beam,
                                       traj_method=TRAJECTORY_METHOD_ANALYTIC, formula=2,
                                       Nb_pts_trajectory=nb_pts[0])
        sim_test_ODE.trajectory=sim_test_ODE.trajectory_fact.create_from_source(source=sim_test_ODE.source,
                                                                            t=sim_test_analy.trajectory.t)
        rad_old=sim_test_ODE.radiation
        sim_test_ODE.radiation.intensity=sim_test_ODE.radiation_fact.calculate_radiation_intensity(source=sim_test_ODE.source,
                                trajectory=sim_test_ODE.trajectory,X_array=rad_old.X,Y_array=rad_old.Y,
                                                                                     distance=rad_old.distance)



        error_traj[:,i]=sim_test_analy.trajectory.error_max(sim_test_ODE.trajectory)
        b_a = sim_test_analy.trajectory.average_vz()
        b_o = sim_test_ODE.trajectory.average_vz()
        rad_max_analy[i]=sim_test_analy.radiation.max()
        rad_max_ODE[i]=sim_test_ODE.radiation.max()
        average_error[i]=np.abs(b_a-b_o)

    theoric=np.zeros(len(nb_pts))+sim_test_ODE.source.flux_on_axis_theoric(omega=sim_test_ODE.radiation_fact.omega)
    plt.plot(nb_pts,theoric)
    plt.plot(nb_pts,rad_max_analy)
    plt.plot(nb_pts,rad_max_ODE)
    plt.show()

    rel_error_average=average_error/beta
    plt.plot(nb_pts,rel_error_average)
    plt.show()
    rel_error_average=average_error/beta
    plt.plot(nb_pts,rel_error_average)
    plt.show()
    error_trajectory=sim_test_analy.trajectory.copy()
    error_trajectory.convert(error_traj)
    error_trajectory.plot()

#nb_pts=np.linspace(200,30000,101)
#erreur_entre_ODE_ANALYTIC_sans_mm_condition_initial_suivant_nb_pts(und=und_test,beam=beam_test,nb_pts=nb_pts)


def erreur_analy_ODE(magn,beam):
    sim_test_analy = create_simulation(magnetic_structure=magn, electron_beam=beam,
                                       traj_method=TRAJECTORY_METHOD_ANALYTIC, formula=1,
                                       Nb_pts_trajectory=20000)
    sim_test_ODE = create_simulation(magnetic_structure=magn, electron_beam=beam,
                                       traj_method=TRAJECTORY_METHOD_ODE, formula=1,
                                       Nb_pts_trajectory=20000)

    sim_test_ODE.trajectory.plot_2_trajectory(sim_test_analy.trajectory)
    sim_test_analy.radiation.plot()
    sim_test_ODE.radiation.plot()

    # nb_pts = np.linspace(200, 30000, 21)
    # error_rad, traj_error_traj=sim_test_analy.error_trajectory_method(method=TRAJECTORY_METHOD_ODE,nb_pts=nb_pts)
    # plt.plot(nb_pts,error_rad)
    # plt.show()
    # traj_error_traj.plot()
    t_milieu=sim_test_analy.trajectory.t
    delta_t=t_milieu[1]-t_milieu[0]
    t_debut=(-magn.length*0.5-magn.period_length*5.)/(sim_test_ODE.source.average_z_speed_in_undulator()*codata.c)
    t_fin = (magn.length * 0.5 + magn.period_length * 5.) / (
    sim_test_ODE.source.average_z_speed_in_undulator() * codata.c)
    t_avt = np.arange(t_debut,t_milieu[0],delta_t)
    indice_t_milieu=len(t_avt)
    t_apr = np.arange(t_milieu[-1]+delta_t, t_fin+delta_t, delta_t)
    new_t=np.concatenate((t_avt,t_milieu,t_apr))
    sim_test_ODE.trajectory = sim_test_ODE.trajectory_fact.create_from_source(source=sim_test_ODE.source,
                                                                              t=new_t)

    sim_test_ODE.trajectory.plot_2_trajectory(sim_test_analy.trajectory)

    i=indice_t_milieu
    print('t milieu')
    print(new_t[i])
    print(t_milieu[0])
    initial_cond=np.array([sim_test_ODE.trajectory.v_x[i],
                               sim_test_ODE.trajectory.v_y[i],
                               sim_test_ODE.trajectory.v_z[i],
                               sim_test_ODE.trajectory.x[i],
                               sim_test_ODE.trajectory.y[i],
                               sim_test_ODE.trajectory.z[i]])*codata.c
    print('inital donc construite')
    print(initial_cond)


#erreur_analy_ODE(magn=ESRFBM,beam=beam_ESRF)




#cas A
#erreur_analy_ODE_BM(BM=ESRFBM,beam=beam_ESRF)

def example_non_error_ODE_vs_analytique():
    sim_test_analy = create_simulation(magnetic_structure=ESRF18, electron_beam=beam_ESRF,
                                       traj_method=TRAJECTORY_METHOD_ANALYTIC, formula=1,
                                       Nb_pts_trajectory=20000)
    sim_test_ODE=sim_test_analy.copy()
    sim_test_ODE.trajectory_fact.method=TRAJECTORY_METHOD_ODE
    sim_test_ODE.trajectory=sim_test_ODE.trajectory_fact.create_from_source(source=sim_test_ODE.source,t=sim_test_analy.trajectory.t)
    sim_test_ODE.radiation.intentsity=sim_test_ODE.radiation_fact.calculate_radiation_intensity(source=sim_test_ODE.source,
                                                                                      trajectory=sim_test_ODE.trajectory,
                                                                                      X_array=sim_test_ODE.radiation.X,
                                                                                      Y_array=sim_test_ODE.radiation.Y,
                                                                                      distance=sim_test_ODE.radiation.distance)

    sim_test_analy.trajectory.plot_2_trajectory(sim_test_ODE.trajectory)
    diff=sim_test_analy.trajectory.relativ_difference_with(sim_test_ODE.trajectory)
    diff.plot()
    diff_rad=sim_test_analy.radiation.relativ_difference_with(sim_test_ODE.radiation)
    diff_rad.plot()

#example_non_error_ODE_vs_analytique()


#cas C
def example_non_error_ODE_vs_int():
    sim_test_analy = create_simulation(magnetic_structure=und_test, electron_beam=beam_test,
                                       traj_method=TRAJECTORY_METHOD_ODE, formula=1,
                                       Nb_pts_trajectory=20000)
    sim_test_ODE=sim_test_analy.copy()
    sim_test_ODE.trajectory_fact.method=TRAJECTORY_METHOD_INTEGRATION
    sim_test_ODE.trajectory=sim_test_ODE.trajectory_fact.create_from_source(source=sim_test_ODE.source,t=sim_test_analy.trajectory.t)
    sim_test_ODE.radiation.intentsity=sim_test_ODE.radiation_fact.calculate_radiation_intensity(source=sim_test_ODE.source,
                                                                                      trajectory=sim_test_ODE.trajectory,
                                                                                      X_array=sim_test_ODE.radiation.X,
                                                                                      Y_array=sim_test_ODE.radiation.Y,
                                                                                      distance=sim_test_ODE.radiation.distance)

    sim_test_analy.trajectory.plot_2_trajectory(sim_test_ODE.trajectory)
    diff=sim_test_analy.trajectory.relativ_difference_with(sim_test_ODE.trajectory)
    diff.plot()
    diff_rad=sim_test_analy.radiation.relativ_difference_with(sim_test_ODE.radiation)
    diff_rad.plot()

#example_non_error_ODE_vs_int()


#cas D
def example_non_error_analitic_vs_int():
    sim_test_analy = create_simulation(magnetic_structure=und_test, electron_beam=beam_test,
                                       traj_method=TRAJECTORY_METHOD_ANALYTIC, formula=1,
                                       Nb_pts_trajectory=210)
    sim_test_ODE=sim_test_analy.copy()
    sim_test_ODE.trajectory_fact.method=TRAJECTORY_METHOD_INTEGRATION
    sim_test_ODE.trajectory=sim_test_ODE.trajectory_fact.create_from_source(source=sim_test_ODE.source,t=None)#t=sim_test_analy.trajectory.t)
    sim_test_ODE.radiation.intentsity=sim_test_ODE.radiation_fact.calculate_radiation_intensity(source=sim_test_ODE.source,
                                                                                      trajectory=sim_test_ODE.trajectory,
                                                                                      X_array=sim_test_ODE.radiation.X,
                                                                                      Y_array=sim_test_ODE.radiation.Y,
                                                                                      distance=sim_test_ODE.radiation.distance)

    sim_test_analy.trajectory.plot_2_trajectory(sim_test_ODE.trajectory)
    diff=sim_test_analy.trajectory.relativ_difference_with(sim_test_ODE.trajectory)
    diff.plot()
    sim_test_ODE.radiation.plot()
    sim_test_analy.radiation.plot()
    diff_rad=sim_test_analy.radiation.difference_with(sim_test_ODE.radiation)
    diff_rad.plot()

#example_non_error_analitic_vs_int()

#cas D
def erreur_analy_int_und(und,beam):
    sim_test_analy = create_simulation(magnetic_structure=und, electron_beam=beam,
                                       traj_method=TRAJECTORY_METHOD_ANALYTIC, formula=1,
                                       Nb_pts_trajectory=20000)
    nb_pts = np.linspace(2000, 30000, 21)
    error_rad, traj_error_traj=sim_test_analy.error_trajectory_method(method=TRAJECTORY_METHOD_INTEGRATION,nb_pts=nb_pts)
    plt.plot(nb_pts,error_rad)
    plt.show()
    traj_error_traj.plot()


#erreur_analy_int_und(und=und_test,beam=beam_test)

def erreur_analy_int_ODE(und,beam):
    sim_test_analy = create_simulation(magnetic_structure=und, electron_beam=beam,
                                       traj_method=TRAJECTORY_METHOD_ODE, formula=1,
                                       Nb_pts_trajectory=20000)
    nb_pts = np.linspace(2000, 30000, 21)
    error_rad, traj_error_traj=sim_test_analy.error_trajectory_method(method=TRAJECTORY_METHOD_INTEGRATION,nb_pts=nb_pts)
    plt.plot(nb_pts,error_rad)
    plt.show()
    traj_error_traj.plot()

# z=np.linspace(-0.99999997,0.99999997,10000)
# k_u=2.*np.pi/(0.035)
# f=np.sin(k_u*z)*np.sin(2.*k_u*z)-np.cos(3.*k_u*z)
# plt.plot(z,f)
# plt.show()
#erreur_analy_int_und(und=und_test,beam=beam_test)


def erreur_near_ff(und,beam):
    sim_test_analy = create_simulation(magnetic_structure=und, electron_beam=beam,
                                       traj_method=TRAJECTORY_METHOD_ANALYTIC, formula=2,NB_pts_rad=51,
                                       Nb_pts_trajectory=10000,distance=100)
    rad_max=sim_test_analy.radiation.max()
    print(rad_max)
    rad_max_theo=sim_test_analy.source.theorical_flux_on_axis(photon_frequency=sim_test_analy.radiation_fact.omega)
    print(rad_max_theo)
    d = np.linspace(50,100, 26)
    error_rad=sim_test_analy.error_radiation_method_distance(method=RADIATION_METHOD_NEAR_FIELD,D=d)
    plt.plot(d,error_rad)
    plt.xlabel("distance")
    plt.ylabel("error")
    plt.title('absolut error between two radiation method')
    plt.show()
    plt.plot(d,error_rad/rad_max)
    plt.xlabel("distance")
    plt.ylabel("error")
    plt.title('relative error between two radiation method')
    plt.show()



#erreur_near_ff(und=ESRF18,beam=beam_ESRF)

def example_erreur_near_ff(und,beam):
    sim_test_analy = create_simulation(magnetic_structure=und, electron_beam=beam,
                                       traj_method=TRAJECTORY_METHOD_ANALYTIC, formula=1,
                                       Nb_pts_trajectory=20000, distance=20)
    rad_max = sim_test_analy.radiation.max()
    print(rad_max)
    sim_test_2=sim_test_analy.copy()
    sim_test_2.change_radiation_method(RADIATION_METHOD_NEAR_FIELD)
    diff=sim_test_analy.radiation.relativ_difference_with(sim_test_2.radiation)
    sim_test_analy.radiation.plot()
    sim_test_2.radiation.plot()
    diff.plot()

#example_erreur_near_ff(ESRF18,beam_ESRF)




def erreur_near_ff_central_cone(und,beam):
    sim_test_analy = create_simulation(magnetic_structure=und, electron_beam=beam,
                                       traj_method=TRAJECTORY_METHOD_ANALYTIC, formula=1,
                                       Nb_pts_trajectory=20000,distance=100)
    #sim_test_analy.calculate_on_central_cone()
    sim_test_analy.change_XY_radiation(X=sim_test_analy.radiation.X / 5., Y=sim_test_analy.radiation.Y / 5.)
    rad_max=sim_test_analy.radiation.max()
    print(rad_max)
    rad_max_theo=sim_test_analy.source.theorical_flux_on_axis(photon_frequency=sim_test_analy.radiation_fact.omega)
    print(rad_max_theo)
    d = np.linspace(5,50, 21)
    error_rad=sim_test_analy.error_radiation_method_distance(method=RADIATION_METHOD_NEAR_FIELD,D=d)
    plt.plot(d,error_rad)
    plt.show()
    plt.plot(d,error_rad/rad_max)
    plt.show()

#erreur_near_ff_central_cone(und=ESRF18,beam=beam_ESRF)


def example_erreur_near_ff_central_cone(und,beam):

    sim_test_analy = create_simulation(magnetic_structure=und, electron_beam=beam,
                                       traj_method=TRAJECTORY_METHOD_ANALYTIC, formula=1,
                                       Nb_pts_trajectory=20000, distance=10)
    sim_test_analy.change_XY_radiation(X=sim_test_analy.radiation.X/5.,Y=sim_test_analy.radiation.Y/5.)


    rad_max = sim_test_analy.radiation.max()
    print(rad_max)
    sim_test_2=sim_test_analy.copy()
    sim_test_2.change_radiation_method(RADIATION_METHOD_NEAR_FIELD)
    diff=sim_test_analy.radiation.relativ_difference_with(sim_test_2.radiation)
    sim_test_analy.radiation.plot()
    sim_test_2.radiation.plot()
    diff.plot()


#example_erreur_near_ff_central_cone(ESRF18,beam_ESRF)




def erreur_nb_pts_trajectory(und,beam):
    X=np.array([0.0])
    Y=np.array([0.0])
    sim_test_analy = create_simulation(magnetic_structure=und, electron_beam=beam,
                                       traj_method=TRAJECTORY_METHOD_ANALYTIC, formula=1,
                                       Nb_pts_trajectory=200, distance=10,X=X,Y=Y)
    #sim_test_analy.calculate_until_wave_number()
    rad_theo=sim_test_analy.source.theorical_flux_on_axis(photon_frequency=sim_test_analy.radiation_fact.omega)
    nb_pts=np.arange(1000,5000,100)
    res =np.zeros(len(nb_pts))
    for i in range(len(nb_pts)):
        print(i)
        sim_test_analy.change_Nb_pts_trajectory(nb_pts[i])
        res[i]=sim_test_analy.radiation.max()

    #plt.plot(nb_pts,np.ones(len(nb_pts))*rad_theo)
    print(rad_theo)
    print(sim_test_analy.radiation.max())
    print(res)
    print(res-rad_theo)

    plt.show()


    plt.show()


    dif_REL=np.abs(res-rad_theo)/rad_theo
    plt.plot(nb_pts, res)
    plt.xlabel("nb point")
    plt.ylabel("error")
    plt.title('relative error between two radiation method')
    plt.show()
    plt.plot(nb_pts, dif_REL)
    plt.ylabel("error")
    plt.title('absolut error between two radiation method')
    plt.show()


#erreur_nb_pts_trajectory(und_test,beam_test)
# #sim_test.print_parameters()
#
# sim_test.change_omega(sim_test.source.critical_frequency()*0.5)
#
# print('omega')
# print(sim_test.radiation_fact.omega)
# print('omega critic 2')
# print(sim_test.source.critical_frequency2())
# print('arc length')
# print(sim_test.source.arc_length()/codata.c)
# print('physical length')
# print((sim_test.trajectory.z[-1]-sim_test.trajectory.z[0]))
# #
#
#
# #sim_test.trajectory.plot()
# print(type(sim_test.radiation))
# sim_test.trajectory.plot_3D()
# sim_test.radiation.plot()
# # #
# print('intensity[0][0]/1e13')
# print(sim_test.radiation.intensity[0,0]/1e13)
# print('radiation.max()/1e13')
# print(sim_test.radiation.max()/1e13)
# omega=sim_test.radiation_fact.omega
# print('flux_on_axis_theoric(omega=omega)/1e13')
# print(sim_test.source.flux_on_axis_theoric(omega=omega)/1e13)
#
# rad_theoric=sim_test.create_theoric_radiation()
# print(type(rad_theoric))
# print(sim_test.radiation.XY_are_like_in(rad_theoric))
# print('rad_theoric.intensity[0,0]')
# print(rad_theoric.intensity[0,0]/1e13)
# rad_theoric.plot()
#
#





#
# sim_test=create_simulation(magnetic_structure=und_test,electron_beam=beam_test,traj_method=TRAJECTORY_METHOD_ODE,
#                            rad_method=RADIATION_METHOD_APPROX_FARFIELD,Nb_pts_trajectory=10000)
#
# #print(sim_test.source.magnetic_field_strength())
# z=sim_test.trajectory.z*codata.c
# x=sim_test.trajectory.x*codata.c
#
# plt.plot(z,x)
# plt.xlabel('Z')
# plt.ylabel('X')
# plt.show()
#
# # observation_angle=np.linspace(0.0,sim_test.source.angle_wave_number(1,2),401)
# # sim_test.calculate_for_observation_angles(observation_angle=observation_angle)
#
# #sim_test.change_distance(D=None)
# # nb_period=np.linspace(10,50,5)
# # print(nb_period)
# # error=sim_test.error_radiation_method_nb_period(method=RADIATION_METHOD_APPROX,nb_period=nb_period)
# #
# # plt.plot(nb_period,error)
# # plt.show()
#
#
# # sim_test.plot_magnetic_field_along_Z()
# # sim_test.trajectory.plot()
# # sim_test.trajectory.plot_3D()
# # sim_test.calculate_until_wave_number(harmonic_number=1,wave_number=2)
# sim_test.radiation.plot()
# # sim_test.radiation.plot_wave(Nb_pts=sim_test.radiation_fact.Nb_pts)
#
# # X_max=sim_test.radiation.X.max()
# # Y_max=sim_test.radiation.Y.max()
#
# sim_test.plot_spectre_on_axis()
#sim_test.plot_spectre_central_cone()

