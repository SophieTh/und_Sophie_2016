import numpy as np
import scipy.constants as codata
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import scipy.integrate as integrate
import time
from pySRU.MagneticStructureUndulatorPlane import MagneticStructureUndulatorPlane as Undulator
from pySRU.MagneticStructureBendingMagnet import MagneticStructureBendingMagnet  as BM
from pySRU.ElectronBeam import ElectronBeam
from pySRU.Source import Source
from pySRU.Simulation import Simulation ,create_simulation
from pySRU.TrajectoryFactory import TRAJECTORY_METHOD_ANALYTIC,TRAJECTORY_METHOD_ODE,\
                                        TRAJECTORY_METHOD_INTEGRATION
from pySRU.RadiationFactory import RADIATION_METHOD_NEAR_FIELD, RADIATION_METHOD_APPROX,\
                                RADIATION_METHOD_FARFIELD, RADIATION_METHOD_APPROX_FARFIELD


eV_to_J=1.602176487e-19
######################################################

#BM_test=MagneticStructureBendingMagnet(Bo=0.8, div=5e-3, R=25.0)

E_1=7876.0

#omega=E_1*eV_to_J/codata.hbar

# # recuperation des donnees de B en array en fonction de z
# reference=np.load("../tests/x_ray_booklet_field.npz")
# Z=reference['ct']
# Z -= (Z[-1])/2.0
# By=reference['B_y']


#
# # # minimum  :
#
beam_test=ElectronBeam(Electron_energy=1.3, I_current=1.0)
beam_ESRF=ElectronBeam(Electron_energy=6.0, I_current=0.2)
und_test=Undulator(  K = 1.87,  period_length= 0.035, length=0.035 * 14)
ESRF18=Undulator( K = 1.68, period_length = 0.018, length=2.0)
ESRFBM=BM(Bo=0.8,L=0.1249994791673177)


# TODO Tout metre en GeV !
# essai=(3./2.)*(codata.h/(2.*np.pi))*((E/(codata.m_e*codata.c))**2)*B*codata.e
# print(essai)
#
# print(3./(codata.m_e**3*codata.c**5*3.3))
# print(5.59e3/2.9e9)
#
# print(codata.m_e*codata.c/(np.sqrt(codata.h*codata.e)))



#
# #
# sim_test = create_simulation(magnetic_structure=ESRFBM, electron_beam=beam_ESRF,
#                              traj_method=TRAJECTORY_METHOD_ANALYTIC,
#                              distance=100.,X=X,Y=Y,Nb_pts_trajectory=10000)
#
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
# # # TODO peut etre faire deux classe source plutot que magne_structure .....yes
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
# #TODO simps ou trapz pour rad ?
# #TODO faire varier erreur rad avec method de radiation et nb period augmente






sim_test=create_simulation(magnetic_structure=und_test,electron_beam=beam_test,traj_method=TRAJECTORY_METHOD_ANALYTIC,
                           rad_method=RADIATION_METHOD_APPROX_FARFIELD)

# sim_test.plot_magnetic_field_along_Z()
# sim_test.trajectory.plot()
# sim_test.trajectory.plot_3D()
#sim_test.radiation.plot()

X_max=sim_test.radiation.X.max()
Y_max=sim_test.radiation.Y.max()

X=np.linspace(-X_max,X_max,101)
Y=np.linspace(-Y_max,Y_max,105)
sim_test.change_XY_radiation(X=X,Y=Y)


sim_test.radiation.plot()

