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
from pySRU.Source import Source
from pySRU.Simulation import Simulation ,create_simulation
from pySRU.TrajectoryFactory import TRAJECTORY_METHOD_ANALYTIC,TRAJECTORY_METHOD_ODE
from pySRU.RadiationFactory import RADIATION_METHOD_NEAR_FIELD, RADIATION_METHOD_APPROX_FARFIELD,RADIATION_METHOD_APPROX


eV_to_J=1.602176487e-19
######################################################

#BM_test=MagneticStructureBendingMagnet(Bo=0.8, div=5e-3, R=25.0)

E_1=7876.0

beam_test=ElectronBeam(Electron_energy=1.3, I_current=1.0)
beam_ESRF=ElectronBeam(Electron_energy=6.0, I_current=0.2)
und_test=Undulator(  K = 1.87,  period_length= 0.035, length=0.035 * 14)
ESRF18=Undulator( K = 1.68, period_length = 0.018, length=2.0)
ESRFBM=BM(Bo=0.8,horizontale_divergeance=0.005,electron_energy=6.0)



vx= 2e-4
vz= np.sqrt(beam_test.electron_speed()**2-vx**2)*codata.c

# initial_cond=np.array([ vx*codata.c,  0.00000000e+00 ,vz , 0.0 , 0.0 ,-0.42,])
#X=np.linspace(-0.02,0.02,150)
#Y=np.linspace(-0.02,0.02,150)
sim_test = create_simulation(magnetic_structure=ESRF18, electron_beam=beam_ESRF, traj_method=TRAJECTORY_METHOD_ANALYTIC,
                             rad_method=RADIATION_METHOD_APPROX_FARFIELD, Nb_pts_trajectory=5000,distance=20)
#sim_test.change_energy_eV(E=20000)
#sim_test.change_harmonic_number(5.5)
#print(sim_test.source.choose_nb_pts_trajectory())
#print(sim_test.source.choose_distance_automatic())
sim_test.print_parameters()


sim_test.trajectory.plot_3D(title="Analytical electron trajectory in bending magnet")
#sim_test.trajectory.plot(title="Analytical electron trajectory in undulator")
#sim_test.calculate_on_central_cone()
#sim_test.calculate_spectrum_on_axis()
# omega1=sim_test.source.harmonic_frequency(1)
# omega_array=np.arange(.9*omega1,3.06*omega1,0.01*omega1)
# sim_test.calculate_spectrum_central_cone(abscissas_array=omega_array)
# sim_test.change_harmonic_number(3.)
print (sim_test.radiation.max())
# sim_test.radiation.plot(title="radiation intensity from ODE trajectory solution ")
sim_test.radiation.plot(title="radiation intensity , far field formula")