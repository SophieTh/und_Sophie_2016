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
from pySRU.RadiationFactory import RADIATION_METHOD_NEAR_FIELD, RADIATION_METHOD_APPROX_FARFIELD


eV_to_J=1.602176487e-19
######################################################

#BM_test=MagneticStructureBendingMagnet(Bo=0.8, div=5e-3, R=25.0)

E_1=7876.0

beam_test=ElectronBeam(Electron_energy=1.3, I_current=1.0)
beam_ESRF=ElectronBeam(Electron_energy=6.0, I_current=0.2)
und_test=Undulator(  K = 1.87,  period_length= 0.035, length=0.035 * 14)
ESRF18=Undulator( K = 1.68, period_length = 0.018, length=2.0)
ESRFBM=BM(Bo=0.8,length=0.1249994791673177)



vx= 2e-4
vz= np.sqrt(beam_test.electron_speed()**2-vx**2)*codata.c

initial_cond=np.array([ vx*codata.c,  0.00000000e+00 ,vz , 0.0 , 0.0 ,-0.42,])
X=np.linspace(-0.1,0.1,100)
Y=np.linspace(-0.1,0.1,100)
sim_test = create_simulation(magnetic_structure=und_test, electron_beam=beam_test, traj_method=TRAJECTORY_METHOD_ODE,
                             rad_method=RADIATION_METHOD_APPROX_FARFIELD, Nb_pts_trajectory=10000,distance=100,
                             initial_condition=initial_cond,X=X,Y=Y)


sim_test.print_parameters()
# sim_test.plot_everything()

sim_test.trajectory.plot_3D()
sim_test.radiation.plot()

