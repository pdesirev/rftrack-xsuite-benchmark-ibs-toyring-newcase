import numpy as np
from cpymad.madx import Madx
import xpart as xp
import xtrack as xt
import xfields as xf
import xcoll as xc
import matplotlib.pyplot as plt
import xobjects as xo

mad = Madx()

mad.input('''!===========================================================
! 16-Cell Geometrically Closed FODO Ring — 100 MeV electrons
! Fully flattened for MAD-X 5.09
!===========================================================

TITLE, "16-cell FODO ring — flattened, geometrically closed";

!-----------------------------
! Beam
!-----------------------------
BEAM, PARTICLE=electron, ENERGY=0.5,
EX=3.11E-3, EY=3.11E-3, SIGT=1.0E-3, SIGE=1.0E-3;
deg = pi/180.0;

!-----------------------------
! Elements
!-----------------------------
! Focusing quadrupoles
kqf := 2.12;
kqd := -1.666;
QF: QUADRUPOLE, L=0.3, K1=kqf;
QD: QUADRUPOLE, L=0.3, K1=kqd;
D1: DRIFT, L=1.0;
DC: DRIFT, L=0.8;
B1: SBEND, L=2.4, ANGLE=22.5*deg, E1=0, E2=0;
CAV: RFCAVITY, L=0.2, VOLT=1.0, FREQ=3.74, HARMON=1, LAG=0.5;

!-----------------------------
! Flattened 16-cell ring
!-----------------------------
RING: SEQUENCE, L=80.0, REFER=ENTRY;

! Cell 1
  QF, AT=0.0;
  D1, AT=0.3;
  QD, AT=1.3;
  D1, AT=1.6;
  B1, AT=2.6;

! Cell 2
  QF, AT=5.0;
  D1, AT=5.3;
  QD, AT=6.3;
  D1, AT=6.6;
  B1, AT=7.6;

! Cell 3
  QF, AT=10.0;
  D1, AT=10.3;
  QD, AT=11.3;
  D1, AT=11.6;
  B1, AT=12.6;

! Cell 4
  QF, AT=15.0;
  D1, AT=15.3;
  QD, AT=16.3;
  D1, AT=16.6;
  B1, AT=17.6;

! Cell 5
  QF, AT=20.0;
  D1, AT=20.3;
  QD, AT=21.3;
  D1, AT=21.6;
  B1, AT=22.6;

! Cell 6
  QF, AT=25.0;
  D1, AT=25.3;
  QD, AT=26.3;
  D1, AT=26.6;
  B1, AT=27.6;

! Cell 7
  QF, AT=30.0;
  D1, AT=30.3;
  QD, AT=31.3;
  D1, AT=31.6;
  B1, AT=32.6;

! Cell 8
  QF, AT=35.0;
  D1, AT=35.3;
  QD, AT=36.3;
  D1, AT=36.6;
  B1, AT=37.6;

! Cell 9
  QF, AT=40.0;
  D1, AT=40.3;
  QD, AT=41.3;
  D1, AT=41.6;
  B1, AT=42.6;

! Cell 10
  QF, AT=45.0;
  D1, AT=45.3;
  QD, AT=46.3;
  D1, AT=46.6;
  B1, AT=47.6;

! Cell 11
  QF, AT=50.0;
  D1, AT=50.3;
  QD, AT=51.3;
  D1, AT=51.6;
  B1, AT=52.6;

! Cell 12
  QF, AT=55.0;
  D1, AT=55.3;
  QD, AT=56.3;
  D1, AT=56.6;
  B1, AT=57.6;

! Cell 13
  QF, AT=60.0;
  D1, AT=60.3;
  QD, AT=61.3;
  D1, AT=61.6;
  B1, AT=62.6;

! Cell 14
  QF, AT=65.0;
  D1, AT=65.3;
  QD, AT=66.3;
  D1, AT=66.6;
  B1, AT=67.6;

! Cell 15
  QF, AT=70.0;
  D1, AT=70.3;
  QD, AT=71.3;
  D1, AT=71.6;
  B1, AT=72.6;

! Cell 16
  QF, AT=75.0;
  D1, AT=75.3;
  QD, AT=76.3;
  DC, AT=76.6;
  CAV,AT=77.4;
  B1, AT=77.6;



ENDSEQUENCE;

!-----------------------------
! Activate the ring
!-----------------------------
USE, SEQUENCE=RING;

!-----------------------------
! Match periodic optics
!-----------------------------
MATCH, SEQUENCE=RING;
  VARY, NAME=kqf, STEP=0.001;
  VARY, NAME=kqd, STEP=0.001;

  ! enforce periodic alpha + dispersion
  !CONSTRAINT, RANGE=#S/#E, ALFX=0., ALFY=0., DX=0., DPX=0., TUNEX=3.09857, TUNEY=2.47505;
  !CONSTRAINT, SEQUENCE=RING, RANGE=#E, BETX=6.6, ALFX=-0.8, BETY=3.5, ALFY=0.2;
  CONSTRAINT, SEQUENCE=RING, RANGE=#E, BETX=6.6, BETY=3.5;
  LMDIF, CALLS=1500, TOLERANCE=1.0E-12;
ENDMATCH;

PRINT, TEXT="After Step 1: FMINS =", VALUE=FMINS;
PRINT, TEXT="K1(QF) =", VALUE=QF->K1;
PRINT, TEXT="K1(QD) =", VALUE=QD->K1;
TWISS, FILE="twiss_withCavity2.tws";
''')
line = xt.Line.from_madx_sequence(
    sequence=mad.sequence.ring,
    allow_thick=True,
    enable_align_errors=True,
    deferred_expressions=True,
)
context = xo.ContextCpu()

line.set_particle_ref('electron', p0c=500e6) # the p0c is in eV
tw = line.twiss(method="6d")

# Emittance monitors
monitor = xc.EmittanceMonitor.install(
    line=line,
      name="monitor",  # Create a unique name for each monitor
      at_s= 0.001,       # Position scaled by the index
      stop_at_turn=1000,        # Number of turns to track
)

monitor2 = xt.ParticlesMonitor(
    start_at_turn=1,
    stop_at_turn=1000,
    num_particles=10000,
    auto_to_numpy=True
)

ibs_kick = xf.IBSKineticKick(num_slices=100)
line.configure_intrabeam_scattering(
    element=ibs_kick, name="ibskick", index=1, update_every=5
)
line.build_tracker(_context=context)
import json
# Load particles from json file to selected context
with open("particles_25pC.json", 'r') as fid:
    particles= xp.Particles.from_dict(json.load(fid), _context=context)

line.track(particles, num_turns=1000, with_progress=1, turn_by_turn_monitor=True)


nemitt_x_values = monitor.nemitt_x
nemitt_y_values = monitor.nemitt_y
nemitt_z_values = monitor.nemitt_zeta
gemitt_z_values = monitor.gemitt_zeta
at_s_points = [i for i in range(100)]
nemitt_x_values = np.array(nemitt_x_values)
nemitt_y_values = np.array(nemitt_y_values)
nemitt_z_values = np.array(nemitt_z_values)
gemitt_z_values = np.array(gemitt_z_values)

rec = line.record_last_track
x = rec.x
px = rec.px
y = rec.y
py = rec.py
zeta = rec.zeta
beta0 = rec.beta0
delta = rec.delta
p0c = rec.p0c
for turn in range(1000):
  N = x.shape[0]
  P = 1.0 + delta[:, turn]
  Pz = np.sqrt(P**2 - px[:, turn]**2 - py[:, turn]**2)
  phase_space = np.zeros((N, 9))
  phase_space[:, 0] = x[:, turn] * 1e3
  phase_space[:, 1] = px[:, turn] * 1e3 / Pz
  phase_space[:, 2] = y[:, turn] * 1e3
  phase_space[:, 3] = py[:, turn] * 1e3 / Pz
  phase_space[:, 4] = -zeta[:, turn] * 1e3 / beta0[:, turn]
  phase_space[:, 5] = p0c[:, turn] * P / 1e6
  np.savetxt(
        f"Results/PS_XSuite_track/ps_{turn+1}.dat",
        phase_space,
        fmt="%.6e"
    )
    
line.discard_tracker()

file_path_x = "Results/XSuite_25pC_nemitt_x_1000.dat"
file_path_y = "Results/XSuite_25pC_nemitt_y_1000.dat"
file_path_z = "Results/XSuite_25pC_nemitt_z_1000.dat"
file_path_gz = "Results/XSuite_25pC_gemitt_z_1000.dat"

#np.savetxt(file_path_x, nemitt_x_values, fmt="%.6e")  # Adjust format if needed
#np.savetxt(file_path_y, nemitt_y_values, fmt="%.6e")  # Adjust format if needed
#np.savetxt(file_path_z, nemitt_z_values, fmt="%.6e")  # Adjust format if needed
#np.savetxt(file_path_gz, gemitt_z_values, fmt="%.6e")  # Adjust format if needed
