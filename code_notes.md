# code notes — phase 1

lecture code = two files.
- LJ_gas.py — all the physics (classes + functions)
- LJ_gas_run_MD.py — the script that runs a simulation and makes the plots

## where epsilon and sigma are set
LJ_gas_run_MD.py, ~line 72
- sigma_argon = 0.34 nm (distance where the potential crosses zero)
- epsilon_argon = 120*R*1e-3 kJ/mol (depth of the attractive well)
- stored per particle with ps.set_parameters(), but every particle gets
  the same values, so it's one kind of atom (argon)
- minimum of the potential sits at r_min = 2^(1/6)*sigma ≈ 1.122*sigma

## LJ potential — potential_energy(), LJ_gas.py line 146
- computes all pairwise distances at once with numpy broadcasting
- applies minimum image convention for periodic boundaries:
  rij_matrix -= L * np.rint(rij_matrix / L)
- takes only unique pairs (upper triangle, i < j)
- clips tiny distances to rij_min so (sigma/r)^12 doesn't explode
- the formula: sr6 = (sigma/r)^6, then E = 4*eps*(sr6^2 - sr6)
  that IS the 12-6 potential, sr6^2 = (sigma/r)^12
- returns total E_pot in kJ/mol

## forces — calculate_force(), LJ_gas.py line 275
- same distance setup as potential_energy
- dV_dr = 24*eps/r * (-2*sr6^2 + sr6)  <- derivative of the potential
- force = negative gradient, pointed along the vector between the pair
- newton's third law: each pair force goes -f to particle i, +f to particle j
- updates ps.force in place, doesn't return anything

## integrator — built from three pieces, LJ_gas.py
- A_step (line 338): move positions. r += v*dt
- B_step (line 363): kick velocities with the force. v += (F/m)*dt
- O_step (line 391): the langevin/thermostat piece (see below)
- simulate_NVE_step (line 430) = B A B = velocity verlet. energy conserved.
- simulate_NVT_step (line 460) = B A O A B = langevin integrator.
  temperature controlled. this is the one the run script uses (NVT = True)

## thermostat — O_step, LJ_gas.py line 391
- langevin thermostat. two effects on velocities every step:
  1. friction: v gets multiplied by exp(-xi*dt), slows everything down
  2. noise: gaussian random kicks, scaled to the TARGET temperature
- together they drag the system toward sim.temperature
- xi = 1/tau_thermostat (set in SimulationParameters, line 74)
- run script uses tau_thermostat = 1 ps. small tau = strong coupling = fast
  relaxation to target T
- IMPORTANT FOR PROJECT: cooling = lowering sim.temperature. the O_step is
  what actually pulls the system down to the new temperature

## energies stored per step
LJ_gas_run_MD.py main loop, ~line 150
- energy_trajectory column 0 = E_pot, 1 = E_kin, 2 = T, 3 = P
- recorded every single step, plus positions in position_trajectory
- at the end saved as .npy (binary) and .dat (text), plots saved as .png

## initial setup
LJ_gas.py line 108 + 112, values in LJ_gas_run_MD.py ~line 70
- positions: uniform random in a cubic box (initialize_positions)
- velocities: maxwell-boltzmann at temperature T, then center-of-mass
  motion removed so the box doesn't drift (initialize_velocities)
- defaults: N = 200 particles, box = 100 nm, dt = 0.1 ps, 1000 steps, T = 300 K

## periodic boundary conditions
- yes, used in two places:
  1. minimum image convention inside potential_energy and calculate_force
     (interactions wrap around the box edges)
  2. apply_periodic_boundary (line 497), called at the end of every step,
     wraps positions back into [0, L) with np.mod
- so a particle leaving the right side comes back in on the left

## observation from first run
- default box is SO dilute (200 particles in 100^3 nm^3) that E_pot is
  basically 0 the whole time. spikes in the E_pot plot = rare events where
  two particles get close. dips = attractive well, big positive spike =
  head-on collision hitting the repulsive wall
- for condensation we probably need a smaller box / higher density

## main loop summary (maggie)
- every step the code does one round of BAOAB (simulate_NVT_step):

1. B: velocities get a half kick from the current forces. v += (F/m) * dt/2
2. A: positions move half a step using those velocities. r += v * dt/2
3. O: thermostat. velocities get shrunk by exp(-xi*dt) (friction) and
   gaussian noise gets added, scaled so the system heads toward the
   target temperature. this is the langevin part. only randomness in the loop.
4. A: positions move the second half step
5. forces recalculated from the new positions (calculate_force)
6. B: second half kick with the new forces
7. positions wrapped back into the box (periodic boundaries)

then back in LJ_gas_run_MD.py the new positions go into position_trajectory
and E_pot, E_kin, T, P get computed and stored in energy_trajectory.
repeat n_steps times.

short version: kick, move, thermostat shuffles velocities toward target T,
move, new forces, kick, wrap. without step 3 it would just be velocity
verlet (NVE) and energy would be conserved instead of T being controlled.
## main loop summary (soner)
- (write yours here, 3-4 sentences, own words)



## tunable parameters (defaults from baseline run)
all set in LJ_gas_run_MD.py, parameters section ~line 68-85

- n_particles = 200
- mass_argon = 39.95 u
- sigma_argon = 0.34 nm         <- LJ zero crossing
- epsilon_argon = 120*R*1e-3 ~ 0.998 kJ/mol   <- LJ well depth
- dt = 0.1 ps                   <- timestep
- n_steps = 1000                <- 100 ps total
- temperature = 300 K           <- thermostat target
- box_length = 100 nm           <- controls density. density = N and box together
- tau_thermostat = 1 ps         <- coupling. small = faster relaxation to target T
- rij_min = 0.01 nm             <- numerical safety cutoff, not physics
- NVT = True                    <- True = langevin thermostat, False = NVE
- file_name_base                <- output file names

## baseline run stats
- density 1.327e-05 g/cm^3. liquid argon is ~1.4 g/cm^3, so we are 100000x
  too dilute. explains the boring E_pot
- runtime 359 s for 1000 steps (0.36 s/step). bottleneck = python for loop
  in calculate_force

  ## reduced units
- LJ papers use reduced units: energy in eps, length in sigma, T* = kB*T/eps
- our eps corresponds to eps/kB = 120 K (thats where the 120 in the run script
  comes from). so T* = T_kelvin / 120
- baseline 300 K = T* 2.5 (hot gas)
- hot run 600 K = T* 5.0
- cold run 30 K = T* 0.25
- useful reference points from the literature: LJ triple point T* ~ 0.69,
  critical point T* ~ 1.3. so condensation experiments should target
  T* below ~1, i.e. below ~120-160 K for our argon numbers
- lengths: box = 100 nm = ~294 sigma. enormous. liquid-ish density would be
  more like box ~ 6-7 sigma for 200 particles