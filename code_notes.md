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
- (write yours here, 3-4 sentences, own words)

## main loop summary (soner)
- (write yours here, 3-4 sentences, own words)