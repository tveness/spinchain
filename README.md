# Classical spin chain simulation

![Build](https://github.com/tveness/spinchain/workflows/build-test/badge.svg)

Classical spin-chain simulation for project investigating periodic driving of a
spin chain coupled to a bath.

Numerically simulates the time-evolution of the classical spin-chain given by
$$
H = - \sum_j  {\bf S_j} J_j {\bf S_{j+1}} + \sum_j {\bf B_j} \cdot {\bf S_j}
$$
where $J={\rm diag}(1,1,\lambda)$ plus noise on the diagonal drawn from a
normal distribution, and the local magnetic field is a static field plus noise
drawn from a separate normal distribution.


# Configuration features

In config.toml (if not present, default will be generated when running), there
are the following parameters
```
hsize: is the size of the system
ssize: size of subsystem (driven part)
t: total time of simulation (in units of J)
dt: time-step of numerical evolution
runs: how many parallel runs to simulate
trel: time offset (start at t=-trel)
tau: period of driving
lambda: anisotropy of J-coupling
hfield: static field applied to all sites
jvar: variance in J-coupling
hvar: variance in static magnetic field
method: method for Suzuki-Trotter
ednsty: target energy density for chain initialisation
file: logging filename e.g. "log" would send runs to "log0.dat", "log1.dat", ...
strob: evaluate stroboscopically? true/false
offset: for 0, start logging at file "log0.dat"
drive: drive on/off (true/false)
pub beta: inverse temperature for Monte Carlo
```
