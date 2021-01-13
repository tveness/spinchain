# Classical spin chain simulation

Linux binary [available](https://github.com/tveness/spinchain/releases)

![Build](https://github.com/tveness/spinchain/workflows/build-test/badge.svg)

## How to build the rust project

1. Clone the repository `git clone https://github.com/tveness/spinchain.git`
2. Ensure the rust toolchain is installed, see https://rustup.rs/
3. Run `cargo build --release` to compile a fast version of the code, it may
   take a couple of minutes
4. The binary is now located at target/release/sc
5. `sc --help` will product the help dialogue listed below
6. `sc` will run a simple chain, and write the default configuration to `config.toml`



Classical spin-chain simulation for project investigating periodic driving of a
spin chain coupled to a bath.

Numerically simulates the time-evolution of the classical spin-chain given by
$$
H = - \sum_j  {\bf S_j} J_j {\bf S_{j+1}} + \sum_j {\bf B_j} \cdot {\bf S_j}
$$
where $J={\rm diag}(1,1,\lambda)$ plus noise on the diagonal drawn from a
normal distribution, and the local magnetic field is a static field plus noise
drawn from a separate normal distribution.

## Usage (output of `sc --help`)
```
spinchain 0.1.21
Thomas Veness <thomas.veness@nottingham.ac.uk>
Run classical spin chain simulation

USAGE:
    sc [FLAGS] [OPTIONS]

FLAGS:
    -a, --avg            Calculate average of runs
        --config-desc    Print description of config file
        --help           Prints help information
    -m, --monte-carlo    Calculate an average quantitity in Monte-Carlo
    -V, --version        Prints version information

OPTIONS:
    -b, --beta <BETA>             Overrides config beta
    -h, --histogram <POINTS>      Generate histograms via Monte-Carlo (hist_mc.dat) and via time-evolution
                                  (hist_dyn.dat) [default: 8000]
        --magnus-hist <POINTS>    Generate histogram via Monte-Carlo (hist_mc_magnus.dat) for first-order Magnus
                                  expansion, and print averages [default: 8000]
```



# Configuration features

In config.toml (if not present, default will be generated when running), there
are the following parameters  (obtains by running `sc --config-desc`)
```
# Default config and description of options
hsize = 512    # size of the system
ssize = 8      # size of subsystem (driven part)
t = 256        # final time of simulation
dt = 0.02      # time-step for simulation
runs = 2       # number runs to perform
threads = 2    # number of parallel threads
trel = 0       # minus initial time
tau = 10       # period of drive
lambda = 1     # value of J_z coupling
hfield = [0.0, 0.0, 0.0] # constant h field on entire system
hs = [0.0, 0.0, 0.0]     # constant h field on subsystem
jvar = 0.001   # variance in J couplings (x, y, z independent)
hvar = 0       # variance in field
method = 2     # method for numerical integration (2=2nd order Suzuki-Trotter)
ednsty = -0.66 # energy-density of initial state
file = "log"   # pattern for log files i.e. log0.dat, log1.dat
strob = false  # stroboscopic evaluation
offset = 0     # first file i.e. log0.dat
drive = true   # drive enabled
beta = 2.88    # beta (determined from ednsty)
drivetype = "xyplane" # type of driving, can be "xyplane", "uniaxial", "none"
```

# Running the program

By default when run, the program will produce a set of time-evolved data from
an initial state.

It will log energy density of entire system, and then energy density for the
subsystem, as well as subsystem magnetisations.

# Histograms

Comparing histograms for Monte Carlo and dynamical evolution to check that the system thermalises correctly.
Plots here for beta = 2.89, and B-field = [1.0,0.0,0.0]

![Energy density](examples/hist_e.png)
![M_x density](examples/hist_mx.png)
![M_y density](examples/hist_my.png)
![M_z density](examples/hist_mz.png)


# To do

Write Monte Carlo for first order Magnus expansion result


