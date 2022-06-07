use self::config::{Config, DriveType};
use std::fs::{File, OpenOptions};
use std::env;

fn main(){
let c: Config = Config::default();
let out_dir = env::var_os("OUT_DIR").unwrap();
let dest_path = Path::new(&out_dir).join("config-desc.toml");
let config_desc: String = format!("
    # Default config and description of options\\
    hsize = {}    # size of the system\\
    ssize = {}      # size of subsystem (driven part)\\
    t = {}        # final time of simulation\\
    dt = {}      # time-step for simulation\\
    runs = {}       # number runs to perform\\
    threads = {}    # number of parallel threads\\
    trel = {}       # minus initial time\\
    tau = {}       # period of drive\\
    lambda = {}     # value of J_z coupling\\
    hfield = {:?}     # constant h field on entire system\\
    hs = {:?}     # constant h field on subsystem\\
    jvar = {}   # variance in J couplings (x, y, z independent)\\
    hvar = {}       # variance in field\\
    method = {}     # method for numerical integration (2=2nd order Suzuki-Trotter)\\
    ednsty = {} # energy-density of initial state\\
    file = \"{}\"   # pattern for log files i.e. log0.dat, log1.dat\\
    strob = {}  # stroboscopic evaluation\\
    offset = {}     # first file i.e. log0.dat\\
    beta = {}    # beta (determined from ednsty)\\
    e = {}    # eccentricity (scale of y-drive if drivetype=\"xyelliptic\"),
    drivetype = \"{:?}\" # type of driving, can be \"xyplane\", \"uniaxial\", \"xyelliptic\"
", c.hsize, c.ssize, c.t, c.dt, c.runs, c.threads, c.trel, c.tau, c.lambda, c.hfield, c.hs, c.jvar, c.hvar, c.method, c.ednsty, c.file, c.strob, c.offset, c.beta, c.e, c.drive);
let mut output = File::create(dest_path);
write!(output, &config_desc).unwrap();

}
