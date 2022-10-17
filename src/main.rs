use indicatif::{MultiProgress, ProgressBar, ProgressStyle};
use std::fs;
use std::fs::{File, OpenOptions};
use std::io::{prelude::*, BufReader};
use std::sync::mpsc;
use std::thread;
use threadpool::ThreadPool;

mod macros;

#[allow(unused_imports)]
use std::f64::consts::PI;

use sc::config::{Config, DriveType, CS};
use sc::spin::Spin;
use sc::spinchain::SpinChain;
use sc::spinchain_langevin::SpinChainLangevin;

mod writers;
use writers::{profile_in_bg, write_in_bg};

mod obs;
use obs::{
    average, find_beta, find_beta2, find_beta_maglab, find_beta_rot, find_betah, gen_hist,
    gen_hist_dynamics, gen_hist_magnus, gen_single_traj, get_obs, get_obs_maglab, get_obs_rot,
    get_obsh, run_dyn_profile, run_ext, run_langevin, run_mc, run_mc_magnus, run_mc_profile,
    run_sim, run_sim_response, run_tau, trajectory_mean_e,
};
//use obs::{estim_e, estim_e2, estim_e_maglab, estim_e_rot, estim_eh}

//use clap::{crate_version, App, Arg};
use clap::Parser;

/// Creates struct for different options to run to
#[derive(Parser, Debug)]
#[clap(author,version,about, long_about=None)]
struct Args {
    /// Override tau from config
    #[clap(short = 't', long = "tau")]
    tau: Option<f64>,
    /// Calculate average of runs
    #[clap(short = 'a', long = "avg", parse(from_flag))]
    avg: bool,
    ///Override config beta
    #[clap(short = 'b', long = "beta")]
    beta: Option<f64>,
    ///Calculate an average quantitiy in Monte-Carlo
    #[clap(short = 'm', long = "monte-carlo", parse(from_flag))]
    mc: bool,
    ///Calculate an average quantitiy in Monte-Carlo, spatially resolved
    #[clap(short = 'p', long = "monte-carlo-profile", parse(from_flag))]
    mc_profile: bool,
    ///Calculate an average quantitity in dynamical runs, spatially resolved
    #[clap(short = 'P', long = "dynamic-profile", parse(from_flag))]
    dyn_profile: bool,
    ///Calculate an average quantitity in Monte-Carlo, Magnus expansion to second order in omega^{-1}
    #[clap(short, long = "monte-carlo-magnus", parse(from_flag))]
    mc_full: bool,
    ///Generate histograms via Monte-Carlo (hist_mc.dat) and via time-evolution (hist_dyn.dat) (default POINTS=8000)
    #[clap(short = 'h', long = "histogram", value_name = "POINTS")]
    hist: Option<Option<usize>>,
    ///Generate time-evolution histogram (hist_dyn.dat) (default POINTS=1000)
    #[clap(short = 'd', long = "dynamic-histogram", value_name = "POINTS")]
    dynhist: Option<Option<usize>>,
    ///Generate time-evolution histogram (hist_sj.dat) (default POINTS=1000)
    #[clap(short = 's', long = "trajectory", value_name = "POINTS")]
    singletraj: Option<Option<usize>>,
    ///Directly simulate Langevin dynamics on the system proper (default GAMMA=1.0)
    #[clap(short = 'l', long = "langevin", value_name = "GAMMA")]
    langevin: Option<Option<f64>>,

    ///Generate histogram via Monte-Carlo (hist_mc_magnus.dat) for first-order Magnus expansion, and print averages (default POINTS=8000)
    #[clap(long = "magnus-hist", value_name = "POINTS")]
    magnus_hist: Option<Option<usize>>,
    ///Print description of config file
    #[clap(long = "config-desc")]
    config_desc: bool,

    ///Generate single-shot time-evolution at different drive periods tau
    #[clap(long = "n-steps", short, value_name = "TAU1,TAU2,...")]
    n_steps_taus: Option<Vec<f64>>,

    ///Find adiabatic ensemble H-\\omega S^z
    #[clap(long = "adiabatic", short = 'A', value_name = "TAU")]
    adiabatic: Option<f64>,

    /// Find adiabatic ensemble H-\\omega S^z inhomogeneous
    #[clap(long = "adiab_full", short = 'D', value_name = "TAU")]
    adiab2: Option<f64>,

    /// Fit temperature of Monte-Carlo ensemble with energy density E in rotating frame
    #[clap(
        long = "rot-frame",
        short,
        value_name = "E",
        allow_hyphen_values = true
    )]
    rot: Option<f64>,

    /// Produce profile of energy density as a function of t
    #[clap(short = 'e', long)]
    ext: bool,

    /// Fit temperature of Monte-Carlo ensemble with energy density E in lab frame, leading Magnus
    #[clap(short = 'E', long = "magnus-fit", value_name = "E")]
    maglab: Option<f64>, //E

    //Set step limit when doing n-steps
    #[clap(long)]
    steps: Option<u32>,

    ///Fit temperature of Monte Carlo ensemble with energy density E
    #[clap(short = 'f', allow_hyphen_values = true, value_name = "E")]
    fit: Option<f64>, //E

    ///Calculate effective ensemble temperature via conserved quantity arguments
    #[clap(short = 'F', allow_hyphen_values = true, value_name = "E")]
    fit_effective: Option<f64>, //E

    ///Calculate response function
    #[clap(long)]
    response: bool, //E

    ///Calculate effective ensemble with initial temp and first-order Magnus
    #[clap(short = 'H', long, allow_hyphen_values = true, value_name = "E")]
    high_freq: Option<f64>,
}

use glob::glob;

fn print_config_description() {
    println!("{}", CS);
}

fn main() {
    let args = Args::parse();

    let mut default = true;

    //Want to read num from file
    let mut conf: Config = SpinChain::read_config("config.toml");
    if let Some(tau) = args.tau {
        conf.tau = tau;
    }

    match args.config_desc {
        false => {}
        true => {
            print_config_description();
            default = false;
        }
    }
    if args.ext {
        println!("Generating ext");
        let mut conf = conf.clone();
        run_ext(&mut conf);
        default = false;
    }

    let steps = args.steps.unwrap_or(1000);

    if let Some(taus) = args.n_steps_taus {
        //        println!("Overriding to: {}", b);
        /*         let tau_vec: Vec<f64> = taus
            .split(',')
            .map(|word| word.parse::<f64>().unwrap())
            .collect();
        println!("{:?}", tau_vec);
        */
        run_tau(taus, steps, &mut conf);
        default = false;
    }

    if let Some(tau) = args.adiab2 {
        println!("Adiab 2");
        conf.tau = tau;

        conf.hs = vec![1.0, 0.0, 0.0];
        conf.hfield = vec![0.0, 0.0, -2.0 * PI / conf.tau];

        let result: f64 = find_beta2(&mut conf);

        println!("Optimal beta: {}", result);

        let beta_file = OpenOptions::new()
            .create(true)
            .append(true)
            .open("beta.dat")
            .unwrap();
        writeln!(&beta_file, "{} {} {}", conf.tau, conf.beta, result).unwrap();

        default = false;
    }
    if let Some(targete) = args.rot {
        println!("Evaluating with first-order correction in rot frame");
        //Start with no extra fields paplied to chain
        conf.ednsty = targete;

        // find_beta_rot does Monte-Carlo updates with additional fields, and also evaluates energy
        // with additional fields
        let result: f64 = find_beta_rot(&mut conf);
        conf.beta = result;

        let betat_file = OpenOptions::new()
            .create(true)
            .append(true)
            .open("betat-rot.dat")
            .unwrap();

        let results: [f64; 4] = get_obs_rot(&mut conf);

        writeln!(
            &betat_file,
            "{} {} {} {} {} {} {}",
            conf.tau, targete, conf.beta, results[0], results[1], results[2], results[3]
        )
        .unwrap();

        default = false;
    }
    if let Some(targete) = args.maglab {
        println!("Evaluating with first-order correction in lab frame");
        //Start with no extra fields paplied to chain
        conf.ednsty = targete;

        // find_beta_rot does Monte-Carlo updates with additional fields, and also evaluates energy
        // with additional fields
        let result: f64 = find_beta_maglab(&mut conf);
        conf.beta = result;

        let betat_file = OpenOptions::new()
            .create(true)
            .append(true)
            .open("betat-maglab.dat")
            .unwrap();

        let results: [f64; 4] = get_obs_maglab(&mut conf);

        writeln!(
            &betat_file,
            "{} {} {} {} {} {} {}",
            conf.tau, targete, conf.beta, results[0], results[1], results[2], results[3]
        )
        .unwrap();

        default = false;
    }

    if let Some(targete) = args.fit {
        conf.ednsty = targete;
        conf.hs = vec![1.0, 0.0, 0.0];
        conf.hfield = vec![0.0, 0.0, -2.0 * PI / conf.tau];
        conf.hsize = conf.ssize;

        let result: f64 = find_beta(&mut conf);
        conf.beta = result;

        //Also calculate observables from this

        let results: [f64; 4] = get_obs(&mut conf);

        let betat_file = OpenOptions::new()
            .create(true)
            .append(true)
            .open("betat.dat")
            .unwrap();
        writeln!(
            &betat_file,
            "{} {} {} {} {} {} {}",
            conf.tau, targete, conf.beta, results[0], results[1], results[2], results[3]
        )
        .unwrap();

        default = false;
    }

    if let Some(targete) = args.fit_effective {
        conf.ednsty = targete;
        conf.hs = vec![1.0, 0.0, 0.0];
        conf.hfield = vec![0.0, 0.0, -2.0 * PI / conf.tau];
        //        conf.hsize = conf.ssize;

        let result: f64 = find_betah(&mut conf);
        conf.beta = result;

        //Also calculate observables from this
        conf.hs = vec![1.0, 0.0, 0.0];

        let results: [f64; 4] = get_obsh(&mut conf);

        let betat_file = OpenOptions::new()
            .create(true)
            .append(true)
            .open("betat.dat")
            .unwrap();
        writeln!(
            &betat_file,
            "{} {} {} {} {} {} {}",
            conf.tau, targete, conf.beta, results[0], results[1], results[2], results[3]
        )
        .unwrap();

        default = false;
    }
    if let Some(targete) = args.high_freq {
        conf.ednsty = targete;

        conf.hfield = vec![0.0, 0.0, -conf.tau / (4.0 * PI)];
        //        conf.hsize = conf.ssize;

        let results: [f64; 4] = get_obsh(&mut conf);

        let betat_file = OpenOptions::new()
            .create(true)
            .append(true)
            .open("betat.dat")
            .unwrap();
        writeln!(
            &betat_file,
            "{} {} {} {} {} {} {}",
            conf.tau, targete, conf.beta, results[0], results[1], results[2], results[3]
        )
        .unwrap();

        default = false;
    }

    if let Some(tau) = args.adiabatic {
        //Want to search for beta such that
        // < H_f e^{-\beta H_f} > = e_0
        // Can't think of a clever way to do this,
        // so will set up a grid of beta, and produce estimates for
        // the mean energy. Can then do a bisection search
        conf.tau = tau;
        println!("tau={}", conf.tau);

        // First need to sample a single trajectory to get the target mean energy density
        // Collect 250 sample points from 4 separate trajectories to estimate the mean energy
        // density

        let mean_e: f64 = trajectory_mean_e(&mut conf.clone());
        println!("Mean energy: {}", mean_e);

        // Now set this as the target
        // and calculate the optimal effective beta
        conf.ednsty = mean_e;

        conf.hs = vec![1.0, 0.0, -2.0 * PI / conf.tau];
        //        conf.hs = vec![0.0, 0.0, 0.0];

        let initial_beta: f64 = conf.beta;
        let result: f64 = find_beta(&mut conf);
        println!("Optimal beta: {}", result);

        //Write this to a file

        conf.beta = result;

        // Get e, mx, my, mz for this beta
        let results: [f64; 4] = get_obs(&mut conf);

        let beta_file = OpenOptions::new()
            .create(true)
            .append(true)
            .open("beta.dat")
            .unwrap();
        writeln!(
            &beta_file,
            "{} {} {} {} {} {} {} {}",
            conf.tau, initial_beta, result, results[0], results[1], results[2], results[3], mean_e
        )
        .unwrap();

        default = false;
    }

    if let Some(beta) = args.beta {
        //        println!("Overriding to: {}", b);
        conf.beta = beta;
    } else {
        conf.beta = SpinChain::solve_beta(conf.ednsty);
    }

    if let Some(points) = args.hist {
        let points = points.unwrap_or(8000_usize);

        println!("Running hist");
        default = false;
        gen_hist(&mut conf, points);
    }

    if let Some(sample_num) = args.dynhist {
        let sample_num: usize = sample_num.unwrap_or(1000_usize);
        println!("Running dynamic hist");
        default = false;
        gen_hist_dynamics(&mut conf, sample_num);
    }

    if let Some(points) = args.magnus_hist {
        let points = points.unwrap_or(8000_usize);
        println!("Running Monte-Carlo for first-order Magnus expansion");
        default = false;
        gen_hist_magnus(&mut conf, points);
    }

    if let Some(single_samples) = args.singletraj {
        let single_samples = single_samples.unwrap_or(1000_usize);
        println!("Running single trajectory");
        default = false;
        gen_single_traj(&mut conf, single_samples);
    }

    if args.avg {
        default = false;
        average(&mut conf);
    }

    if args.mc {
        default = false;
        run_mc(&mut conf);
    }

    if args.mc_profile {
        default = false;
        run_mc_profile(&mut conf);
    }

    if args.dyn_profile {
        default = false;
        run_dyn_profile(&mut conf);
    }

    if args.mc_full {
        default = false;
        run_mc_magnus(&mut conf);
    }

    if let Some(langevin) = args.langevin {
        let gamma: f64 = langevin.unwrap_or(1.0);
        println!(
            "Running Langevin dynamics on the system proper with
                gamma={}",
            gamma
        );
        conf.hsize = conf.ssize;
        run_langevin(&mut conf, gamma);
        default = false;
    }

    match args.response {
        true => {
            run_sim_response(&mut conf);
            default = false;
        }
        false => (),
    }

    match default {
        true => run_sim(&mut conf),
        false => (),
    }
}
