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

use sc::config::{cs, Config, DriveType};
use sc::spin::Spin;
use sc::spinchain::SpinChain;
use sc::spinchain_langevin::SpinChainLangevin;

mod writers;
use writers::{profile_in_bg, write_in_bg};

mod obs;
use obs::{
    average, estim_e, estim_e2, estim_e_maglab, estim_e_rot, estim_eh, find_beta, find_beta2,
    find_beta_maglab, find_beta_rot, find_betah, gen_hist, gen_hist_dynamics, gen_hist_magnus,
    gen_single_traj, get_obs, get_obs_maglab, get_obs_rot, get_obsh, run_dyn_profile, run_ext,
    run_langevin, run_mc, run_mc_magnus, run_mc_profile, run_sim, run_tau, trajectory_mean_e,
};

use clap::{crate_version, App, Arg};

use glob::glob;

fn print_config_description() {
    println!("{}", cs);
}

fn main() {
    //Create spin chain with parameters in file "config.toml"

    let matches = App::new("spinchain")
        .version(crate_version!())
        .author("Thomas Veness <thomas.veness@nottingham.ac.uk>")
        .arg(
            Arg::with_name("tau")
                .short("t")
                .value_name("TAU")
                .takes_value(true)
                .long("tau")
                .help("Override tau of config"),
        )
        .about("Run classical spin chain simulation")
        .arg(
            Arg::with_name("avg")
                .short("a")
                .long("avg")
                .help("Calculate average of runs"),
        )
        .arg(
            Arg::with_name("beta")
            .short("b")
            .long("beta")
            .value_name("BETA")
            .takes_value(true)
            .help("Overrides config beta")
            )
        .arg(
            Arg::with_name("mc")
                .short("m")
                .long("monte-carlo")
                .help("Calculate an average quantitity in Monte-Carlo"),
        )
        .arg(
            Arg::with_name("mc-profile")
                .short("p")
                .long("monte-carlo-profile")
                .help("Calculate an average quantitity in Monte-Carlo, spatially resolved"),
        )
        .arg(
            Arg::with_name("dyn-profile")
                .short("P")
                .long("dynamic-profile")
                .help("Calculate an average quantitity in dynamical runs, spatially resolved"),
        )
        .arg(
            Arg::with_name("mc-full")
                .long("monte-carlo-magnus")
                .help("Calculate an average quantitity in Monte-Carlo, Magnus expansion to second order in omega^{-1}"),
        )
        .arg(
            Arg::with_name("hist")
            .short("h")
            .value_name("POINTS")
            .long("histogram")
            .takes_value(true)
            .default_value("8000")
            .help("Generate histograms via Monte-Carlo (hist_mc.dat) and via time-evolution (hist_dyn.dat)"),
            )
        .arg(
            Arg::with_name("dynhist")
            .short("d")
            .value_name("POINTS")
            .long("dynamic-histogram")
            .takes_value(true)
            .default_value("1000")
            .help("Generate time-evolution histogram (hist_dyn.dat)"),
            )
        .arg(
            Arg::with_name("singletraj")
            .short("s")
            .value_name("POINTS")
            .long("single-trajectory")
            .takes_value(true)
            .default_value("1000")
            .help("Generate time-evolution histogram (hist_sj.dat)"),
            )
        .arg(
            Arg::with_name("langevin")
            .short("l")
            .value_name("GAMMA")
            .long("langevin")
            .takes_value(true)
            .default_value("1.0")
            .help("Directly simulate Langevin dynamics on the system proper"),
            )
        .arg(
            Arg::with_name("magnus-hist")
            .value_name("POINTS")
            .long("magnus-hist")
            .takes_value(true)
            .default_value("8000")
            .help("Generate histogram via Monte-Carlo (hist_mc_magnus.dat) for first-order Magnus expansion, and print averages"),
            )
        .arg(
            Arg::with_name("config-description")
            .long("config-desc")
            .help("Print description of config file"),
            )
        .arg(
            Arg::with_name("n-steps")
            .value_name("TAUS")
            .takes_value(true)
            .long("n-steps")
            .help("Generate single-shot time-evolution at different drive periods tau"),
            )
        .arg(
            Arg::with_name("adiabatic")
            .short("A")
            .value_name("TAU")
            .takes_value(true)
            .long("adiabatic")
            .help("Find adiabatic ensemble H-\\omega S^z"),
            )
        .arg(
            Arg::with_name("adiab2")
            .short("D")
            .value_name("TAU")
            .takes_value(true)
            .long("adiab-full")
            .help("Find adiabatic ensemble H-\\omega S^z inhomogeneous"),
            )
        .arg(
            Arg::with_name("rot").allow_hyphen_values(true)
            .value_name("E")
            .takes_value(true)
            .long("rot-frame")
            .help("Fit temperature of Monte-Carlo ensemble with energy density E in rotating frame")
            )
        .arg(
            Arg::with_name("maglab").allow_hyphen_values(true)
            .value_name("E")
            .takes_value(true)
            .long("magnus-fit")
            .help("Fit temperature of Monte-Carlo ensemble with energy density E in lab frame, leading Magnus")
            )
        .arg(
            Arg::with_name("steps")
            .value_name("STEPS")
            .long("steps")
            .takes_value(true)
            .help("Set step limit when doing n-steps")
            )
        .arg(
            Arg::with_name("fit").allow_hyphen_values(true)
            .value_name("E")
            .short("f")
            .takes_value(true)
            .help("Fit temperature of Monte Carlo ensemble with energy density E")
            )
        .arg(
            Arg::with_name("fit-effective").allow_hyphen_values(true)
            .value_name("E")
            .short("F")
            .takes_value(true)
            .help("Calculate effective ensemble temperature via conserved quantity arguments")
            )
        .arg(
            Arg::with_name("high-freq").allow_hyphen_values(true)
            .value_name("E")
            .short("H")
            .takes_value(true)
            .help("Calculate effective ensemble with initial temp and first-order Magnus")
            )
	.arg(
	    Arg::with_name("ext")
	    .short("e")
            .help("Produce profile of energy density as a function of t"),
)
        .get_matches();
    let mut default = true;
    //Want to read num from file
    let mut conf: Config = SpinChain::read_config("config.toml");
    if let Some(tau) = matches.value_of("tau") {
        conf.tau = tau.parse::<f64>().unwrap();
    }

    match matches.occurrences_of("config-description") {
        0 => {}
        _ => {
            print_config_description();
            default = false;
        }
    }
    match matches.occurrences_of("ext") {
        0 => {}
        _ => {
            println!("Generating ext");
            let mut conf = conf.clone();
            run_ext(&mut conf);
            default = false;
        }
    }

    let steps: u32 = match matches.value_of("steps") {
        Some(x) => x.parse::<u32>().unwrap(),
        None => 1000,
    };

    if let Some(taus) = matches.value_of("n-steps") {
        //        println!("Overriding to: {}", b);
        let tau_vec: Vec<f64> = taus
            .split(',')
            .map(|word| word.parse::<f64>().unwrap())
            .collect();
        println!("{:?}", tau_vec);
        run_tau(tau_vec, steps, &mut conf);
        default = false;
    }

    if let Some(tau) = matches.value_of("adiab2") {
        println!("Adiab 2");
        conf.tau = tau.parse::<f64>().unwrap();

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
    if let Some(targete) = matches.value_of("rot") {
        println!("Evaluating with first-order correction in rot frame");
        //Start with no extra fields paplied to chain
        conf.ednsty = targete.parse::<f64>().unwrap();

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
    if let Some(targete) = matches.value_of("maglab") {
        println!("Evaluating with first-order correction in lab frame");
        //Start with no extra fields paplied to chain
        conf.ednsty = targete.parse::<f64>().unwrap();

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

    if let Some(targete) = matches.value_of("fit") {
        conf.ednsty = targete.parse::<f64>().unwrap();
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

    if let Some(targete) = matches.value_of("fit-effective") {
        conf.ednsty = targete.parse::<f64>().unwrap();
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
    if let Some(targete) = matches.value_of("high-freq") {
        conf.ednsty = targete.parse::<f64>().unwrap();

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

    if let Some(tau) = matches.value_of("adiabatic") {
        //Want to search for beta such that
        // < H_f e^{-\beta H_f} > = e_0
        // Can't think of a clever way to do this,
        // so will set up a grid of beta, and produce estimates for
        // the mean energy. Can then do a bisection search
        conf.tau = tau.parse::<f64>().unwrap();
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

    if let Some(b) = matches.value_of("beta") {
        //        println!("Overriding to: {}", b);
        conf.beta = b.parse::<f64>().unwrap();
    } else {
        conf.beta = SpinChain::solve_beta(conf.ednsty);
    }

    let points: usize = match matches.value_of("hist") {
        Some(x) => x.parse::<usize>().unwrap(),
        None => 8000_usize,
    };
    let sample_num: usize = match matches.value_of("dynhist") {
        Some(x) => x.parse::<usize>().unwrap(),
        None => 1000_usize,
    };
    let single_samples: usize = match matches.value_of("singletraj") {
        Some(x) => x.parse::<usize>().unwrap(),
        None => 1000_usize,
    };

    match matches.occurrences_of("magnus-hist") {
        0 => (),
        _ => {
            println!("Running Monte-Carlo for first-order Magnus expansion");
            default = false;
            gen_hist_magnus(&mut conf, points);
        }
    }

    match matches.occurrences_of("hist") {
        0 => (),
        _ => {
            println!("Running hist");
            default = false;
            gen_hist(&mut conf, points);
        }
    }
    match matches.occurrences_of("dynhist") {
        0 => (),
        _ => {
            println!("Running dynamic hist");
            default = false;
            gen_hist_dynamics(&mut conf, sample_num);
        }
    }

    match matches.occurrences_of("singletraj") {
        0 => (),
        _ => {
            println!("Running single trajectory");
            default = false;
            gen_single_traj(&mut conf, single_samples);
        }
    }

    match matches.is_present("avg") {
        true => {
            default = false;
            average(&mut conf);
        }
        false => (),
    };

    match matches.is_present("mc") {
        true => {
            default = false;
            run_mc(&mut conf);
        }
        false => (),
    };
    match matches.is_present("mc-profile") {
        true => {
            default = false;
            run_mc_profile(&mut conf);
        }
        false => (),
    };

    match matches.is_present("dyn-profile") {
        true => {
            default = false;
            run_dyn_profile(&mut conf);
        }
        false => (),
    };

    match matches.is_present("mc-full") {
        true => {
            default = false;
            run_mc_magnus(&mut conf);
        }
        false => (),
    };

    match matches.occurrences_of("langevin") {
        0 => (),
        _ => {
            let gamma: f64 = matches
                .value_of("langevin")
                .unwrap()
                .parse::<f64>()
                .unwrap();
            println!(
                "Running Langevin dynamics on the system proper with gamma={}",
                gamma
            );
            conf.hsize = conf.ssize;

            run_langevin(&mut conf, gamma);
            default = false;
        }
    }

    match default {
        true => run_sim(&mut conf),
        false => (),
    }
}
