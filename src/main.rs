use indicatif::{MultiProgress, ProgressBar, ProgressStyle};
use std::fs;
use std::fs::{File, OpenOptions};
use std::io::{prelude::*, BufReader};
use threadpool::ThreadPool;
mod config;
use self::config::{Config, DriveType};
use std::sync::mpsc;
use std::thread;

mod macros;

mod spin;
use self::spin::Spin;
#[allow(unused_imports)]
use std::f64::consts::PI;

mod spinchain;
use self::spinchain::SpinChain;

mod spinchain_langevin;
use self::spinchain_langevin::SpinChainLangevin;

use clap::{crate_version, App, Arg};

use glob::glob;

fn write_in_bg(file: File, rx: std::sync::mpsc::Receiver<String>) {
    thread::spawn(move || {
        for received in rx {
            writeln!(&file, "{}", &received).unwrap();
        }
    });
}

fn profile_in_bg(
    file: File,
    rx: std::sync::mpsc::Receiver<String>,
    size: usize,
    pts: usize,
) -> thread::JoinHandle<()> {
    thread::spawn(move || {
        let mut e_vec_mu: Vec<f64> = vec![0.0; size];
        let mut mx_vec_mu: Vec<f64> = vec![0.0; size];
        let mut my_vec_mu: Vec<f64> = vec![0.0; size];
        let mut mz_vec_mu: Vec<f64> = vec![0.0; size];
        // Issues with hanging to join threads: how to fix?
        //        for received in rx {
        for _ in 0..(size * pts) {
            let received = rx.recv().unwrap();
            let v: Vec<&str> = received.split(' ').collect();
            //            println!("{:?}",v);
            let index: usize = v[0].parse::<usize>().unwrap();
            let e: f64 = v[1].parse::<f64>().unwrap();
            let sx: f64 = v[2].parse::<f64>().unwrap();
            let sy: f64 = v[3].parse::<f64>().unwrap();
            let sz: f64 = v[4].parse::<f64>().unwrap();

            e_vec_mu[index] += e;
            mx_vec_mu[index] += sx;
            my_vec_mu[index] += sy;
            mz_vec_mu[index] += sz;
        }
        println!("Now writing data");
        for i in 0..size {
            let mu_e: f64 = e_vec_mu[i] / (pts as f64);
            let mu_x: f64 = mx_vec_mu[i] / (pts as f64);
            let mu_y: f64 = my_vec_mu[i] / (pts as f64);
            let mu_z: f64 = mz_vec_mu[i] / (pts as f64);
            writeln!(&file, "{} {} {} {} {}", i, mu_e, mu_x, mu_y, mu_z).unwrap();
        }
    })
}

fn print_config_description() {
    let c: Config = Config::default();
    println!("# Default config and description of options");
    println!("hsize = {}    # size of the system", c.hsize);
    println!("ssize = {}      # size of subsystem (driven part)", c.ssize);
    println!("t = {}        # final time of simulation", c.t);
    println!("dt = {}      # time-step for simulation", c.dt);
    println!("runs = {}       # number runs to perform", c.runs);
    println!("threads = {}    # number of parallel threads", c.threads);
    println!("trel = {}       # minus initial time", c.trel);
    println!("tau = {}       # period of drive", c.tau);
    println!("lambda = {}     # value of J_z coupling", c.lambda);
    println!(
        "hfield = {:?} # constant h field on entire system",
        c.hfield
    );
    println!("hs = {:?}     # constant h field on subsystem", c.hs);
    println!(
        "jvar = {}   # variance in J couplings (x, y, z independent)",
        c.jvar
    );
    println!("hvar = {}       # variance in field", c.hvar);
    println!(
        "method = {}     # method for numerical integration (2=2nd order Suzuki-Trotter)",
        c.method
    );
    println!("ednsty = {} # energy-density of initial state", c.ednsty);
    println!(
        "file = \"{}\"   # pattern for log files i.e. log0.dat, log1.dat",
        c.file
    );
    println!("strob = {}  # stroboscopic evaluation", c.strob);
    println!("offset = {}     # first file i.e. log0.dat", c.offset);
    println!("beta = {}    # beta (determined from ednsty)", c.beta);
    println!(
        "e = {}    # eccentricity (scale of y-drive if drivetype=\"xyelliptic\")",
        c.e
    );
    println!(
        "drivetype = \"{:?}\" # type of driving, can be \"xyplane\", \"uniaxial\", \"xyelliptic\" ",
        c.drive
    );
}

fn estim_e(conf: &mut Config, beta: f64) -> f64 {
    let mut sc: SpinChain = SpinChain::new(conf.clone(), 0);

    sc.vars.beta = beta;
    let samples: usize = 1000;
    let mut e_samples: Vec<f64> = vec![0.0; samples];
    let omega: f64 = 2.0 * PI / conf.tau;

    for i in 0..samples {
        for _ in 0..1000 {
            sc.metropolis_update();
        }
        e_samples[i] = sc.mc_system_energy() + omega * sc.m()[2] / sc.vars.ssize as f64;
    }

    //    println!("{:?}", e_samples);
    e_samples.iter().fold(0.0, |acc, x| acc + x) / (samples as f64)
}

fn estim_e_maglab(conf: &mut Config, beta: f64) -> f64 {
    let mut sc: SpinChain = SpinChain::new(conf.clone(), 0);

    sc.vars.beta = beta;
    let samples: usize = 1000;
    let mut e_samples: Vec<f64> = vec![0.0; samples];
    //let omega: f64 = 2.0 * PI / conf.tau;

    for i in 0..samples {
        for _ in 0..10000 {
            sc.metropolis_update_maglab();
        }
        e_samples[i] = sc.mc_system_energy_maglab(); // + omega * sc.m()[2] / sc.vars.ssize as f64;
    }

    //    println!("{:?}", e_samples);
    e_samples.iter().fold(0.0, |acc, x| acc + x) / (samples as f64)
}

fn estim_e_rot(conf: &mut Config, beta: f64) -> f64 {
    let mut sc: SpinChain = SpinChain::new(conf.clone(), 0);

    sc.vars.beta = beta;
    let samples: usize = 1000;
    let mut e_samples: Vec<f64> = vec![0.0; samples];
    //let omega: f64 = 2.0 * PI / conf.tau;

    for i in 0..samples {
        for _ in 0..10000 {
            sc.metropolis_update_rot_magnus();
        }
        e_samples[i] = sc.mc_system_energy_rot(); // + omega * sc.m()[2] / sc.vars.ssize as f64;
    }

    //    println!("{:?}", e_samples);
    e_samples.iter().fold(0.0, |acc, x| acc + x) / (samples as f64)
}
fn estim_e2(conf: &mut Config, beta: f64) -> f64 {
    let mut sc: SpinChain = SpinChain::new(conf.clone(), 0);
    //    println!("Field: {:?}", sc.static_h);

    //let omega: f64 = 2.0 * PI / conf.tau;
    sc.vars.beta = beta;
    let samples: usize = 1000;
    let mut e_samples: Vec<f64> = vec![0.0; samples];

    for i in 0..samples {
        for _ in 0..1e3 as usize {
            sc.metropolis_update();
        }
        //        e_samples[i] = sc.mc_total_energy();
        e_samples[i] = sc.mc_total_energy();
    }

    //    println!("{:?}", e_samples);
    e_samples.iter().fold(0.0, |acc, x| acc + x) / (samples as f64)
}
fn estim_eh(conf: &mut Config, beta: f64) -> f64 {
    // let omega: f64 = 2.0 * PI / conf.tau;
    let mut sc: SpinChain = SpinChain::new(conf.clone(), 0);
    //    println!("Field: {:?}", sc.static_h);

    //let omega: f64 = 2.0 * PI / conf.tau;
    sc.vars.beta = beta;
    let samples: usize = 1000;
    let mut e_samples: Vec<f64> = vec![0.0; samples];

    for i in 0..samples {
        for _ in 0..1e3 as usize {
            sc.metropolis_update();
        }
        //        e_samples[i] = sc.mc_total_energy();
        e_samples[i] = sc.mc_total_energy(); // + omega*sc.m()[2]/ sc.vars.ssize as f64;
    }

    //    println!("{:?}", e_samples);
    e_samples.iter().fold(0.0, |acc, x| acc + x) / (samples as f64)
}

fn find_betah(conf: &mut Config) -> f64 {
    let e_target = conf.ednsty;
    let mut beta_low: f64 = 0.01;
    let mut beta_high: f64 = 5.0;
    let mut beta_mid: f64 = 0.5 * (beta_low + beta_high);
    //    println!("Target? {}", estim_e(conf,2.88));

    let mut e_low = estim_eh(conf, beta_low);
    let mut e_mid = estim_eh(conf, beta_mid);
    let mut e_high = estim_eh(conf, beta_high);
    println!("Target: {}", e_target);

    println!("low: {}, mid: {}, high: {}", beta_low, beta_mid, beta_high);
    println!("elow: {}, emid:{}, ehigh: {}", e_low, e_mid, e_high);
    while beta_high - beta_low > 0.001 {
        //If middle energy is too high, then beta_mid is too low (i.e. beta_mid is hot)
        if e_mid > e_target {
            beta_low = beta_mid;
            e_low = e_mid;
        }
        //Otherwise, energy is too low, increase temp i.e. lower beta_high
        else {
            beta_high = beta_mid;
            e_high = e_mid;
        }
        beta_mid = 0.5 * (beta_low + beta_high);
        e_mid = estim_eh(conf, beta_mid);

        println!("low: {}, mid: {}, high: {}", beta_low, beta_mid, beta_high);
        println!("elow: {}, emid:{}, ehigh: {}", e_low, e_mid, e_high);
    }

    beta_mid
}

fn find_beta2(conf: &mut Config) -> f64 {
    let e_target = conf.ednsty;
    let mut beta_low: f64 = 0.01;
    let mut beta_high: f64 = 5.0;
    let mut beta_mid: f64 = 0.5 * (beta_low + beta_high);
    //    println!("Target? {}", estim_e(conf,2.88));

    let mut e_low = estim_e2(conf, beta_low);
    let mut e_mid = estim_e2(conf, beta_mid);
    let mut e_high = estim_e2(conf, beta_high);
    println!("Target: {}", e_target);

    println!("low: {}, mid: {}, high: {}", beta_low, beta_mid, beta_high);
    println!("elow: {}, emid:{}, ehigh: {}", e_low, e_mid, e_high);
    while beta_high - beta_low > 0.001 {
        //If middle energy is too high, then beta_mid is too low (i.e. beta_mid is hot)
        if e_mid > e_target {
            beta_low = beta_mid;
            e_low = e_mid;
        }
        //Otherwise, energy is too low, increase temp i.e. lower beta_high
        else {
            beta_high = beta_mid;
            e_high = e_mid;
        }
        beta_mid = 0.5 * (beta_low + beta_high);
        e_mid = estim_e2(conf, beta_mid);

        println!("low: {}, mid: {}, high: {}", beta_low, beta_mid, beta_high);
        println!("elow: {}, emid:{}, ehigh: {}", e_low, e_mid, e_high);
    }

    beta_mid
}

fn find_beta_maglab(conf: &mut Config) -> f64 {
    let e_target = conf.ednsty;
    let mut beta_low: f64 = 0.0001;
    let mut beta_high: f64 = 5.0;
    let mut beta_mid: f64 = 0.5 * (beta_low + beta_high);
    //    conf.ssize = conf.hsize;
    //    println!("Target? {}", estim_e(conf,2.88));

    let mut e_low = estim_e_maglab(conf, beta_low);
    let mut e_mid = estim_e_maglab(conf, beta_mid);
    let mut e_high = estim_e_maglab(conf, beta_high);
    println!("Target: {}", e_target);

    println!("low: {}, mid: {}, high: {}", beta_low, beta_mid, beta_high);
    println!("elow: {}, emid:{}, ehigh: {}", e_low, e_mid, e_high);
    while beta_high - beta_low > 0.001 {
        //If middle energy is too high, then beta_mid is too low (i.e. beta_mid is hot)
        if e_mid > e_target {
            beta_low = beta_mid;
            e_low = e_mid;
        }
        //Otherwise, energy is too low, increase temp i.e. lower beta_high
        else {
            beta_high = beta_mid;
            e_high = e_mid;
        }
        beta_mid = 0.5 * (beta_low + beta_high);
        e_mid = estim_e_maglab(conf, beta_mid);

        println!("low: {}, mid: {}, high: {}", beta_low, beta_mid, beta_high);
        println!("elow: {}, emid:{}, ehigh: {}", e_low, e_mid, e_high);
    }

    beta_mid
}

fn find_beta_rot(conf: &mut Config) -> f64 {
    let e_target = conf.ednsty;
    let mut beta_low: f64 = 0.0001;
    let mut beta_high: f64 = 5.0;
    let mut beta_mid: f64 = 0.5 * (beta_low + beta_high);
    //    conf.ssize = conf.hsize;
    //    println!("Target? {}", estim_e(conf,2.88));

    let mut e_low = estim_e_rot(conf, beta_low);
    let mut e_mid = estim_e_rot(conf, beta_mid);
    let mut e_high = estim_e_rot(conf, beta_high);
    println!("Target: {}", e_target);

    println!("low: {}, mid: {}, high: {}", beta_low, beta_mid, beta_high);
    println!("elow: {}, emid:{}, ehigh: {}", e_low, e_mid, e_high);
    while beta_high - beta_low > 0.001 {
        //If middle energy is too high, then beta_mid is too low (i.e. beta_mid is hot)
        if e_mid > e_target {
            beta_low = beta_mid;
            e_low = e_mid;
        }
        //Otherwise, energy is too low, increase temp i.e. lower beta_high
        else {
            beta_high = beta_mid;
            e_high = e_mid;
        }
        beta_mid = 0.5 * (beta_low + beta_high);
        e_mid = estim_e_rot(conf, beta_mid);

        println!("low: {}, mid: {}, high: {}", beta_low, beta_mid, beta_high);
        println!("elow: {}, emid:{}, ehigh: {}", e_low, e_mid, e_high);
    }

    beta_mid
}

fn find_beta(conf: &mut Config) -> f64 {
    let e_target = conf.ednsty;
    let mut beta_low: f64 = 0.01;
    let mut beta_high: f64 = 5.0;
    let mut beta_mid: f64 = 0.5 * (beta_low + beta_high);
    conf.ssize = conf.hsize;
    //    println!("Target? {}", estim_e(conf,2.88));

    let mut e_low = estim_e(conf, beta_low);
    let mut e_mid = estim_e(conf, beta_mid);
    let mut e_high = estim_e(conf, beta_high);
    println!("Target: {}", e_target);

    println!("low: {}, mid: {}, high: {}", beta_low, beta_mid, beta_high);
    println!("elow: {}, emid:{}, ehigh: {}", e_low, e_mid, e_high);
    while beta_high - beta_low > 0.001 {
        //If middle energy is too high, then beta_mid is too low (i.e. beta_mid is hot)
        if e_mid > e_target {
            beta_low = beta_mid;
            e_low = e_mid;
        }
        //Otherwise, energy is too low, increase temp i.e. lower beta_high
        else {
            beta_high = beta_mid;
            e_high = e_mid;
        }
        beta_mid = 0.5 * (beta_low + beta_high);
        e_mid = estim_e(conf, beta_mid);

        println!("low: {}, mid: {}, high: {}", beta_low, beta_mid, beta_high);
        println!("elow: {}, emid:{}, ehigh: {}", e_low, e_mid, e_high);
    }

    beta_mid
}

fn get_obs_maglab(conf: &mut Config) -> [f64; 4] {
    //    println!("Target? {}", estim_e(conf,2.88));
    //
    conf.ssize = conf.hsize;

    let mut sc: SpinChain = SpinChain::new(conf.clone(), 0);

    let samples: usize = 8000;
    let mut e_samples: Vec<f64> = vec![0.0; samples];
    let mut mx_samples: Vec<f64> = vec![0.0; samples];
    let mut my_samples: Vec<f64> = vec![0.0; samples];
    let mut mz_samples: Vec<f64> = vec![0.0; samples];

    for i in 0..samples {
        for _ in 0..10000 {
            sc.metropolis_update_maglab();
        }

        e_samples[i] = sc.mc_system_energy() + sc.m()[0] / sc.vars.ssize as f64;
        mx_samples[i] = sc.m()[0] / sc.vars.ssize as f64;
        my_samples[i] = sc.m()[1] / sc.vars.ssize as f64;
        mz_samples[i] = sc.m()[2] / sc.vars.ssize as f64;
    }

    //    println!("{:?}", e_samples);
    [
        e_samples.iter().fold(0.0, |acc, x| acc + x) / (samples as f64),
        mx_samples.iter().fold(0.0, |acc, x| acc + x) / (samples as f64),
        my_samples.iter().fold(0.0, |acc, x| acc + x) / (samples as f64),
        mz_samples.iter().fold(0.0, |acc, x| acc + x) / (samples as f64),
    ]
}

///Get observations by MC-sampling from Magnus enemble, but evaluating regular observables
fn get_obs_rot(conf: &mut Config) -> [f64; 4] {
    //    println!("Target? {}", estim_e(conf,2.88));
    //
    conf.ssize = conf.hsize;

    let mut sc: SpinChain = SpinChain::new(conf.clone(), 0);

    let samples: usize = 8000;
    let mut e_samples: Vec<f64> = vec![0.0; samples];
    let mut mx_samples: Vec<f64> = vec![0.0; samples];
    let mut my_samples: Vec<f64> = vec![0.0; samples];
    let mut mz_samples: Vec<f64> = vec![0.0; samples];

    for i in 0..samples {
        for _ in 0..10000 {
            sc.metropolis_update_rot_magnus();
        }

        //e_samples[i] = sc.mc_system_energy_rot();
        e_samples[i] = sc.mc_system_energy() + sc.m()[0] / sc.vars.ssize as f64;
        mx_samples[i] = sc.m()[0] / sc.vars.ssize as f64;
        my_samples[i] = sc.m()[1] / sc.vars.ssize as f64;
        mz_samples[i] = sc.m()[2] / sc.vars.ssize as f64;
    }

    //    println!("{:?}", e_samples);
    [
        e_samples.iter().fold(0.0, |acc, x| acc + x) / (samples as f64),
        mx_samples.iter().fold(0.0, |acc, x| acc + x) / (samples as f64),
        my_samples.iter().fold(0.0, |acc, x| acc + x) / (samples as f64),
        mz_samples.iter().fold(0.0, |acc, x| acc + x) / (samples as f64),
    ]
}

fn get_obsh(conf: &mut Config) -> [f64; 4] {
    conf.ssize = conf.hsize;
    //    println!("Target? {}", estim_e(conf,2.88));

    let mut sc: SpinChain = SpinChain::new(conf.clone(), 0);

    let samples: usize = 8000;
    let mut e_samples: Vec<f64> = vec![0.0; samples];
    let mut mx_samples: Vec<f64> = vec![0.0; samples];
    let mut my_samples: Vec<f64> = vec![0.0; samples];
    let mut mz_samples: Vec<f64> = vec![0.0; samples];

    for i in 0..samples {
        for _ in 0..1000 {
            sc.metropolis_update();
        }
        e_samples[i] = sc.mc_system_energy();
        mx_samples[i] = sc.m()[0] / sc.vars.ssize as f64;
        my_samples[i] = sc.m()[1] / sc.vars.ssize as f64;
        mz_samples[i] = sc.m()[2] / sc.vars.ssize as f64;
    }

    //    println!("{:?}", e_samples);
    [
        e_samples.iter().fold(0.0, |acc, x| acc + x) / (samples as f64),
        mx_samples.iter().fold(0.0, |acc, x| acc + x) / (samples as f64),
        my_samples.iter().fold(0.0, |acc, x| acc + x) / (samples as f64),
        mz_samples.iter().fold(0.0, |acc, x| acc + x) / (samples as f64),
    ]
}

fn get_obs(conf: &mut Config) -> [f64; 4] {
    conf.ssize = conf.hsize;
    //    println!("Target? {}", estim_e(conf,2.88));

    let mut sc: SpinChain = SpinChain::new(conf.clone(), 0);

    let samples: usize = 8000;
    let mut e_samples: Vec<f64> = vec![0.0; samples];
    let mut mx_samples: Vec<f64> = vec![0.0; samples];
    let mut my_samples: Vec<f64> = vec![0.0; samples];
    let mut mz_samples: Vec<f64> = vec![0.0; samples];
    let omega: f64 = 2.0 * PI / conf.tau;

    for i in 0..samples {
        for _ in 0..1000 {
            sc.metropolis_update();
        }
        e_samples[i] = sc.mc_system_energy() + omega * sc.m()[2] / sc.vars.ssize as f64;
        mx_samples[i] = sc.m()[0] / sc.vars.ssize as f64;
        my_samples[i] = sc.m()[1] / sc.vars.ssize as f64;
        mz_samples[i] = sc.m()[2] / sc.vars.ssize as f64;
    }

    //    println!("{:?}", e_samples);
    [
        e_samples.iter().fold(0.0, |acc, x| acc + x) / (samples as f64),
        mx_samples.iter().fold(0.0, |acc, x| acc + x) / (samples as f64),
        my_samples.iter().fold(0.0, |acc, x| acc + x) / (samples as f64),
        mz_samples.iter().fold(0.0, |acc, x| acc + x) / (samples as f64),
    ]
}

///Generate ```sample_num``` samples via Monte Carlo using first-order Magnus expansion system
///indeed thermalises via time-evolution
fn gen_hist_magnus(conf: &mut Config, sample_num: usize) {
    conf.drive = DriveType::none;
    conf.file = "x".to_string();
    let mut sc: SpinChain = SpinChain::new(conf.clone(), 0);
    let s: f64 = sc.vars.ssize as f64;

    let mc_hist = OpenOptions::new()
        .create(true)
        .append(true)
        .open("hist_mc_magnus.dat")
        .unwrap();

    println!("Generating Monte-Carlo histogram");

    let pb = ProgressBar::new(sample_num as u64);
    pb.set_style(
        ProgressStyle::default_bar()
            .template(
                "{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} ({eta})",
            )
            .progress_chars("#>-"),
    );
    let mut mx_avg: f64 = 0.0;
    let mut my_avg: f64 = 0.0;
    let mut mz_avg: f64 = 0.0;
    let mut ed_avg: f64 = 0.0;

    // First generate Monte-Carlo histogram
    // Generate 1000 points, then bin them
    for i in 0..sample_num {
        //Do Metropolis updates
        for _ in 0..1e4 as usize {
            sc.metropolis_update_magnus();
        }

        //Get a sample
        let mx: f64 = sc.m()[0] / sc.vars.ssize as f64;
        let my: f64 = sc.m()[1] / sc.vars.ssize as f64;
        let mz: f64 = sc.m()[2] / sc.vars.ssize as f64;
        let ed: f64 = sc.system_energy();
        let sm: [f64; 3] = sc.spins[(s / 4.0) as usize].dir;
        let sp: [f64; 3] = sc.spins[(3.0 * s / 4.0) as usize].dir;

        mx_avg += mx;
        my_avg += my;
        mz_avg += mz;
        ed_avg += ed;

        writeln!(
            &mc_hist,
            "{} {} {} {} {} {} {} {} {} {} {}",
            i, mx, my, mz, ed, sm[0], sm[1], sm[2], sp[0], sp[1], sp[2]
        )
        .unwrap();

        pb.inc(1);
    }

    pb.finish_with_message("Done");
    println!("mx: {}", mx_avg / sample_num as f64);
    println!("my: {}", my_avg / sample_num as f64);
    println!("mz: {}", mz_avg / sample_num as f64);
    println!("ed: {}", ed_avg / sample_num as f64);
}

fn trajectory_mean_e(conf: &mut Config) -> f64 {
    let (tx, rx) = mpsc::channel();

    let m = MultiProgress::new();
    let sty = ProgressStyle::default_bar()
        .template(
            "{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] t={pos}/{len} ({eta}) {msg}",
        )
        .progress_chars("#>-");

    let pool = ThreadPool::new(conf.threads);

    println!("Generating single trajectory energy estimate");

    for i in 0..((conf.runs) as usize) {
        let txc = mpsc::Sender::clone(&tx);
        let pb = m.add(ProgressBar::new((1000.0 * conf.tau) as u64));
        let steps: usize = (1000.0 * conf.tau / conf.dt) as usize;
        let steps_in_cycle: usize = (conf.tau / conf.dt) as usize;
        pb.set_style(sty.clone());

        let mut sc: SpinChain = SpinChain::new(conf.clone(), 0);

        pool.execute(move || {
            pb.set_message(&format!("run {i}", i = i));
            //Initialise chain
            for _ in 0..2e7 as usize {
                sc.metropolis_update();
            }
            //Do dynamical updates
            pb.reset_eta();
            for cntr in 0..steps {
                sc.update();
                if cntr % steps_in_cycle == 0 {
                    pb.set_position(sc.t as u64);
                    if cntr > 99 * steps_in_cycle {
                        let ed: f64 = sc.system_energy();
                        txc.send(ed).unwrap();
                    }
                }
            }

            pb.finish_and_clear();
        });
    }

    let mut running_total: f64 = 0.0;
    let mut running_num: f64 = 0.0;

    m.join().unwrap();

    //Close channel so we can average
    drop(tx);

    for received in rx {
        running_num += 1.0;
        running_total += received;
        //        println!("e: {}", received);
    }
    running_total / running_num

    // Now process txc
}

fn gen_single_traj(conf: &mut Config, sample_num: usize) {
    let hist_name: String = format!("hist-st-t-{}.dat", conf.t);

    let dyn_hist = OpenOptions::new()
        .create(true)
        .append(true)
        .open(&hist_name)
        .unwrap();

    //To do samples, run for 30 cycles to equilibrate
    //Then sample every 10 cycles?
    let offset: f64 = 30.0;
    let final_t: f64 = (sample_num as f64 + offset) * conf.tau;
    let mut sc: SpinChain = SpinChain::new(conf.clone(), 0);

    let pb = ProgressBar::new(final_t as u64);
    let sty = ProgressStyle::default_bar()
        .template(
            "{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] t={pos}/{len} ({eta}) {msg}",
        )
        .progress_chars("#>-");
    pb.set_style(sty);

    //First do some MC updates in order to get on the correct trajectory
    //for _ in 0..100 {
    //    sc.metropolis_update();
    //}

    let init_steps = (offset * conf.tau / conf.dt) as usize;
    let steps = (final_t / conf.dt) as usize;
    println!(
        "Generating dynamics histogram for {} steps i.e. t={}",
        steps,
        (steps as f64) * conf.dt
    );
    let mut step_count: u64 = 0;
    let cycle_steps: u64 = (conf.tau / conf.dt) as u64;

    println!("Sample num: {}", sample_num);
    for _ in 0..init_steps {
        sc.update();
        pb.set_position(sc.t as u64);
        step_count += 1;
    }
    while sc.t < final_t {
        sc.update();
        step_count += 1;

        pb.set_position(sc.t as u64);
        if step_count % cycle_steps == 0 {
            //Get a sample
            let mx: f64 = sc.m()[0] / sc.vars.ssize as f64;
            let my: f64 = sc.m()[1] / sc.vars.ssize as f64;
            let mz: f64 = sc.m()[2] / sc.vars.ssize as f64;
            let ed: f64 = sc.system_energy();

            let formatted_string: String = format!("{} {} {} {} {}", sc.t, mx, my, mz, ed);

            writeln!(&dyn_hist, "{}", formatted_string).unwrap();
        }
    }
}

///Runs single-shot trajectories and logs observables at the end of the trajectory
fn gen_hist_dynamics(conf: &mut Config, sample_num: usize) {
    let hist_name: String = format!("dyn-hist-t-{}.dat", conf.t);
    let dyn_hist = OpenOptions::new()
        .create(true)
        .append(true)
        .open(&hist_name)
        .unwrap();

    let (tx, rx) = mpsc::channel();
    write_in_bg(dyn_hist, rx);
    let m = MultiProgress::new();
    let sty = ProgressStyle::default_bar()
        .template(
            "{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] t={pos}/{len} ({eta}) {msg}",
        )
        .progress_chars("#>-");

    let steps = (conf.t / conf.dt) as usize;
    println!(
        "Generating dynamics histogram for {} steps i.e. t={}",
        steps,
        (steps as f64) * conf.dt
    );

    let threads = match conf.threads {
        _ if conf.threads >= sample_num => sample_num,
        _ => conf.threads,
    };

    let pool = ThreadPool::new(threads);

    let f = File::create("log0.dat").unwrap();

    println!("Sample num: {}", sample_num);
    for i in 0..sample_num {
        let txc = mpsc::Sender::clone(&tx);
        let pb = m.add(ProgressBar::new(conf.t as u64));
        pb.set_style(sty.clone());

        let mut sc: SpinChain = SpinChain::new(conf.clone(), 0);
        //        sc.file = f.try_clone().unwrap();

        pool.execute(move || {
            pb.set_message(&format!("run {i}", i = i));
            //Do dynamical updates
            pb.reset_eta();
            for _ in 0..steps {
                sc.update();
                if steps % 100 == 0 {
                    pb.set_position(sc.t as u64);
                }
            }

            //Get a sample
            let mx: f64 = sc.m()[0] / sc.vars.ssize as f64;
            let my: f64 = sc.m()[1] / sc.vars.ssize as f64;
            let mz: f64 = sc.m()[2] / sc.vars.ssize as f64;
            let ed: f64 = sc.system_energy();

            let formatted_string: String = format!("{} {} {} {} {}", i, mx, my, mz, ed);
            txc.send(formatted_string).unwrap();

            pb.finish_and_clear();
        });
    }
    m.join().unwrap();
}

///Generate ```sample_num``` samples via Monte Carlo and dynamical evolution to check that the
///system indeed thermalises via time-evolution
fn gen_hist(conf: &mut Config, sample_num: usize) {
    conf.drive = DriveType::none;
    conf.file = "x".to_string();

    let mut sc: SpinChain = SpinChain::new(conf.clone(), 0);

    let mc_hist = OpenOptions::new()
        .create(true)
        .append(true)
        .open("hist_mc.dat")
        .unwrap();

    let dyn_hist = OpenOptions::new()
        .create(true)
        .append(true)
        .open("hist_dyn.dat")
        .unwrap();
    //    let mc_hist = File::create("hist_mc.dat").expect("Could not open file hist_mc.dat");
    //    let dyn_hist = File::create("hist_dyn.dat").expect("Could not open file hist_mc.dat");

    println!("Generating Monte-Carlo histogram");
    let pb = ProgressBar::new(sample_num as u64);
    pb.set_style(
        ProgressStyle::default_bar()
            .template(
                "{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} ({eta})",
            )
            .progress_chars("#>-"),
    );

    // First generate Monte-Carlo histogram
    // Generate 1000 points, then bin them
    for i in 0..sample_num {
        //Do Metropolis updates
        for _ in 0..1e4 as usize {
            sc.metropolis_update();
        }

        //Get a sample
        let mx: f64 = sc.m()[0] / sc.vars.ssize as f64;
        let my: f64 = sc.m()[1] / sc.vars.ssize as f64;
        let mz: f64 = sc.m()[2] / sc.vars.ssize as f64;
        let ed: f64 = sc.system_energy();
        let s: f64 = sc.vars.ssize as f64;
        let sm: [f64; 3] = sc.spins[(s / 4.0) as usize].dir;
        let sp: [f64; 3] = sc.spins[(3.0 * s / 4.0) as usize].dir;

        writeln!(
            &mc_hist,
            "{} {} {} {} {} {} {} {} {} {} {}",
            i, mx, my, mz, ed, sm[0], sm[1], sm[2], sp[0], sp[1], sp[2]
        )
        .unwrap();

        pb.inc(1);
    }
    pb.finish_with_message("Done");
    let pb = ProgressBar::new(sample_num as u64);
    pb.set_style(
        ProgressStyle::default_bar()
            .template(
                "{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} ({eta})",
            )
            .progress_chars("#>-"),
    );
    println!("Generating dynamics histogram");

    for i in 0..sample_num {
        //Do dynamical updates
        for _ in 0..1e3 as usize {
            sc.update();
        }

        //Get a sample
        let mx: f64 = sc.m()[0] / sc.vars.ssize as f64;
        let my: f64 = sc.m()[1] / sc.vars.ssize as f64;
        let mz: f64 = sc.m()[2] / sc.vars.ssize as f64;
        //        let ed: f64 = sc.system_energy();
        //        Make it comparable to MC

        let ed: f64 = sc.system_energy();

        writeln!(&dyn_hist, "{} {} {} {} {}", i, mx, my, mz, ed).unwrap();

        pb.inc(1);
    }
    pb.finish_with_message("Done");
}

fn run_tau(taus: Vec<f64>, steps: u32, conf: &mut Config) {
    let m = MultiProgress::new();
    let sty = ProgressStyle::default_bar()
        .template(
            "{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] t={pos}/{len} ({eta}) {msg}",
        )
        .progress_chars("#>-");

    let threads = conf.threads;
    let pool = ThreadPool::new(threads);

    let (tx, rx) = mpsc::channel();

    //Write logs
    let evt_string: &str = &["evt".to_string(), (steps).to_string(), ".dat".to_string()].join("");

    let evt_file = OpenOptions::new()
        .create(true)
        .append(true)
        .open(&evt_string)
        .unwrap();

    write_in_bg(evt_file, rx);

    println!("Energy density: {}", conf.ednsty);
    println!("Effective temperature: {}", conf.beta);

    for (i, tau) in taus.iter().enumerate() {
        let mut spin_chain: SpinChain = SpinChain::new(conf.clone(), i + conf.offset as usize);
        spin_chain.vars.tau = *tau;
        spin_chain.vars.dt = match spin_chain.vars.tau {
            x if x < 10.0 => spin_chain.vars.tau / 100.0,
            _ => 0.1,
        };
        spin_chain.vars.t = steps as f64 * spin_chain.vars.tau;

        // Initialise at a particular temperature, say T=1

        let pb = m.add(ProgressBar::new(spin_chain.vars.t as u64));
        pb.set_style(sty.clone());
        let txc = mpsc::Sender::clone(&tx);

        pool.execute(move || {
            pb.set_message(&format!("tau={tau}", tau = spin_chain.vars.tau));

            for _ in 0..2e7 as usize {
                spin_chain.metropolis_update();
            }
            pb.reset_eta();

            spin_chain.log();

            while spin_chain.t < spin_chain.vars.t {
                spin_chain.update();

                if spin_chain.t.fract() < spin_chain.vars.dt {
                    pb.set_position(spin_chain.t as u64);
                }
            }

            let es = spin_chain.system_energy();
            let e = spin_chain.total_energy2();
            let m: [f64; 3] = spin_chain.m();
            let s: f64 = spin_chain.vars.ssize as f64;
            let formatted_string: String = format!(
                "{tau} {t} {e_total} {e_sub} {mx} {my} {mz}",
                tau = spin_chain.vars.tau,
                t = spin_chain.t,
                e_total = e,
                e_sub = es,
                mx = m[0] / s,
                my = m[1] / s,
                mz = m[2] / s
            );
            txc.send(formatted_string).unwrap();

            pb.finish_and_clear();
            //            pb.finish_with_message(&format!("Done {}",i));
        });
    }

    m.join().unwrap();
}

fn run_langevin(conf: &mut Config, gamma: f64) {
    let m = MultiProgress::new();
    let sty = ProgressStyle::default_bar()
        .template(
            "{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] t={pos}/{len} ({eta}) {msg}",
        )
        .progress_chars("#>-");

    let threads = conf.threads;
    let num: usize = conf.runs as usize;
    let pool = ThreadPool::new(threads);

    println!("Energy density: {}", conf.ednsty);
    println!("Effective temperature: {}", conf.beta);

    for i in 0..num {
        let mut spin_chain = SpinChainLangevin::new(conf.clone(), i + conf.offset as usize, gamma);

        // Initialise at a particular temperature, say T=1

        let pb = m.add(ProgressBar::new(spin_chain.vars.t as u64));
        pb.set_style(sty.clone());

        pool.execute(move || {
            pb.set_message(&format!("Run {}", i));

            for _ in 0..2e7 as usize {
                spin_chain.metropolis_update();
            }
            pb.reset_eta();
            //Round is necessary!
            let tau_steps: u64 = (spin_chain.vars.tau / spin_chain.vars.dt).round() as u64;

            spin_chain.log();

            while spin_chain.t < spin_chain.vars.t {
                spin_chain.update();

                if spin_chain.vars.strob {
                    if ((spin_chain.t / spin_chain.vars.dt) as u64).rem_euclid(tau_steps) == 0 {
                        spin_chain.log();
                    }
                } else {
                    spin_chain.log();
                }

                if spin_chain.t.fract() < spin_chain.vars.dt {
                    pb.set_position(spin_chain.t as u64);
                }
            }
            pb.finish_and_clear();
            //            pb.finish_with_message(&format!("Done {}",i));
        });
    }

    //    m.join_and_clear().unwrap();
    m.join().unwrap();

    //Update config
    conf.offset += num as u32;
    fs::write("config.toml", toml::to_string(&conf).unwrap()).unwrap();
    println!("Finished {} runs", num);
}

fn run_sim(conf: &mut Config) {
    let m = MultiProgress::new();
    let sty = ProgressStyle::default_bar()
        .template(
            "{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] t={pos}/{len} ({eta}) {msg}",
        )
        .progress_chars("#>-");

    let threads = conf.threads;
    let num: usize = conf.runs as usize;
    let pool = ThreadPool::new(threads);

    println!("Energy density: {}", conf.ednsty);
    println!("Effective temperature: {}", conf.beta);

    for i in 0..num {
        let mut spin_chain: SpinChain = SpinChain::new(conf.clone(), i + conf.offset as usize);

        // Initialise at a particular temperature, say T=1

        let pb = m.add(ProgressBar::new(spin_chain.vars.t as u64));
        pb.set_style(sty.clone());

        pool.execute(move || {
            pb.set_message(&format!("Run {}", i));

            //            spin_chain.vars.hs=vec![1.0,0.0,-2.0*PI/spin_chain.vars.tau];

            for _ in 0..2e7 as usize {
                spin_chain.metropolis_update();
            }
            pb.reset_eta();
            let tau_steps: u64 = (spin_chain.vars.tau / spin_chain.vars.dt) as u64;

            spin_chain.log();

            while spin_chain.t < spin_chain.vars.t {
                spin_chain.update();

                if spin_chain.vars.strob {
                    if ((spin_chain.t / spin_chain.vars.dt) as u64).rem_euclid(tau_steps) == 0 {
                        spin_chain.log();
                    }
                } else {
                    spin_chain.log();
                }

                if spin_chain.t.fract() < spin_chain.vars.dt {
                    pb.set_position(spin_chain.t as u64);
                }
            }
            pb.finish_and_clear();
            //            pb.finish_with_message(&format!("Done {}",i));
        });
    }

    //    m.join_and_clear().unwrap();
    m.join().unwrap();

    //Update config
    conf.offset += num as u32;
    fs::write("config.toml", toml::to_string(&conf).unwrap()).unwrap();
    println!("Finished {} runs", num);
}

fn run_dyn_profile(conf: &mut Config) {
    let m = MultiProgress::new();
    let sty = ProgressStyle::default_bar()
        .template(
            "{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] t={pos}/{len} ({eta}) {msg}",
        )
        .progress_chars("#>-");

    let threads = conf.threads;
    let num: usize = conf.runs as usize;
    let pool = ThreadPool::new(threads);
    // Set up to run for 4000 cycles
    // Log for every 500 cycles
    conf.t = conf.tau * 4000.0;

    //

    let file_500 = File::create("dyn_profile-500-tau.dat").unwrap();
    let file_1000 = File::create("dyn_profile-1000-tau.dat").unwrap();
    let file_1500 = File::create("dyn_profile-1500-tau.dat").unwrap();
    let file_2000 = File::create("dyn_profile-2000-tau.dat").unwrap();
    let file_2500 = File::create("dyn_profile-2500-tau.dat").unwrap();
    let file_3000 = File::create("dyn_profile-3000-tau.dat").unwrap();
    let file_3500 = File::create("dyn_profile-3500-tau.dat").unwrap();
    let file_4000 = File::create("dyn_profile-4000-tau.dat").unwrap();
    let (tx_500, rx_500) = mpsc::channel();
    let (tx_1000, rx_1000) = mpsc::channel();
    let (tx_1500, rx_1500) = mpsc::channel();
    let (tx_2000, rx_2000) = mpsc::channel();
    let (tx_2500, rx_2500) = mpsc::channel();
    let (tx_3000, rx_3000) = mpsc::channel();
    let (tx_3500, rx_3500) = mpsc::channel();
    let (tx_4000, rx_4000) = mpsc::channel();
    let writing_thread_500 = profile_in_bg(file_500, rx_500, conf.hsize as usize, num);
    let writing_thread_1000 = profile_in_bg(file_1000, rx_1000, conf.hsize as usize, num);
    let writing_thread_1500 = profile_in_bg(file_1500, rx_1500, conf.hsize as usize, num);
    let writing_thread_2000 = profile_in_bg(file_2000, rx_2000, conf.hsize as usize, num);
    let writing_thread_2500 = profile_in_bg(file_2500, rx_2500, conf.hsize as usize, num);
    let writing_thread_3000 = profile_in_bg(file_3000, rx_3000, conf.hsize as usize, num);
    let writing_thread_3500 = profile_in_bg(file_3500, rx_3500, conf.hsize as usize, num);
    let writing_thread_4000 = profile_in_bg(file_4000, rx_4000, conf.hsize as usize, num);

    println!("Energy density: {}", conf.ednsty);
    println!("Effective temperature: {}", conf.beta);
    //    let mut sc_orig: SpinChain = SpinChain::new(conf.clone(), conf.offset as usize).clone();

    for i in 0..num {
        let mut spin_chain: SpinChain = SpinChain::new(conf.clone(), i + conf.offset as usize);
        //        let mut spin_chain: SpinChain = sc_orig.clone();
        //        let txc = mpsc::Sender::clone(&tx);
        let txc_500 = mpsc::Sender::clone(&tx_500);
        let txc_1000 = mpsc::Sender::clone(&tx_1000);
        let txc_1500 = mpsc::Sender::clone(&tx_1500);
        let txc_2000 = mpsc::Sender::clone(&tx_2000);
        let txc_2500 = mpsc::Sender::clone(&tx_2500);
        let txc_3000 = mpsc::Sender::clone(&tx_3000);
        let txc_3500 = mpsc::Sender::clone(&tx_3500);
        let txc_4000 = mpsc::Sender::clone(&tx_4000);

        // Initialise at a particular temperature, say T=1

        let pb = m.add(ProgressBar::new(spin_chain.vars.t as u64));
        pb.set_style(sty.clone());

        pool.execute(move || {
            pb.set_message(&format!("Run {}", i));

            //            spin_chain.vars.hs=vec![1.0,0.0,-2.0*PI/spin_chain.vars.tau];

            for _ in 0..2e7 as usize {
                spin_chain.metropolis_update();
            }
            pb.reset_eta();
            let tau_steps: u64 = (spin_chain.vars.tau / spin_chain.vars.dt) as u64;

            while spin_chain.t < spin_chain.vars.t {
                spin_chain.update();
                if spin_chain.t.fract() < spin_chain.vars.dt {
                    pb.set_position(spin_chain.t as u64);
                }

                let txc_matched = match spin_chain.t {
                    _ if (spin_chain.t - 500.0 * spin_chain.vars.tau).abs()
                        < spin_chain.vars.dt / 5.0 =>
                    {
                        Some(mpsc::Sender::clone(&txc_500))
                    }
                    _ if (spin_chain.t - 1000.0 * spin_chain.vars.tau).abs()
                        < spin_chain.vars.dt / 5.0 =>
                    {
                        Some(mpsc::Sender::clone(&txc_1000))
                    }
                    _ if (spin_chain.t - 1500.0 * spin_chain.vars.tau).abs()
                        < spin_chain.vars.dt / 5.0 =>
                    {
                        Some(mpsc::Sender::clone(&txc_1500))
                    }
                    _ if (spin_chain.t - 2000.0 * spin_chain.vars.tau).abs()
                        < spin_chain.vars.dt / 5.0 =>
                    {
                        Some(mpsc::Sender::clone(&txc_2000))
                    }
                    _ if (spin_chain.t - 2500.0 * spin_chain.vars.tau).abs()
                        < spin_chain.vars.dt / 5.0 =>
                    {
                        Some(mpsc::Sender::clone(&txc_2500))
                    }
                    _ if (spin_chain.t - 3000.0 * spin_chain.vars.tau).abs()
                        < spin_chain.vars.dt / 5.0 =>
                    {
                        Some(mpsc::Sender::clone(&txc_3000))
                    }
                    _ if (spin_chain.t - 3500.0 * spin_chain.vars.tau).abs()
                        < spin_chain.vars.dt / 5.0 =>
                    {
                        Some(mpsc::Sender::clone(&txc_3500))
                    }
                    _ if (spin_chain.t - 4000.0 * spin_chain.vars.tau).abs()
                        < spin_chain.vars.dt / 5.0 =>
                    {
                        Some(mpsc::Sender::clone(&txc_4000))
                    }
                    _ => None,
                };

                if let Some(txc) = txc_matched {
                    //Now log with correct log
                    //Write data
                    let L: usize = spin_chain.vars.hsize as usize;
                    for i in 0..L {
                        let s_l: &Spin = match i {
                            _ if i == 0 => &spin_chain.spins[L - 1],
                            _ => &spin_chain.spins[i - 1],
                        };
                        let s_c: &Spin = &spin_chain.spins[i];
                        let s_r: &Spin = match i {
                            _ if i == L - 1 => &spin_chain.spins[0],
                            _ => &spin_chain.spins[i + 1],
                        };
                        let j_l: &[f64; 3] = match i {
                            _ if i == 0 => &spin_chain.j_couple[L - 1],
                            _ => &spin_chain.j_couple[i - 1],
                        };

                        let j_r: &[f64; 3] = &spin_chain.j_couple[i];

                        let mut sjs: f64 = 0.0;
                        for k in 0..3 {
                            sjs -= s_l.dir[k] * j_l[k] * s_c.dir[k];
                            sjs -= s_c.dir[k] * j_r[k] * s_r.dir[k];
                        }
                        sjs = sjs / 2.0;

                        let formatted_string: String =
                            format!("{} {} {} {} {}", i, sjs, s_c.dir[0], s_c.dir[1], s_c.dir[2]);
                        txc.send(formatted_string).unwrap();
                    }
                }
            }

            pb.finish_and_clear();
            //            pb.finish_with_message(&format!("Done {}",i));
        });
    }

    //    m.join_and_clear().unwrap();
    m.join().unwrap();
    writing_thread_500.join().unwrap();
    writing_thread_1000.join().unwrap();
    writing_thread_1500.join().unwrap();
    writing_thread_2000.join().unwrap();
    writing_thread_2500.join().unwrap();
    writing_thread_3000.join().unwrap();
    writing_thread_3500.join().unwrap();
    writing_thread_4000.join().unwrap();

    //Update config
    conf.offset += num as u32;
    fs::write("config.toml", toml::to_string(&conf).unwrap()).unwrap();
    println!("Finished {} runs", num);
}

fn run_ext(conf: &mut Config) {
    let m = MultiProgress::new();
    let sty = ProgressStyle::default_bar()
        .template(
            "{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] t={pos}/{len} ({eta}) {msg}",
        )
        .progress_chars("#>-");

    let threads = conf.threads;
    let num: usize = conf.runs as usize;
    let pool = ThreadPool::new(threads);

    println!("Energy density: {}", conf.ednsty);
    println!("Effective temperature: {}", conf.beta);

    for i in 0..num {
        let mut spin_chain: SpinChain = SpinChain::new(conf.clone(), i + conf.offset as usize);

        // Initialise at a particular temperature, say T=1

        let pb = m.add(ProgressBar::new(spin_chain.vars.t as u64));
        pb.set_style(sty.clone());

        pool.execute(move || {
            pb.set_message(&format!("Run {}", i));

            //            spin_chain.vars.hs=vec![1.0,0.0,-2.0*PI/spin_chain.vars.tau];

            for _ in 0..2e7 as usize {
                spin_chain.metropolis_update();
            }
            pb.reset_eta();
            let tau_steps: u64 = (spin_chain.vars.tau / spin_chain.vars.dt) as u64;

            spin_chain.log();

            while spin_chain.t < spin_chain.vars.t {
                spin_chain.update();

                if spin_chain.vars.strob {
                    if ((spin_chain.t / spin_chain.vars.dt) as u64).rem_euclid(tau_steps) == 0 {
                        spin_chain.log_ext();
                    }
                } else {
                    spin_chain.log_ext();
                }

                if spin_chain.t.fract() < spin_chain.vars.dt {
                    pb.set_position(spin_chain.t as u64);
                }
            }
            pb.finish_and_clear();
            //            pb.finish_with_message(&format!("Done {}",i));
        });
    }

    //    m.join_and_clear().unwrap();
    m.join().unwrap();

    //Update config
    conf.offset += num as u32;
    fs::write("config.toml", toml::to_string(&conf).unwrap()).unwrap();
    println!("Finished {} runs", num);
}

fn average(conf: &mut Config) {
    println!("Averaging");
    let mut file_prefix = conf.file.clone();
    file_prefix.push_str("*.dat");
    let mut avg_data: Vec<Vec<f64>> = Vec::with_capacity((conf.t / conf.dt) as usize);
    let mut avg_nums: Vec<f64> = Vec::with_capacity((conf.t / conf.dt) as usize);
    let entries = glob(&file_prefix).expect("Failed to find log files");
    //Can this be done more elegantly?
    let entriesc = glob(&file_prefix).expect("Failed to find log files");
    let pb = ProgressBar::new(entriesc.fold(0, |acc, _| acc + 1));
    pb.set_style(
        ProgressStyle::default_bar()
            .template(
                "{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} ({eta})",
            )
            .progress_chars("#>-"),
    );
    for entry in entries {
        match entry {
            Ok(path) => {
                //Load data into memory
                let file = File::open(&path).expect("Unable to open file");
                //                println!("file: {:?}", path.display());
                pb.inc(1);
                let reader = BufReader::new(file);

                //Run through line-by-line
                //Skip first line
                for (i, line) in reader.lines().skip(1).enumerate() {
                    let line = line.unwrap();

                    //                    println!("Line {}: {}", i, &line);
                    // Read line into appropriate vector

                    let line_data: Vec<f64> = line
                        .split(' ')
                        .flat_map(str::parse::<f64>)
                        .collect::<Vec<f64>>();

                    //Check if we have reached this line-number before, add if need be
                    if i == avg_nums.len() {
                        avg_nums.push(1.0);
                        avg_data.push(line_data);
                    } else {
                        //Increase count of this line
                        avg_nums[i] += 1.0;
                        //Update data in avg_data
                        avg_data[i] = avg_data[i]
                            .iter()
                            .zip(line_data.iter())
                            .map(|(x, y)| x + y)
                            .collect::<Vec<f64>>();
                    }
                }
            }
            Err(e) => println!("Error reading file: {:?}", e),
        };
    }

    pb.finish_with_message("Done");
    // Write data to file, and average as we go
    println!("Writing averaged data");
    let ofile = File::create("avg.dat").unwrap();
    for (i, item) in avg_data.iter().enumerate() {
        //Can this be done better?
        //Should create a string and do it this way
        let mut output: String = "".to_owned();
        for j in item {
            let to_add: String = format!(" {}", j / avg_nums[i]);
            output = output + &to_add;
        }

        writeln!(&ofile, "{}", output).unwrap();
    }
}

fn run_mc_magnus(conf: &mut Config) {
    let file = OpenOptions::new()
        .create(true)
        .append(true)
        .open("mc_live.dat")
        .unwrap();

    let m = MultiProgress::new();

    let sty = ProgressStyle::default_bar()
        .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} ({eta})")
        .progress_chars("#>-");
    let (tx, rx) = mpsc::channel();

    let points: usize = conf.mc_points as usize;
    write_in_bg(file, rx);

    // Create storage for the points
    /*
    let mut e_vec: Vec<f64> = Vec::with_capacity(points);
    let mut es_vec: Vec<f64> = Vec::with_capacity(points);
    let mut mx_vec: Vec<f64> = Vec::with_capacity(points);
    let mut my_vec: Vec<f64> = Vec::with_capacity(points);
    let mut mz_vec: Vec<f64> = Vec::with_capacity(points);
    let mut mxt_vec: Vec<f64> = Vec::with_capacity(points);
    let mut myt_vec: Vec<f64> = Vec::with_capacity(points);
    let mut mzt_vec: Vec<f64> = Vec::with_capacity(points);
    */

    println!("Running Monte-Carlo simulation");

    //What if points <= threads? Then generate 1 point per thread and truncate threads

    let threads = match conf.threads {
        _ if conf.threads >= points => points,
        _ => conf.threads,
    };

    let pool = ThreadPool::new(threads);

    let t_points: u32 = (points / threads) as u32;

    for i in 0..threads {
        let mut sc: SpinChain = SpinChain::new(conf.clone(), i);
        let txc = mpsc::Sender::clone(&tx);

        //Make up points in last thread if not divisble
        let t_points = match i {
            _ if i == threads - 1 => t_points + (points as u32) - (threads as u32) * t_points,
            _ => t_points,
        };

        let pb = m.add(ProgressBar::new(t_points as u64));
        pb.set_style(sty.clone());

        pool.execute(move || {
            for _ in 0..t_points as usize {
                for _ in 0..1e6 as usize {
                    sc.metropolis_update_mag_full();
                }
                let m: [f64; 3] = sc.m();
                let mt: [f64; 3] = sc.m_tot();
                let ms: [f64; 3] = [
                    m[0] / sc.vars.ssize as f64,
                    m[1] / sc.vars.ssize as f64,
                    m[2] / sc.vars.ssize as f64,
                ];
                let _mst: [f64; 3] = [
                    mt[0] / sc.vars.hsize as f64,
                    mt[1] / sc.vars.hsize as f64,
                    mt[2] / sc.vars.hsize as f64,
                ];
                let e: f64 = sc.total_energy2();
                let es: f64 = sc.mc_system_energy();

                let s: f64 = sc.vars.ssize as f64;
                let sm: [f64; 3] = sc.spins[(s / 4.0) as usize].dir;
                let sp: [f64; 3] = sc.spins[(3.0 * s / 4.0) as usize].dir;
                let formatted_string: String = format!(
                    "{} {} {} {} {} {} {} {} {} {} {}",
                    e, es, ms[0], ms[1], ms[2], sm[0], sm[1], sm[2], sp[0], sp[1], sp[2]
                );
                txc.send(formatted_string).unwrap();
                pb.inc(1);
            }
            pb.finish_with_message("Done");
        });
    }
    m.join().unwrap();
}

fn run_mc_profile(conf: &mut Config) {
    // Generate x magnetisation profile

    let file = File::create("mc_profile.dat").unwrap();

    let m = MultiProgress::new();

    let sty = ProgressStyle::default_bar()
        .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} ({eta})")
        .progress_chars("#>-");
    let (tx, rx) = mpsc::channel();

    let points: usize = conf.mc_points as usize;
    let writing_thread = profile_in_bg(file, rx, conf.hsize as usize, points);

    // Create storage for the points
    /*
    let mut e_vec: Vec<f64> = Vec::with_capacity(points);
    let mut es_vec: Vec<f64> = Vec::with_capacity(points);
    let mut mx_vec: Vec<f64> = Vec::with_capacity(points);
    let mut my_vec: Vec<f64> = Vec::with_capacity(points);
    let mut mz_vec: Vec<f64> = Vec::with_capacity(points);
    let mut mxt_vec: Vec<f64> = Vec::with_capacity(points);
    let mut myt_vec: Vec<f64> = Vec::with_capacity(points);
    let mut mzt_vec: Vec<f64> = Vec::with_capacity(points);
    */

    println!("Running Monte-Carlo simulation");

    //What if points <= threads? Then generate 1 point per thread and truncate threads

    let threads = match conf.threads {
        _ if conf.threads >= points => points,
        _ => conf.threads,
    };

    let pool = ThreadPool::new(threads);

    let t_points: u32 = (points / threads) as u32;

    for i in 0..threads {
        let mut sc: SpinChain = SpinChain::new(conf.clone(), i);
        let txc = mpsc::Sender::clone(&tx);

        //Make up points in last thread if not divisble
        let t_points = match i {
            _ if i == threads - 1 => t_points + (points as u32) - (threads as u32) * t_points,
            _ => t_points,
        };

        let pb = m.add(ProgressBar::new(t_points as u64));
        pb.set_style(sty.clone());

        pool.execute(move || {
            for _ in 0..t_points as usize {
                for _ in 0..1e6 as usize {
                    sc.metropolis_update();
                }
                for i in 0..sc.vars.hsize as usize {
                    let L: usize = sc.vars.hsize as usize;
                    let s_l: &Spin = match i {
                        _ if i == 0 => &sc.spins[L - 1],
                        _ => &sc.spins[i - 1],
                    };
                    let s_c: &Spin = &sc.spins[i];
                    let s_r: &Spin = match i {
                        _ if i == L - 1 => &sc.spins[0],
                        _ => &sc.spins[i + 1],
                    };
                    let j_l: &[f64; 3] = match i {
                        _ if i == 0 => &sc.j_couple[L - 1],
                        _ => &sc.j_couple[i - 1],
                    };

                    let j_r: &[f64; 3] = &sc.j_couple[i];

                    let mut sjs: f64 = 0.0;
                    for k in 0..3 {
                        sjs -= s_l.dir[k] * j_l[k] * s_c.dir[k];
                        sjs -= s_c.dir[k] * j_r[k] * s_r.dir[k];
                    }
                    sjs = sjs / 2.0;

                    let formatted_string: String = format!(
                        "{} {} {} {} {}",
                        i, sjs, sc.spins[i].dir[0], sc.spins[i].dir[1], sc.spins[i].dir[2]
                    );
                    txc.send(formatted_string).unwrap();
                }
                pb.inc(1);
            }
            pb.finish_with_message("Done");
        });
    }
    m.join().unwrap();
    println!("Joined m threads");
    writing_thread.join().unwrap();
}

fn run_mc(conf: &mut Config) {
    let file = OpenOptions::new()
        .create(true)
        .append(true)
        .open("mc_live.dat")
        .unwrap();

    let m = MultiProgress::new();

    let sty = ProgressStyle::default_bar()
        .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} ({eta})")
        .progress_chars("#>-");
    let (tx, rx) = mpsc::channel();

    let points: usize = conf.mc_points as usize;
    write_in_bg(file, rx);

    // Create storage for the points
    /*
    let mut e_vec: Vec<f64> = Vec::with_capacity(points);
    let mut es_vec: Vec<f64> = Vec::with_capacity(points);
    let mut mx_vec: Vec<f64> = Vec::with_capacity(points);
    let mut my_vec: Vec<f64> = Vec::with_capacity(points);
    let mut mz_vec: Vec<f64> = Vec::with_capacity(points);
    let mut mxt_vec: Vec<f64> = Vec::with_capacity(points);
    let mut myt_vec: Vec<f64> = Vec::with_capacity(points);
    let mut mzt_vec: Vec<f64> = Vec::with_capacity(points);
    */

    println!("Running Monte-Carlo simulation");

    //What if points <= threads? Then generate 1 point per thread and truncate threads

    let threads = match conf.threads {
        _ if conf.threads >= points => points,
        _ => conf.threads,
    };

    let pool = ThreadPool::new(threads);

    let t_points: u32 = (points / threads) as u32;

    for i in 0..threads {
        let mut sc: SpinChain = SpinChain::new(conf.clone(), i);
        let txc = mpsc::Sender::clone(&tx);

        //Make up points in last thread if not divisble
        let t_points = match i {
            _ if i == threads - 1 => t_points + (points as u32) - (threads as u32) * t_points,
            _ => t_points,
        };

        let pb = m.add(ProgressBar::new(t_points as u64));
        pb.set_style(sty.clone());

        pool.execute(move || {
            for _ in 0..t_points as usize {
                for _ in 0..1e6 as usize {
                    sc.metropolis_update();
                }
                let m: [f64; 3] = sc.m();
                let mt: [f64; 3] = sc.m_tot();
                let ms: [f64; 3] = [
                    m[0] / sc.vars.ssize as f64,
                    m[1] / sc.vars.ssize as f64,
                    m[2] / sc.vars.ssize as f64,
                ];
                let _mst: [f64; 3] = [
                    mt[0] / sc.vars.hsize as f64,
                    mt[1] / sc.vars.hsize as f64,
                    mt[2] / sc.vars.hsize as f64,
                ];
                let e: f64 = sc.total_energy2();
                let es: f64 = sc.mc_system_energy();

                let s: f64 = sc.vars.ssize as f64;
                let sm: [f64; 3] = sc.spins[(s / 4.0) as usize].dir;
                let sp: [f64; 3] = sc.spins[(3.0 * s / 4.0) as usize].dir;
                let formatted_string: String = format!(
                    "{} {} {} {} {} {} {} {} {} {} {}",
                    e, es, ms[0], ms[1], ms[2], sm[0], sm[1], sm[2], sp[0], sp[1], sp[2]
                );
                txc.send(formatted_string).unwrap();
                pb.inc(1);
            }
            pb.finish_with_message("Done");
        });
    }
    m.join().unwrap();
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
