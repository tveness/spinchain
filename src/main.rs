use indicatif::{MultiProgress, ProgressBar, ProgressStyle};
use std::fs;
use std::fs::{File, OpenOptions};
use std::io::{prelude::*, BufReader};
use threadpool::ThreadPool;
mod config;
use self::config::Config;

mod spin;
use self::spin::Spin;
use std::f64::consts::PI;

mod spinchain;
use self::spinchain::SpinChain;

use clap::{App, Arg};

use glob::glob;

///Generate ```sample_num``` samples via Monte Carlo and dynamical evolution to check that the
///system indeed thermalises via time-evolution
fn gen_hist(conf: &mut Config, sample_num: usize) {
    conf.drive = false;
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

        writeln!(&mc_hist, "{} {} {} {} {}", i, mx, my, mz, ed).unwrap();

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
        let ed: f64 = sc.system_energy();

        writeln!(&dyn_hist, "{} {} {} {} {}", i, mx, my, mz, ed).unwrap();

        pb.inc(1);
    }
    pb.finish_with_message("Done");
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

            for _ in 0..2e7 as usize {
                spin_chain.metropolis_update();
            }
            spin_chain.log();
            while spin_chain.t < spin_chain.vars.t {
                spin_chain.update();
                if spin_chain.vars.strob {
                    if spin_chain.t.fract() < spin_chain.vars.dt / 2.0 {
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

fn average(conf: &mut Config) {
    println!("Averaging");
    let mut file_prefix = conf.file.clone();
    file_prefix.push_str("*.dat");
    let mut avg_data: Vec<Vec<f64>> = Vec::with_capacity((conf.t / conf.dt) as usize);
    let mut avg_nums: Vec<f64> = Vec::with_capacity((conf.t / conf.dt) as usize);
    let entries = glob(&file_prefix).expect("Failed to find log files");
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
        writeln!(
            &ofile,
            "{} {} {} {} {} {} {}",
            item[0] / avg_nums[i],
            item[1] / avg_nums[i],
            item[2] / avg_nums[i],
            item[3] / avg_nums[i],
            item[4] / avg_nums[i],
            item[5] / avg_nums[i],
            item[6] / avg_nums[i]
        )
        .unwrap();
    }
}

fn run_mc(conf: &mut Config) {
    let file_log = OpenOptions::new()
        .create(true)
        .append(true)
        .open("mc.dat")
        .unwrap();
    let file = File::create("mc_live.dat").unwrap();

    let points: usize = conf.mc_points as usize;
    // Gather 100 points
    let mut e_vec: Vec<f64> = Vec::with_capacity(points);
    let mut es_vec: Vec<f64> = Vec::with_capacity(points);
    let mut mx_vec: Vec<f64> = Vec::with_capacity(points);
    let mut my_vec: Vec<f64> = Vec::with_capacity(points);
    let mut mz_vec: Vec<f64> = Vec::with_capacity(points);
    let mut mxt_vec: Vec<f64> = Vec::with_capacity(points);
    let mut myt_vec: Vec<f64> = Vec::with_capacity(points);
    let mut mzt_vec: Vec<f64> = Vec::with_capacity(points);

    let mut sc: SpinChain = SpinChain::new(conf.clone(), 0);
    println!("Running Monte-Carlo simulation");
    println!("beta: {}", sc.vars.beta);

    let pb = ProgressBar::new(points as u64);
    pb.set_style(
        ProgressStyle::default_bar()
            .template(
                "{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} ({eta})",
            )
            .progress_chars("#>-"),
    );

    for i in 0..points {
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
        let mst: [f64; 3] = [
            mt[0] / sc.vars.hsize as f64,
            mt[1] / sc.vars.hsize as f64,
            mt[2] / sc.vars.hsize as f64,
        ];
        let e: f64 = sc.total_energy2();
        let es: f64 = sc.system_energy();
        e_vec.push(e);
        es_vec.push(es);
        mx_vec.push(ms[0]);
        my_vec.push(ms[1]);
        mz_vec.push(ms[2]);
        mxt_vec.push(mst[0]);
        myt_vec.push(mst[1]);
        mzt_vec.push(mst[2]);
        writeln!(&file, "{} {} {} {} {} {}", i, e, es, ms[0], ms[1], ms[2]).unwrap();
        pb.inc(1);
    }
    pb.finish_with_message("Done");
    println!("Average values");
    let e_avg = e_vec.iter().sum::<f64>() / points as f64;
    let es_avg = es_vec.iter().sum::<f64>() / points as f64;
    let mx_avg = mx_vec.iter().sum::<f64>() / points as f64;
    let my_avg = my_vec.iter().sum::<f64>() / points as f64;
    let mz_avg = mz_vec.iter().sum::<f64>() / points as f64;
    let mxt_avg = mxt_vec.iter().sum::<f64>() / points as f64;
    let myt_avg = myt_vec.iter().sum::<f64>() / points as f64;
    let mzt_avg = mzt_vec.iter().sum::<f64>() / points as f64;
    writeln!(
        &file_log,
        "{} {} {} {} {}",
        sc.vars.beta, e_avg, mx_avg, my_avg, mz_avg
    )
    .unwrap();
    println!("e: {}", e_avg);
    println!("es: {}", es_avg);
    println!("mx: {}", mx_avg);
    println!("my: {}", my_avg);
    println!("mz: {}", mz_avg);
    println!("e-omega mz: {}", e_avg - 2.0 * PI / (sc.vars.tau) * mzt_avg);
    println!(
        "es-omega mz: {}",
        es_avg - 2.0 * PI / (sc.vars.tau) * mz_avg
    );
}

fn main() {
    //Create spin chain with parameters in file "config.toml"

    //Want to read num from file
    let mut conf: Config = SpinChain::read_config("config.toml");

    let matches = App::new("spinchain")
        .version("0.1")
        .author("Thomas Veness <thomas.veness@nottingham.ac.uk>")
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
            Arg::with_name("hist")
            .short("h")
            .value_name("POINTS")
            .long("histogram")
            .takes_value(true)
            .default_value("8000")
            .help("Generate histograms via Monte-Carlo (hist_mc.dat) and via time-evolution (hist_dyn.dat)"),
            )
        .get_matches();
    let mut default = true;

    if let Some(b) = matches.value_of("beta") {
        //        println!("Overriding to: {}", b);
        conf.beta = b.parse::<f64>().unwrap();
    } else {
        conf.beta = SpinChain::solve_beta(conf.ednsty);
    }

    let points: usize = match matches.value_of("hist") {
        Some(x) => x.parse::<usize>().unwrap(),
        None => 8000 as usize,
    };

    match matches.occurrences_of("hist") {
        0 => (),
        _ => {
            println!("Running hist");
            default = false;
            gen_hist(&mut conf, points);
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

    match default {
        true => run_sim(&mut conf),
        false => (),
    }
}
