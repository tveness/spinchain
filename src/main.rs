use indicatif::{MultiProgress, ProgressBar, ProgressStyle};
use std::fs;
use std::fs::File;
use std::io::{self, prelude::*, BufReader};
use threadpool::ThreadPool;
mod config;
use self::config::Config;

mod spin;
use self::spin::Spin;

mod spinchain;
use self::spinchain::SpinChain;

use clap::{App, Arg};

use glob::glob;

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

    conf.beta = SpinChain::solve_beta(conf.ednsty);
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
                spin_chain.log();
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
    let mut avg_data: Vec<Vec<f64>> = vec![];
    let mut avg_nums: Vec<f64> = vec![];

    for entry in glob(&file_prefix).expect("Failed to find log files") {
        match entry {
            Ok(path) => {
                //Load data into memory
                let file = File::open(&path).expect("Unable to open file");
                println!("file: {:?}", path.display());
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

    // Write data to file, and average as we go
    println!("Writing averaged data");
    let ofile = File::create("avg.dat").unwrap();
    for (i, item) in avg_data.iter().enumerate() {
        writeln!(
            &ofile,
            "{} {} {} {} {} {}",
            item[0] / avg_nums[i],
            item[1] / avg_nums[i],
            item[2] / avg_nums[i],
            item[3] / avg_nums[i],
            item[4] / avg_nums[i],
            item[5] / avg_nums[i]
        )
        .unwrap();
    }
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
        .get_matches();

    match matches.is_present("avg") {
        true => average(&mut conf),
        false => run_sim(&mut conf),
    };
}
