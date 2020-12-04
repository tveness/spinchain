use indicatif::{MultiProgress, ProgressBar, ProgressStyle};
use std::fs;
mod config;
use self::config::Config;

mod spin;
use self::spin::Spin;

mod spinchain;
use self::spinchain::SpinChain;

use std::thread;

fn main() {
    //Create spin chain with parameters in file "config.toml"

    let m = MultiProgress::new();
    let sty = ProgressStyle::default_bar()
        .template(
            "{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] t={pos}/{len} ({eta}) {msg}",
        )
        .progress_chars("#>-");

    //Want to read num from file
    let mut conf: Config = SpinChain::read_config("config.toml");

    let num: usize = conf.runs as usize;

    conf.beta = SpinChain::solve_beta(conf.ednsty);
    println!("Energy density: {}", conf.ednsty);
    println!("Effective temperature: {}", conf.beta);

    for i in 0..num {
        let mut spin_chain: SpinChain = SpinChain::new(conf.clone(), i + conf.offset as usize);

        // Initialise at a particular temperature, say T=1

        let pb = m.add(ProgressBar::new(spin_chain.vars.t as u64));
        pb.set_style(sty.clone());

        pb.set_message(&format!("Run {}", i));

        let _ = thread::spawn(move || {
            /*
            for _ in 0..2e7 as usize {
                spin_chain.metropolis_update();
            }
            */
            spin_chain.log();
            while spin_chain.t < spin_chain.vars.t {
                spin_chain.update();
                spin_chain.log();
                pb.set_position(spin_chain.t as u64);
            }
            pb.finish_with_message("Done");
        });
    }

    m.join_and_clear().unwrap();

    //Update config
    conf.offset += num as u32;
    fs::write("config.toml", toml::to_string(&conf).unwrap()).unwrap();
    println!("Finished {} runs", num);
}
