use indicatif::{MultiProgress, ProgressBar, ProgressStyle};
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
    let conf: Config = SpinChain::read_config("config.toml");

    let num: usize = conf.runs as usize;

    for i in 0..num {
        let mut spin_chain: SpinChain = SpinChain::new(Some("config.toml"), i);
        let pb = m.add(ProgressBar::new(spin_chain.vars.t as u64));
        pb.set_style(sty.clone());

        pb.set_message(&format!("Run {}", i));

        let _ = thread::spawn(move || {
            while spin_chain.t < spin_chain.vars.t {
                spin_chain.update(true);
                spin_chain.log();
                pb.set_position(spin_chain.t as u64);
            }
            pb.finish_with_message("Done");
        });
    }

    m.join_and_clear().unwrap();
    println!("Finished {} runs", num);



}
