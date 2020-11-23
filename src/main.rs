use indicatif::{ProgressBar, ProgressStyle};
mod config;
use self::config::Config;

mod spin;
use self::spin::Spin;

mod spinchain;
use self::spinchain::{Dir, SpinChain};

fn main() {
    //Create spin chain with parameters in file "config.toml"
    let mut spin_chain: SpinChain = SpinChain::new(Some("config.toml"));
    //    let mut spin_chain: SpinChain = SpinChain::new(None);
    println!("Running with configuration:");
    println!("{}", toml::to_string(&spin_chain.vars).unwrap());

    println!("Mx: {}", spin_chain.m(Dir::X));
    println!("My: {}", spin_chain.m(Dir::Y));
    println!("Mz: {}", spin_chain.m(Dir::Z));

    let pb = ProgressBar::new(spin_chain.vars.t as u64);
    pb.set_style(
        ProgressStyle::default_bar()
            .template(
                "{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] t={pos}/{len} ({eta})",
            )
            .progress_chars("#>-"),
    );
//    println!("H field: {:?}",spin_chain.static_h);

    while spin_chain.t < spin_chain.vars.t {
        spin_chain.update(true);
        spin_chain.log();
        pb.set_position(spin_chain.t as u64);
    }
    pb.finish_with_message("Done");
    println!("Energy density: {}", spin_chain.total_energy(true));
}
