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

    //println!("{:?}", spin_chain);

    println!("Mx: {}", spin_chain.m(Dir::X));
    println!("My: {}", spin_chain.m(Dir::Y));
    println!("Mz: {}", spin_chain.m(Dir::Z));

    //spin_chain.start_log("data.dat");
    while spin_chain.t<spin_chain.vars.t {
        spin_chain.update();
        spin_chain.log();
        //spin_chain.log();
    }
    println!("Energy density: {}", spin_chain.total_energy());
    //spin_chain.end_log();
    //    println!("Spins: {:?}", spin_chain.spins);
}
