#![allow(dead_code)]
#![allow(unused_variables)]

use super::*;
use rand::Rng;
use rand_distr::{Distribution, Normal, Uniform};
use std::f64::consts::PI;
use std::fs;
use std::fs::File;
use std::io::Write;

///Cartesian direction
pub enum Dir {
    X,
    Y,
    Z,
}

#[derive(Debug)]
///Stores spins, couplings, fields, and allows dynamical updates via a Suzuki-Trotter decomposition
pub struct SpinChain {
    ///Coupling matrices for spins
    j_couple: Vec<[f64; 3]>,
    ///Static field acting on the system
    pub static_h: Vec<[f64; 3]>,
    ///Spins
    pub spins: Vec<Spin>,
    ///Current time
    pub t: f64,
    //    spin_b: Vec<Spin>,
    ///Configuration variables for the system
    pub vars: Config,
    ///Log file
    file: std::fs::File,
}

impl SpinChain {
    /// Reads a configuration from a file.
    /// If the file cannot be opened, write a config.toml with default values and use default
    /// values
    pub fn read_config(filename: &str) -> Config {
        //Read files
        match fs::read_to_string(filename) {
            //If okay, parse
            Ok(s) => match toml::from_str(&s) {
                Ok(x) => x,
                Err(e) => {
                    //Panic if file is invalid
                    panic!("Invalid configuration: {:?}", e);
                    //println!("Invalid configuration: {:?}", e);
                    //println!("Using default configuration instead");
                    //Config::default()
                }
            },
            Err(e) => {
                eprintln!(
                    "Problem reading {}, written default configuration to config.toml and using this",
                    filename
                );
                fs::write("config.toml", toml::to_string(&Config::default()).unwrap()).unwrap();
                Config::default()
            }
        }
    }

    pub fn solve_beta(e: f64) -> f64 {
        let mut beta: f64 = 20.0;
        let mut beta2: f64 = 10.0;

        while (beta - beta2).abs() > 0.01 {
            beta2 = beta;
            beta = 1.0 / (e + 1.0 / beta.tanh());
        }
        beta
    }

    /// Creates a new spin chain by reading the configuration variables from a toml file
    ///
    /// The spins
    /// are then initialised in a random configuration at the correct energy density for the case
    /// of isotropic coupling and zero magnetic field
    pub fn new(conf: Config, num: usize) -> SpinChain {
        //Read configuration file
        //Can make this refer to a default if error and then write a config file
        let mut r = rand::thread_rng();

        //Initialise spins
        //Draw from distribution with particular energy density
        //        let spins: Vec<Spin> = (0..conf.hsize).map(|x| Spin::new()).collect();
        let beta_eff: f64 = conf.beta;
        let b_e: f64 = beta_eff.exp();
        let b_ei: f64 = 1.0 / (b_e);

        let pi = std::f64::consts::PI;
        let mut theta: f64 = 0.01;
        let mut phi: f64 = PI / 2.0;

        let unif = Uniform::from(0.0..1.0);
        let spins: Vec<Spin> = (0..conf.hsize as usize)
            .map(|x| {
                // Choose rotation direction randomly from sphere

                let theta_perp: f64 = 2.0 * PI * unif.sample(&mut r);
                let phi_perp: f64 = (1.0 - 2.0 * unif.sample(&mut r)).acos();

                let rot_dir: [f64; 3] = [
                    theta_perp.cos() * phi_perp.sin(),
                    theta_perp.sin() * phi_perp.sin(),
                    phi_perp.cos(),
                ];

                // Inverse transform sample \theta from the distribution P(\theta) =  \sin \theta e^{\beta\cos\theta}

                let rot_angle: f64 =
                    ((b_e - unif.sample(&mut r) * (b_e - b_ei)).ln() / beta_eff).acos();

                // Rotate by \theta about random direction
                let n: [f64; 2] = Spin::rot_perp(&[theta, phi], &rot_dir, rot_angle);

                theta = n[0];
                phi = n[1];

                Spin::new_xyz(&[theta.cos() * phi.sin(), theta.sin() * phi.sin(), phi.cos()])
            })
            .collect();

        //Initialise J as [1+N(jvar),1+N(jvar),lambda+N(jvar)]
        let j_normal = Normal::new(0.0, conf.jvar).unwrap();

        let j_couple: Vec<[f64; 3]> = (0..conf.hsize)
            .map(|x| {
                [
                    1.0 + j_normal.sample(&mut r),
                    1.0 + j_normal.sample(&mut r),
                    conf.lambda + j_normal.sample(&mut r),
                ]
            })
            .collect();

        //Initialise static h
        let h_normal = Normal::new(0.0, conf.hvar).unwrap();
        let static_h: Vec<[f64; 3]> = (0..conf.hsize)
            .map(|x| match x {
                _ if x < conf.ssize => [
                    conf.hfield[0] + conf.hs[0] + h_normal.sample(&mut r),
                    conf.hfield[1] + conf.hs[1] + h_normal.sample(&mut r),
                    conf.hfield[2] + conf.hs[2] + h_normal.sample(&mut r),
                ],
                _ => [
                    conf.hfield[0] + h_normal.sample(&mut r),
                    conf.hfield[1] + h_normal.sample(&mut r),
                    conf.hfield[2] + h_normal.sample(&mut r),
                ],
            })
            .collect();
        let conf_file_name: &str =
            &[conf.file.clone(), num.to_string(), ".dat".to_string()].join("");

        //Log file
        let file = File::create(&conf_file_name).unwrap();
        writeln!(&file, "#t E_total E_sub M_x M_y M_z").unwrap();
        let t0: f64 = -conf.trel;

        SpinChain {
            vars: conf,
            spins,
            j_couple,
            static_h,
            t: t0,
            file,
        }
    }

    /// Does a Metropolis update for a single spin with inverse temperature beta
    pub fn metropolis_update(&mut self) {
        let mut rng = rand::thread_rng();
        let unif = Uniform::from(0.0..1.0);

        // Select a spin at random from the chain
        let s: usize = self.vars.hsize as usize;
        let index: usize = rng.gen_range(0, s) as usize;

        let local_field: [f64; 3] = self.static_h[index];

        //Calculate initial energy
        //First get left and right spins and couplings
        let s_c: &Spin = &self.spins[index];
        let s_l: &Spin = match index {
            _ if index == 0 => &self.spins[s - 1],
            _ => &self.spins[index - 1],
        };

        let s_r: &Spin = match index {
            _ if index == s - 1 => &self.spins[0],
            _ => &self.spins[index + 1],
        };

        let j_l: &[f64; 3] = match index {
            _ if index == 0 => &self.j_couple[s - 1],
            _ => &self.j_couple[index - 1],
        };

        let j_r: &[f64; 3] = &self.j_couple[index];

        let mut ei: f64 = 0.0;
        for k in 0..3 {
            ei -= s_l.dir[k] * j_l[k] * s_c.dir[k] + s_r.dir[k] * j_r[k] * s_c.dir[k];
            ei += s_c.dir[k] * local_field[k];
        }

        // Select random angles
        let theta: f64 = 2.0 * PI * unif.sample(&mut rng);
        let phi: f64 = (1.0 - 2.0 * unif.sample(&mut rng)).acos();
        let new_spin: Spin = Spin::new_from_angles(theta, phi);

        let mut de: f64 = -ei;
        for k in 0..3 {
            de -= s_l.dir[k] * j_l[k] * new_spin.dir[k] + s_r.dir[k] * j_r[k] * new_spin.dir[k];
            de += new_spin.dir[k] * local_field[k];
        }

        // Evaluate energy difference

        // Accept/reject
        // Accept if <0, or w/ prob e^{-\beta \Delta E} otherwise
        if de < 0.0 || unif.sample(&mut rng) < (-self.vars.beta * de).exp() {
            self.spins[index].dir = new_spin.dir;
        }
    }

    //Could possibly do this with a function pointer for the energy, instead?
    //For now, just get it working
    //This calculates the first correction of the Magnus expansion for a monochromatic drive with
    //period omega=2*pi/tau
    //The first term is obviously the time-average i.e. only time-independent terms
    //
    //The next term is
    //
    //-1/(2 omega) * \sum_{j=1}^\ell S^z_j 2 * y .

    pub fn metropolis_update_magnus(&mut self) {
        let mut rng = rand::thread_rng();
        let unif = Uniform::from(0.0..1.0);

        // Select a spin at random from the chain
        let s: usize = self.vars.hsize as usize;
        let index: usize = rng.gen_range(0, s) as usize;

        let local_field: [f64; 3] = self.static_h[index];

        //Calculate initial energy
        //First get left and right spins and couplings
        let s_j: &Spin = &self.spins[index];

        let s_jm1: &Spin = match index {
            _ if index == 0 => &self.spins[s - 1],
            _ => &self.spins[index - 1],
        };

        let s_jp1: &Spin = match index {
            _ if index == s - 1 => &self.spins[0],
            _ => &self.spins[index + 1],
        };

        let s_jp2: &Spin = match index {
            _ if index == s - 2 => &self.spins[0],
            _ if index == s - 1 => &self.spins[1],
            _ => &self.spins[index + 2],
        };
        let s_jm2: &Spin = match index {
            _ if index == 0 => &self.spins[s - 2],
            _ if index == 1 => &self.spins[s - 1],
            _ => &self.spins[index - 2],
        };

        let j_jm1: &[f64; 3] = match index {
            _ if index == 0 => &self.j_couple[s - 1],
            _ => &self.j_couple[index - 1],
        };

        let j_jp1: &[f64; 3] = match index {
            _ if index == s - 1 => &self.j_couple[0],
            _ => &self.j_couple[index + 1],
        };
        let j_jm2: &[f64; 3] = match index {
            _ if index == 0 => &self.j_couple[s - 2],
            _ if index == 1 => &self.j_couple[s - 1],
            _ => &self.j_couple[index - 2],
        };

        let j_j: &[f64; 3] = &self.j_couple[index];

        //Calculate energy before randomising a spin
        let mut ei: f64 = 0.0;
        for k in 0..3 {
            // S J S terms, i.e. \Omega_1
            // Because j->j+1, flipping j effects *two* terms in the energy
            ei -= s_jm1.dir[k] * j_jm1[k] * s_j.dir[k] + s_j.dir[k] * j_j[k] * s_jp1.dir[k];
        }

        //-(1)/(2 \omega) (  S^z + 2y . (J_j S_{j+1} + J_{j-1} S_{j-1}) x S_j )
        // Because j-1->j->j+1, flipping j effects *three* terms in the energy
        // But only if the field is non-zero at site j
        if index < self.vars.ssize as usize {
            ei -= PI / self.vars.tau
                * (s_j.dir[2]
                    + 2.0
                        * ((j_j[2] * s_jp1.dir[2] + j_jm1[2] * s_jm1.dir[2]) * s_j.dir[0]
                            - (j_j[0] * s_jp1.dir[0] + j_jm1[0] * s_jm1.dir[0]) * s_j.dir[2]
                            + match index {
                                _ if index == 0 => 0.0,
                                _ => {
                                    (j_jm1[2] * s_j.dir[2] + j_jm2[2] * s_jm2.dir[2]) * s_jm1.dir[0]
                                        - (j_jm1[0] * s_j.dir[0] + j_jm2[0] * s_jm2.dir[0])
                                            * s_jm1.dir[2]
                                }
                            }
                            + match index {
                                _ if index == (self.vars.ssize - 1) as usize => 0.0,
                                _ => {
                                    (j_jp1[2] * s_jp2.dir[2] + j_j[2] * s_j.dir[2]) * s_jp1.dir[0]
                                        - (j_jp1[0] * s_jp2.dir[0] + j_j[0] * s_j.dir[0])
                                            * s_jp1.dir[2]
                                }
                            }));
        }

        // Select random angles
        let theta: f64 = 2.0 * PI * unif.sample(&mut rng);
        let phi: f64 = (1.0 - 2.0 * unif.sample(&mut rng)).acos();
        let new_spin: Spin = Spin::new_from_angles(theta, phi);

        //Calculate difference in energy with randomised spin
        let mut de: f64 = -ei;
        for k in 0..3 {
            de -=
                s_jm1.dir[k] * j_jm1[k] * new_spin.dir[k] + new_spin.dir[k] * j_j[k] * s_jp1.dir[k];
        }

        if index < self.vars.ssize as usize {
            de -= PI / self.vars.tau
                * (s_j.dir[2]
                    + 2.0
                        * ((j_j[2] * s_jp1.dir[2] + j_jm1[2] * s_jm1.dir[2]) * new_spin.dir[0]
                            - (j_j[0] * s_jp1.dir[0] + j_jm1[0] * s_jm1.dir[0]) * new_spin.dir[2]
                            + match index {
                                _ if index == 0 => 0.0,
                                _ => {
                                    (j_jm1[2] * new_spin.dir[2] + j_jm2[2] * s_jm2.dir[2])
                                        * s_jm1.dir[0]
                                        - (j_jm1[0] * new_spin.dir[0] + j_jm2[0] * s_jm2.dir[0])
                                            * s_jm1.dir[2]
                                }
                            }
                            + match index {
                                _ if index == (self.vars.ssize - 1) as usize => 0.0,
                                _ => {
                                    (j_jp1[2] * s_jp2.dir[2] + j_j[2] * new_spin.dir[2])
                                        * s_jp1.dir[0]
                                        - (j_jp1[0] * s_jp2.dir[0] + j_j[0] * new_spin.dir[0])
                                            * s_jp1.dir[2]
                                }
                            }));
        }
        // Evaluate energy difference

        // Accept/reject
        // Accept if <0, or w/ prob e^{-\beta \Delta E} otherwise
        if de < 0.0 || unif.sample(&mut rng) < (-self.vars.beta * de).exp() {
            self.spins[index].dir = new_spin.dir;
        }
    }

    pub fn log(&self) {
        let e: f64 = self.system_energy();
        let es: f64 = self.total_energy2();
        //        let de: f64 = self.de();
        let m: [f64; 3] = self.m();
        let mt: [f64; 3] = self.m_tot();
        let s: f64 = self.vars.ssize as f64;

        //Also calculate boundary terms
        let boundary_term: f64 = self.boundary_term();

        writeln!(
            &self.file,
            "{} {} {} {} {} {} {} {}",
            self.t,
            es,
            e,
            m[0] / s,
            m[1] / s,
            m[2] / s,
            mt[2] / self.vars.hsize as f64,
            boundary_term
        )
        .unwrap();
    }

    /// Returns the energy of two spins coupled with j i.e. spin1.j.spin2
    fn sjs_energy(spin1: &Spin, spin2: &Spin, j: &[f64]) -> f64 {
        // Calculate s1 J s2
        let s1_xyz = spin1.xyz();
        let s2_xyz = spin2.xyz();
        let mut sjs: f64 = 0.0;
        for i in 0..3 {
            sjs += s1_xyz[i] * j[i] * s2_xyz[i];
        }
        sjs
    }

    #[allow(clippy::needless_range_loop)]
    ///Returns the energy of a spin in a field i.e. spin.field
    fn sh_energy(spin: &Spin, field: &[f64; 3]) -> f64 {
        let mut sh: f64 = 0.0;
        for i in 0..3 {
            sh += spin.dir[i] * field[i];
        }
        sh
    }

    pub fn total_energy2(&self) -> f64 {
        let mut e: f64 = 0.0;

        let l: usize = self.vars.hsize as usize;

        let h_e: [f64; 3] = self.h_ext(self.t);

        let hfield: Vec<[f64; 3]> = (0..self.vars.hsize as usize)
            .map(|x| match x {
                x if x < self.vars.ssize as usize => Self::sum_vec(&self.static_h[x], &h_e),
                _ => self.static_h[x],
            })
            .collect();
        //        println!("h field: {:?} at time t: {}", hfield, self.t);

        for j in 0..l - 1 {
            for k in 0..3 {
                e -= self.spins[j].dir[k] * self.j_couple[j][k] * self.spins[j + 1].dir[k];
                e += self.spins[j].dir[k] * hfield[j][k];
            }
        }

        for k in 0..3 {
            e -= self.spins[l - 1].dir[k] * self.j_couple[l - 1][k] * self.spins[0].dir[k];
            e += self.spins[l - 1].dir[k] * hfield[l - 1][k];
        }

        e / (self.vars.hsize as f64)
    }

    ///Calculate the total energy (with periodic boundary conditions) of the spin chain at the
    ///current time
    pub fn total_energy(&self) -> f64 {
        //E = - S_n J_n S_{n+1} + B_n S_n

        //Size is hsize
        let s: usize = self.vars.hsize as usize;
        let ss: usize = self.vars.ssize as usize;

        //Get energy of s_n J_n s_{n+1} terms
        let mut sjs: f64 = 0.0;
        for i in 0..(s - 1) {
            sjs -= SpinChain::sjs_energy(&self.spins[i], &self.spins[i + 1], &self.j_couple[i]);
            //            for k in 0..3 {
            //                sjs -=  self.spins[i].xyz()[k] * self.j_couple[i][k] * self.spins[i+1].xyz()[k];
            //            }
        }

        //Periodic term
        sjs -= SpinChain::sjs_energy(&self.spins[s - 1], &self.spins[0], &self.j_couple[s - 1]);
        //        for k in 0..3 {
        //            sjs -=  self.spins[s-1].xyz()[k] * self.j_couple[s-1][k] * self.spins[0].xyz()[k];
        //        }

        //Magnetic field
        // Field is static field + external field

        //Static field on whole system
        let mut sh: f64 = 0.0;
        for i in 0..s {
            sh += SpinChain::sh_energy(&self.spins[i], &self.static_h[i]);
            //            for k in 0..3 {
            //                sh+=self.spins[i].xyz()[k] * self.static_h[i][k];
            //            }
        }

        //Driving field on sub-system
        let h_e: [f64; 3] = self.h_ext(self.t);
        for i in 0..ss {
            sh += SpinChain::sh_energy(&self.spins[i], &h_e);
            //            for k in 0..3 {
            //                sh+=self.spins[i].xyz()[k] * h_e[k];
            //            }
        }

        (sjs + sh) / s as f64
        //        (sjs + sh/2.0) / s as f64 //Kay's code
    }

    ///Add two vectors element-wise
    fn sum_vec(a: &[f64], b: &[f64]) -> [f64; 3] {
        [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
    }

    ///Calculate boundary term to test Fokker-Planck evolution
    pub fn boundary_term(&self) -> f64 {
        let ell: usize = self.vars.ssize as usize;
        let s1sq: f64 = self.spins[0].dir.iter().fold(0.0, |acc, x| acc + x * x);
        let slsq: f64 = self.spins[ell - 1]
            .dir
            .iter()
            .fold(0.0, |acc, x| acc + x * x);
        let h_e: [f64; 3] = self.h_ext(self.t);

        //Use macro to do these things later, for now just use vec
        let om_1: Vec<f64> = (0..3)
            .map(|x| self.j_couple[0][x] * self.spins[1].dir[x] - h_e[x])
            .collect();
        let om_ell: Vec<f64> = (0..3)
            .map(|x| self.j_couple[ell - 2][x] * self.spins[ell - 2].dir[x] - h_e[x])
            .collect();

        let om_1_sq: f64 = om_1.iter().fold(0.0, |acc, x| acc + x * x);
        let om_ell_sq: f64 = om_ell.iter().fold(0.0, |acc, x| acc + x * x);
        let oms1: f64 = (0..3).fold(0.0, |acc, x| acc + om_1[x] * self.spins[0].dir[x]);
        let omsell: f64 = (0..3).fold(0.0, |acc, x| acc + om_ell[x] * self.spins[ell - 1].dir[x]);

        s1sq * om_1_sq - oms1 * oms1 + slsq * om_ell_sq - omsell * omsell
    }

    ///Calculate the energy of just the system proper, ignore boundary terms
    pub fn system_energy(&self) -> f64 {
        //Size is ssize
        let s: usize = self.vars.ssize as usize;

        //
        let mut e: f64 = 0.0;

        //Get energy of magnetic field terms
        let h_e: [f64; 3] = self.h_ext(self.t);

        // Field is static field + external field
        let h: Vec<[f64; 3]> = (0..s)
            .map(|x| SpinChain::sum_vec(&self.static_h[x], &h_e))
            .collect();

        for i in 0..s - 1 {
            for k in 0..3 {
                e -= self.spins[i].dir[k] * self.j_couple[i][k] * self.spins[i + 1].dir[k];
                e += self.spins[i].dir[k] * h[i][k];
            }
        }

        e / (s - 1) as f64
    }

    ///Calculate the magnetisation in direction ```Dir``` for the system proper
    pub fn m(&self) -> [f64; 3] {
        //Determine spins
        let s: usize = self.vars.ssize as usize;
        let spins = self.spins[..s].iter();

        //Now calculate magnetisation for this class of spins
        let result: [f64; 3] = self.spins[..s as usize]
            .iter()
            .fold([0.0, 0.0, 0.0], |acc, x| {
                [acc[0] + x.dir[0], acc[1] + x.dir[1], acc[2] + x.dir[2]]
            });
        result
    }

    pub fn m_tot(&self) -> [f64; 3] {
        //Determine spins
        let s: usize = self.vars.hsize as usize;
        let spins = self.spins.iter();

        //Now calculate magnetisation for this class of spins
        let result: [f64; 3] = self.spins[..s as usize]
            .iter()
            .fold([0.0, 0.0, 0.0], |acc, x| {
                [acc[0] + x.dir[0], acc[1] + x.dir[1], acc[2] + x.dir[2]]
            });
        result
    }

    ///Calculate the driving field at the current time
    fn h_ext(&self, t: f64) -> [f64; 3] {
        //            return vec![0.0, 0.0, 0.0];
        //        if t < 0.0 {
        //            return [0.0, 0.0, 0.0];
        //        }
        match self.vars.drive {
            DriveType::xyplane => {
                let pi = std::f64::consts::PI;
                let phi: f64 = 2.0 * pi * t / self.vars.tau;

                [phi.cos(), phi.sin(), 0.0]
            }
            DriveType::uniaxial => {
                let pi = std::f64::consts::PI;
                let phi: f64 = 2.0 * pi * t / self.vars.tau;

                [phi.cos(), 0.0, 0.0]
            }
            DriveType::xyelliptic => {
                let pi = std::f64::consts::PI;
                let phi: f64 = 2.0 * pi * t / self.vars.tau;

                [phi.cos(), self.vars.e * phi.sin(), 0.0]
            }
            DriveType::none => [0.0, 0.0, 0.0],
        }
    }

    fn j_s(j: &[f64; 3], s: &Spin) -> [f64; 3] {
        [j[0] * s.dir[0], j[1] * s.dir[1], j[2] * s.dir[2]]
    }

    ///Update one time-step using Suzuki-Trotter evolution
    pub fn update(&mut self) {
        // Suzuki-Trotter proceeds in three steps:
        // \Omega_n = J_{n-1} S_{n-1} + J_n S_{n+1} - B_n

        // Calculate field here
        let h_ext: [f64; 3] = self.h_ext(self.t + self.vars.dt);

        let s: usize = self.vars.ssize as usize;

        let h: Vec<[f64; 3]> = self
            .static_h
            .iter()
            .enumerate()
            .map(|(i, item)| match i {
                i if i < s => Self::sum_vec(item, &h_ext),
                _ => *item,
            })
            .collect();

        self.rotate_even(&self.even_omega(&h, self.vars.dt / 2.0), self.vars.dt / 2.0);
        self.rotate_odd(&self.odd_omega(&h, self.vars.dt / 2.0), self.vars.dt);
        self.rotate_even(&self.even_omega(&h, self.vars.dt / 2.0), self.vars.dt / 2.0);

        self.t += self.vars.dt;
    }

    ///Rotates every even site by the field for a time dt
    fn rotate_even(&mut self, field: &[[f64; 3]], dt: f64) {
        for (i, item) in field.iter().enumerate() {
            self.spins[2 * i].rotate(&item, dt);
        }
    }

    ///Rotates every odd site by the field for a time dt
    fn rotate_odd(&mut self, field: &[[f64; 3]], dt: f64) {
        let l: usize = self.vars.hsize as usize / 2;
        for (i, item) in field.iter().enumerate() {
            self.spins[2 * i + 1].rotate(&item, dt);
        }
    }

    fn ssh_sum(l: &[f64], r: &[f64], h: &[f64]) -> [f64; 3] {
        /*
        l.iter()
        .zip(r.iter().zip(h.iter()))
        .map(|(x, (y, z))| x + y - z)
        .collect()
        */
        [l[0] + r[0] - h[0], l[1] + r[1] - h[1], l[2] + r[2] - h[2]]
    }

    fn even_omega(&self, h: &[[f64; 3]], delta_t: f64) -> Vec<[f64; 3]> {
        let l: usize = self.spins.len() as usize / 2;
        let mut result: Vec<[f64; 3]> = Vec::with_capacity(l);
        // J_{2n-1} S_{2n-1} + J_{2n} S_{2n+1} - B
        //Do n=0 explicitly
        {
            let left_s: [f64; 3] = Self::j_s(&self.j_couple[2 * l - 1], &self.spins[2 * l - 1]);

            let right_s: [f64; 3] = Self::j_s(&self.j_couple[0], &self.spins[1]);
            result.push(SpinChain::ssh_sum(&left_s, &right_s, &h[0]));
        }

        for n in 1..l {
            let left_s: [f64; 3] = Self::j_s(&self.j_couple[2 * n - 1], &self.spins[2 * n - 1]);
            let right_s: [f64; 3] = Self::j_s(&self.j_couple[2 * n], &self.spins[2 * n + 1]);

            result.push(SpinChain::ssh_sum(&left_s, &right_s, &h[2 * n]));
        }

        result
    }
    fn odd_omega(&self, h: &[[f64; 3]], delta_t: f64) -> Vec<[f64; 3]> {
        let l: usize = self.spins.len() as usize / 2;
        let mut result: Vec<[f64; 3]> = Vec::with_capacity(l);
        // J_{2n} S_{2n} + J_{2n+1} S_{2n+2}

        for n in 0..l - 1 {
            let left_s: [f64; 3] = Self::j_s(&self.j_couple[2 * n], &self.spins[2 * n]);
            let right_s: [f64; 3] = Self::j_s(&self.j_couple[2 * n + 1], &self.spins[2 * n + 2]);

            result.push(SpinChain::ssh_sum(&left_s, &right_s, &h[2 * n + 1]));
        }
        //Do last site explicitly
        {
            let left_s: [f64; 3] = Self::j_s(&self.j_couple[2 * (l - 1)], &self.spins[2 * (l - 1)]);
            let right_s: [f64; 3] = Self::j_s(&self.j_couple[2 * (l - 1) + 1], &self.spins[0]);
            result.push(SpinChain::ssh_sum(&left_s, &right_s, &h[2 * (l - 1) + 1]));
        }

        result
    }

    fn pbc_index_odd(n: usize, l: usize) -> usize {
        match n {
            0 => 2 * l - 1,
            _ => 2 * n - 1,
        }
    }
    fn pbc_index_even(n: usize, l: usize) -> usize {
        match n {
            x if x < l => 2 * x,
            _ => 0,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn vec_sum() {
        let v1: [f64; 3] = [1.0, 2.0, 3.0];
        let v2: [f64; 3] = [-1.2, 2.0, 3.3];
        let v3: [f64; 3] = [-0.2, 4.0, 6.3];

        let sum_v1_v2: [f64; 3] = SpinChain::sum_vec(&v1, &v2);

        let abs_diff: f64 = v3
            .iter()
            .zip(sum_v1_v2.iter())
            .map(|(x, y)| (x - y).abs())
            .sum();
        println!("Abs diff: {}", abs_diff);

        assert!(abs_diff < 0.00001);
    }

    #[test]
    fn mag_energy() {
        let f: [f64; 3] = [-2., -3., -2.2];

        let pi: f64 = std::f64::consts::PI;
        //Spin along z axis
        let s1: Spin = Spin::new_xyz(&vec![0.0, 0.0, 1.0]);

        //Spin along x axis
        let mut s2: Spin = Spin::new();
        s2.set_angles(0.0, pi / 2.0);

        //Spin along y axis
        let mut s3: Spin = Spin::new();
        s3.set_angles(pi / 2.0, pi / 2.0);

        let s1doth: f64 = SpinChain::sh_energy(&s1, &f);
        let s2doth: f64 = SpinChain::sh_energy(&s2, &f);
        let s3doth: f64 = SpinChain::sh_energy(&s3, &f);

        assert_almost_eq!(s1doth, f[2], 0.0001);
        assert_almost_eq!(s2doth, f[0], 0.0001);
        assert_almost_eq!(s3doth, f[1], 0.0001);
    }

    #[test]
    fn initialise_energy() {
        let mut e_total: f64 = 0.0;
        let e_target: f64 = Config::default().ednsty;
        let num: f64 = 10.0;
        for i in 0..num as usize {
            let sc: SpinChain = SpinChain::new(Config::default(), 0);
            e_total += sc.total_energy2();
        }
        let e_avg: f64 = e_total / num;
        println!("e_target: {}", e_target);
        println!("e_calc:   {}", e_avg);
        assert_almost_eq!(e_target, e_avg, 0.02);
    }

    #[test]
    fn interactions_off() {
        let mut sc: SpinChain = SpinChain::new(Config::default(), 0);
        //Turn off interactions
        sc.j_couple = vec![[0.0, 0.0, 0.0]; sc.j_couple.len()];
        //Point all spins along the x axis
        sc.spins = vec![Spin::new_xyz(&[1.0, 0.0, 0.0]); sc.spins.len()];
        sc.static_h = vec![[0.0, 0.0, 1.0]; sc.static_h.len()];
        sc.vars.drive = DriveType::none;
        //Update with just static field, and check that each spin evolves correctly
        for i in 0..100 {
            let x_true: f64 = sc.spins[7].dir[0];
            let x_assert: f64 = (sc.vars.dt * i as f64).cos();
            assert_almost_eq!(x_true, x_assert, 0.001);
            sc.update();
        }
    }

    #[test]
    fn static_field_normal() {
        let sc: SpinChain = SpinChain::new(Config::default(), 0);
        let hx: f64 = sc.static_h.iter().map(|x| x[0]).sum();
        let hy: f64 = sc.static_h.iter().map(|x| x[1]).sum();
        let hz: f64 = sc.static_h.iter().map(|x| x[2]).sum();
        let hh: f64 = (hx * hx + hy * hy + hz * hz).sqrt() / (sc.vars.hsize as f64);
        let hvar_n: f64 = sc.vars.hvar / (sc.vars.hsize as f64).sqrt();
        println!("hh: {}", hh);
        println!("hv: {}", hvar_n);

        assert!(hh <= 2.0 * 3.0 * hvar_n);
    }

    //This test is slow: can we make it better?
    //Also, can we make it always pass?
    #[test]
    fn metropolis_low_t() {
        let mut sc: SpinChain = SpinChain::new(Config::default(), 0);
        //Make the chain clean
        sc.j_couple = vec![[1.0, 1.0, 1.0]; sc.vars.hsize as usize];
        sc.static_h = vec![[0.0, 0.0, 0.0]; sc.vars.hsize as usize];
        sc.vars.beta = 1000.0;
        sc.vars.drive = DriveType::none;
        let num: usize = 1_000_000;

        for i in 0..num as usize {
            sc.metropolis_update();
        }
        //Expect energy to be approx -L - B l
        //        let lL: f64 = sc.vars.ssize as f64 / sc.vars.hsize as f64;
        let e: f64 = sc.total_energy2();
        //        println!("Exp: {}", -1.0 - lL);
        let e_exact: f64 = 1.0 / sc.vars.beta - 1.0 / sc.vars.beta.tanh();
        println!("Exp: {}", e);
        println!("Actual: {}", e_exact);
        //        let ediff: f64 = (-1.0 - lL - e).abs();

        assert_almost_eq!(e, e_exact, 30.0 / ((num as f64).sqrt()));
    }
}
