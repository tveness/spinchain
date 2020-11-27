#![allow(dead_code)]
#![allow(unused_variables)]

use super::*;
use rand_distr::{Distribution, Normal};
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
        match fs::read_to_string(filename) {
            Ok(s) => toml::from_str(&s).unwrap(),
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

    /// Creates a new spin chain by reading the configuration variables from a toml file
    ///
    /// The spins
    /// are then initialised in a random configuration at the correct energy density for the case
    /// of isotropic coupling and zero magnetic field
    pub fn new(filename: Option<&str>, num: usize) -> SpinChain {
        //Read configuration file
        //Can make this refer to a default if error and then write a config file
        let mut r = rand::thread_rng();
        let conf: Config = match filename {
            Some(x) => Self::read_config(x),
            None => Config::default(),
        };

        //Initialise spins
        //Draw from distribution with particular energy density
        //        let spins: Vec<Spin> = (0..conf.hsize).map(|x| Spin::new()).collect();

        let pi = std::f64::consts::PI;
        let spin_norm_var: f64 = (1.0 - conf.ednsty * conf.ednsty) / (2.0_f64).powi(12);
        let spin_norm_mean: f64 =
            (conf.ednsty * (spin_norm_var * spin_norm_var / 2.0).exp()).acos();
        let mut theta: f64 = 0.0;

        let spin_normal = Normal::new(spin_norm_mean, spin_norm_var).unwrap();
        let spins: Vec<Spin> = (0..conf.hsize as usize)
            .map(|x| {
                theta += spin_normal.sample(&mut r);
                match x {
                    _ if x % 2 == 0 => Spin::new_xyz(&[theta.cos(), theta.sin(), 0.0]),
                    _ => Spin::new_xyz(&[-theta.cos(), -theta.sin(), 0.0]),
                }
            })
            .collect();

        //Initialise J as [1+N(jvar),1+N(jvar),lamba+N(jvar)]
        let j_normal = Normal::new(0.0, conf.jvar).unwrap();
        //        let v = normal.sample(&mut r);
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
            .map(|x| {
                [
                    conf.hfield[0] + h_normal.sample(&mut r),
                    conf.hfield[1] + h_normal.sample(&mut r),
                    conf.hfield[2] + h_normal.sample(&mut r),
                ]
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
    pub fn log(&self) {
        let e: f64 = self.system_energy();
        let es: f64 = self.total_energy2(true);
        //        let de: f64 = self.de();
        let m: [f64; 3] = self.m();
        let s: f64 = self.vars.ssize as f64;
        writeln!(
            &self.file,
            "{} {} {} {} {} {}",
            self.t,
            es,
            e,
            m[0] / s,
            m[1] / s,
            m[2] / s
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

    pub fn total_energy2(&self, drive: bool) -> f64 {
        let mut e: f64 = 0.0;

        let l: usize = self.vars.hsize as usize / 2;

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
                e -= self.spins[2 * j + 1].xyz()[k]
                    * self.j_couple[2 * j + 1][k]
                    * self.spins[2 * j + 2].xyz()[k];
            }
        }

        for j in 0..l {
            for k in 0..3 {
                e -= self.spins[2 * j].xyz()[k]
                    * self.j_couple[2 * j][k]
                    * self.spins[2 * j + 1].xyz()[k];
                e += self.spins[2 * j].xyz()[k] * hfield[2 * j][k];
                e += self.spins[2 * j + 1].xyz()[k] * hfield[2 * j + 1][k];
            }
        }

        for k in 0..3 {
            e -= self.spins[2 * l - 1].xyz()[k]
                * self.j_couple[2 * l - 1][k]
                * self.spins[0].xyz()[k];
        }

        e / (self.vars.hsize as f64)
    }

    ///Calculate the total energy (with periodic boundary conditions) of the spin chain at the
    ///current time
    pub fn total_energy(&self, drive: bool) -> f64 {
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
        let h_e: [f64; 3] = match drive {
            true => self.h_ext(self.t),
            false => [0.0, 0.0, 0.0],
        };
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

    ///Calculate the energy of just the system proper, ignore boundary terms
    fn system_energy(&self) -> f64 {
        //Size is ssize
        let s: usize = self.vars.ssize as usize;

        //Get energy of s J s terms
        let mut sjs: f64 = 0.0;
        for i in 0..s - 1 {
            sjs -= SpinChain::sjs_energy(&self.spins[i], &self.spins[i + 1], &self.j_couple[i]);
        }

        //Get energy of magnetic field terms
        let h_e: [f64; 3] = self.h_ext(self.t);

        // Field is static field + external field
        let h: Vec<[f64; 3]> = (0..s)
            .map(|x| SpinChain::sum_vec(&self.static_h[x], &h_e))
            .collect();

        let mut sh: f64 = 0.0;
        for (i, item) in h.iter().enumerate() {
            sh += SpinChain::sh_energy(&self.spins[i], &item);
        }

        (sjs + sh) / s as f64
    }

    ///Calculate the magnetisation in direction ```Dir``` for the system proper
    pub fn m(&self) -> [f64; 3] {
        //Determine spins
        let s: usize = self.vars.ssize as usize;
        let spins = self.spins[..s].iter();

        //Now calculate magnetisation for this class of spins
        let result: [f64; 3] = self.spins[..self.vars.ssize as usize]
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
        let pi = std::f64::consts::PI;
        let phi: f64 = 2.0 * pi * t / self.vars.tau;

        [phi.cos(), phi.sin(), 0.0]
    }

    fn j_s(j: &[f64; 3], s: &Spin) -> [f64; 3] {
        [j[0] * s.dir[0], j[1] * s.dir[1], j[2] * s.dir[2]]
    }

    ///Update one time-step using Suzuki-Trotter evolution
    pub fn update(&mut self, drive: bool) {
        // Suzuki-Trotter proceeds in three steps:
        // \Omega_n = J_{n-1} S_{n-1} + J_n S_{n+1} - B_n

        // Calculate field here
        let h_ext: [f64; 3] = match drive {
            true => self.h_ext(self.t + self.vars.dt),
            false => [0.0, 0.0, 0.0],
        };

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

        self.rotate_even(
            &self.even_omega(drive, &h, self.vars.dt / 2.0),
            self.vars.dt / 2.0,
        );
        self.rotate_odd(&self.odd_omega(drive, &h, self.vars.dt / 2.0), self.vars.dt);
        self.rotate_even(
            &self.even_omega(drive, &h, self.vars.dt / 2.0),
            self.vars.dt / 2.0,
        );

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

    fn even_omega(&self, drive: bool, h: &[[f64; 3]], delta_t: f64) -> Vec<[f64; 3]> {
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
    fn odd_omega(&self, drive: bool, h: &[[f64; 3]], delta_t: f64) -> Vec<[f64; 3]> {
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

//extern crate test;
pub fn add_two(a: i32) -> i32 {
    a + 2
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

        assert!((s1doth - f[2]).abs() < 0.0001);
        assert!((s2doth - f[0]).abs() < 0.0001);
        assert!((s3doth - f[1]).abs() < 0.0001);
    }

    #[test]
    fn initialise_energy() {
        let sc: SpinChain = SpinChain::new(None, 0);
        let abs_diff = (sc.vars.ednsty - sc.total_energy(true)).abs();
        assert!(abs_diff < 0.01);
    }

    #[test]
    fn interactions_off() {
        let mut sc: SpinChain = SpinChain::new(None, 0);
        //Turn off interactions
        sc.j_couple = vec![[0.0, 0.0, 0.0]; sc.j_couple.len()];
        //Point all spins along the x axis
        sc.spins = vec![Spin::new_xyz(&[1.0, 0.0, 0.0]); sc.spins.len()];
        sc.static_h = vec![[0.0, 0.0, 1.0]; sc.static_h.len()];
        //Update with just static field, and check that each spin evolves correctly
        for i in 0..100 {
            let abs_diff: f64 = (sc.spins[7].dir[0] - (sc.vars.dt * i as f64).cos()).abs();
            sc.update(false);
            println!("spin[7]: {:?}", sc.spins[7].xyz());
            println!("cos[it]: {}", (sc.vars.dt * i as f64).cos());
            println!("abs diff: {}", abs_diff);
            assert!(abs_diff < 0.001);
        }
    }

    #[test]
    fn static_field_norma() {
        let sc: SpinChain = SpinChain::new(None, 0);
        let hx: f64 = sc.static_h.iter().map(|x| x[0]).sum();
        let hy: f64 = sc.static_h.iter().map(|x| x[1]).sum();
        let hz: f64 = sc.static_h.iter().map(|x| x[2]).sum();
        let hh: f64 = (hx * hx + hy * hy + hz * hz).sqrt() / (sc.vars.hsize as f64);
        let hvar_n: f64 = sc.vars.hvar / (sc.vars.hsize as f64).sqrt();
        println!("hh: {}", hh);
        println!("hv: {}", hvar_n);

        assert!(hh <= 2.0 * 3.0 * hvar_n);
    }
}
