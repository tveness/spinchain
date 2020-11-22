#![allow(dead_code)]
#![allow(unused_variables)]

use super::*;
use rand_distr::{Distribution, Normal};
use std::fs;
use std::fs::File;
use std::io::prelude::*;
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
    j_couple: Vec<Vec<f32>>,
    ///Static field acting on the system
    static_h: Vec<Vec<f32>>,
    ///Spins
    pub spins: Vec<Spin>,
    ///Current time
    pub t: f32,
    //    spin_b: Vec<Spin>,
    ///Configuration variables for the system
    pub vars: Config,
    ///Log file
    file: std::fs::File,
}

impl SpinChain {
    /// Creates a new spin chain by reading the configuration variables from a toml file
    ///
    /// The spins
    /// are then initialised in a random configuration at the correct energy density for the case
    /// of isotropic coupling and zero magnetic field
    pub fn new(filename: Option<&str>) -> SpinChain {
        let mut r = rand::thread_rng();
        //Read configuration file
        //Can make this refer to a default if error and then write a config file
        let conf: Config = match filename {
            Some(x) => {
                toml::from_str(&fs::read_to_string(x).expect("Could not read config file")).unwrap()
            }
            None => Config::default(),
        };

        //Initialise spins
        //Draw from distribution with particular energy density
        //        let spins: Vec<Spin> = (0..conf.hsize).map(|x| Spin::new()).collect();

        let pi = std::f32::consts::PI;
        let spin_norm_var: f32 = (1.0 - conf.ednsty * conf.ednsty) / (2.0_f32).powi(12);
        let spin_norm_mean: f32 =
            (conf.ednsty * (spin_norm_var * spin_norm_var / 2.0).exp()).acos();
        let mut theta: f32 = 0.0;

        let spin_normal = Normal::new(spin_norm_mean, spin_norm_var).unwrap();
        let spins: Vec<Spin> = (0..conf.hsize as usize)
            .map(|x| {
                theta = (theta + spin_normal.sample(&mut r)) % (2.0 * pi);
                match x {
                    _ if x % 2 == 0 => Spin::new_from_angles(theta, pi / 2.0),
                    _ => Spin::new_from_angles(theta + pi, pi / 2.0),
                }
            })
            .collect();

        //Initialise J as [1+N(jvar),1+N(jvar),lamba+N(jvar)]
        let j_normal = Normal::new(0.0, conf.jvar).unwrap();
        //        let v = normal.sample(&mut r);
        let j_couple: Vec<Vec<f32>> = (0..conf.hsize)
            .map(|x| {
                vec![
                    1.0 + j_normal.sample(&mut r),
                    1.0 + j_normal.sample(&mut r),
                    conf.lambda + j_normal.sample(&mut r),
                ]
            })
            .collect();

        //Initialise static h
        let h_normal = Normal::new(0.0, conf.hvar).unwrap();
        let static_h: Vec<Vec<f32>> = (0..conf.hsize)
            .map(|x| {
                vec![
                    conf.hfield[0] + h_normal.sample(&mut r),
                    conf.hfield[1] + h_normal.sample(&mut r),
                    conf.hfield[2] + h_normal.sample(&mut r),
                ]
            })
            .collect();

        //Log file
        let file = File::create(&conf.file).unwrap();

        SpinChain {
            vars: conf,
            spins,
            j_couple,
            static_h,
            t: 0.0,
            file,
        }
    }
    pub fn log(&self) {
        let e: f32 = self.total_energy();
        writeln!(&self.file, "{} {}", self.t, e).unwrap();
    }

    /// Returns the energy of two spins coupled with j i.e. spin1.j.spin2
    fn sjs_energy(spin1: &Spin, spin2: &Spin, j: &[f32]) -> f32 {
        // Calculate s1 J s2
        let s1_xyz = &spin1.xyz();
        let s2_xyz = &spin2.xyz();
        let mut sjs: f32 = 0.0;
        for i in 0..3 {
            sjs += s1_xyz[i] * j[i] * s2_xyz[i];
        }
        sjs
    }

    ///Returns the energy of a spin in a field i.e. spin.field
    fn sh_energy(spin: &Spin, field: &[f32]) -> f32 {
        spin.xyz()
            .iter()
            .zip(field.iter())
            .map(|(x, y)| x * y)
            .sum()
    }

    ///Calculate the total energy (with periodic boundary conditions) of the spin chain at the
    ///current time
    pub fn total_energy(&self) -> f32 {
        //Size is hsize
        let s: usize = self.vars.hsize as usize;
        let ss: usize = self.vars.ssize as usize;
        let spins = &self.spins[..s];
        let js = &self.j_couple[..s];
        let spins1 = &self.spins[1..s];

        //Get energy of s J s terms
        let mut sjs: f32 = 0.0;
        for i in 0..s - 1 {
            sjs += SpinChain::sjs_energy(&spins[i], &spins1[i], &js[i]);
        }

        //Periodic term
        sjs += SpinChain::sjs_energy(&spins[s - 1], &spins[0], &js[s - 1]);

        //Magnetic field
        //
        //Get energy of magnetic field terms
        let h_e: Vec<f32> = self.h_ext();

        // Field is static field + external field
        let h: Vec<Vec<f32>> = (0..s)
            .map(|x| match x {
                index if index < ss => SpinChain::sum_vec(&self.static_h[index], &h_e),
                index => SpinChain::sum_vec(&self.static_h[index], &[0., 0., 0.]),
            })
            .collect();

        let mut sh: f32 = 0.0;
        for i in 0..s {
            sh += SpinChain::sh_energy(&spins[i], &h[i]);
        }
//        println!("sjs: {}", sjs);
 //       println!("sh: {}", sh);
  //      println!("sjs+sh: {}", sjs+sh);
        if sjs.is_nan() {
            panic!("SJS is NAN, spins are: {:?}",self.spins);
            panic!("SJS is NAN");
        }
        -(sjs + sh) / s as f32
    }

    ///Add two vectors element-wise
    fn sum_vec(a: &[f32], b: &[f32]) -> Vec<f32> {
        a.iter().zip(b.iter()).map(|(x, y)| x + y).collect()
    }

    ///Calculate the energy of just the system proper, ignore boundary terms
    fn system_energy(&self) -> f32 {
        //Size is ssize
        let s: usize = self.vars.ssize as usize;
        let spins = &self.spins[..s];
        let js = &self.j_couple[..s];
        let spins1 = &self.spins[1..s];

        //Get energy of s J s terms
        let mut sjs: f32 = 0.0;
        for i in 0..s - 1 {
            sjs += SpinChain::sjs_energy(&spins[i], &spins1[i], &js[i]);
        }

        //Get energy of magnetic field terms
        let h_e: Vec<f32> = self.h_ext();

        // Field is static field + external field
        let h: Vec<Vec<f32>> = (0..s)
            .map(|x| SpinChain::sum_vec(&self.static_h[x], &h_e))
            .collect();

        let mut sh: f32 = 0.0;
        for i in 0..s {
            sh += SpinChain::sh_energy(&spins[i], &h[i]);
        }

        -(sjs + sh) / s as f32
    }

    ///Calculate the magnetisation in direction ```Dir``` for the system proper
    pub fn m(&self, axis: Dir) -> f32 {
        //Determine spins
        let s: usize = self.vars.ssize as usize;
        let spins = self.spins[..s].iter();

        //Now calculate magnetisation for this class of spins
        let result: f32 = match axis {
            Dir::X => spins.map(|s| s.x).sum(),
            Dir::Y => spins.map(|s| s.y).sum(),
            Dir::Z => spins.map(|s| s.z).sum(),
        };
        result / (s as f32)
    }

    ///Calculate the driving field at the current time
    fn h_ext(&self) -> Vec<f32> {
        let pi = std::f32::consts::PI;
        let phi: f32 = 2.0 * pi * (self.t + self.vars.dt / 2.0) / self.vars.tau;
        let result: Vec<f32> = vec![phi.cos(), phi.sin(), 0.0];
        result
    }

    fn j_s(j: &[f32], s: &Spin) -> Vec<f32> {
        j.iter().zip(s.xyz().iter()).map(|(x, y)| x * y).collect()
    }

    ///Update one time-step using Suzuki-Trotter evolution
    pub fn update(&mut self) {

        self.t += self.vars.dt;

        // Suzuki-Trotter proceeds in three steps:
        // \Omega_n = J_{n-1} S_{n-1} + J_n S_{n+1} - B_n

        self.rotate_even(&self.even_omega(), self.vars.dt / 2.0);
        self.rotate_odd(&self.odd_omega(), self.vars.dt);
        self.rotate_even(&self.even_omega(), self.vars.dt / 2.0);
    }

    fn rotate_even(&mut self, field: &[Vec<f32>], dt: f32) {
        for (i, item) in field.iter().enumerate() {
            self.spins[2 * i].rotate(&item, dt);
        }
    }

    fn rotate_odd(&mut self, field: &[Vec<f32>], dt: f32) {
        let l: usize = self.vars.hsize as usize / 2;
        for (i, item) in field.iter().enumerate() {
            self.spins[2*i+1].rotate(&item, dt);
        }
    }

    fn ssh_sum(l: &[f32], r: &[f32], h: &[f32]) -> Vec<f32> {
        l.iter()
            .zip(r.iter().zip(h.iter()))
            .map(|(x, (y, z))| x + y - z)
            .collect()
    }

    fn even_omega(&self) -> Vec<Vec<f32>> {
        let mut result: Vec<Vec<f32>> = vec![];
        let l: usize = self.spins.len() as usize / 2;
        for n in 0..l {
            let left_s: Vec<f32> = match n {
                _ if n == 0 => self.spins[2 * l - 1].xyz(),
                x => self.spins[2 * x - 1].xyz(),
            };
            let right_s: Vec<f32> = self.spins[2 * n + 1].xyz();
            let h: Vec<f32> = self.static_h[2 * n].clone(); //iter().zip( self.h_ext().iter()).map(|(x,y)| x+y).collect();

            result.push(SpinChain::ssh_sum(&left_s, &right_s, &h));
        }

        result
    }
    fn odd_omega(&self) -> Vec<Vec<f32>> {
        let mut result: Vec<Vec<f32>> = vec![];
        let l: usize = self.spins.len() as usize / 2;
        for n in 0..l {
            let left_s: Vec<f32> = self.spins[2 * n].xyz();
            let right_s: Vec<f32> = match n {
                _ if n == l - 1 => self.spins[0].xyz(),
                _ => self.spins[2 * n + 2].xyz(),
            };
            let h: Vec<f32> = self.static_h[2 * n + 1].clone(); //iter().zip( self.h_ext().iter()).map(|(x,y)| x+y).collect();

            result.push(SpinChain::ssh_sum(&left_s, &right_s, &h));
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
        let v1: Vec<f32> = vec![1.0, 2.0, 3.0, -1., -2.2, -6.5];
        let v2: Vec<f32> = vec![-1.2, 2.0, 3.3, -2., -2.2, 6.8];
        let v3: Vec<f32> = vec![-0.2, 4.0, 6.3, -3., -4.4, 0.3];

        let sum_v1_v2: Vec<f32> = SpinChain::sum_vec(&v1, &v2);

        let abs_diff: f32 = v3
            .iter()
            .zip(sum_v1_v2.iter())
            .map(|(x, y)| (x - y).abs())
            .sum();
        println!("Abs diff: {}", abs_diff);

        assert!(abs_diff < 0.00001);
    }

    #[test]
    fn mag_energy() {
        let f: Vec<f32> = vec![-2., -3., -2.2];

        let pi: f32 = std::f32::consts::PI;
        //Spin along z axis
        let s1: Spin = Spin::new_xyz(&vec![0.0,0.0,1.0]);

        //Spin along x axis
        let mut s2: Spin = Spin::new();
        s2.set_angles(0.0, pi / 2.0);

        //Spin along y axis
        let mut s3: Spin = Spin::new();
        s3.set_angles(pi / 2.0, pi / 2.0);

        let s1doth: f32 = SpinChain::sh_energy(&s1, &f);
        let s2doth: f32 = SpinChain::sh_energy(&s2, &f);
        let s3doth: f32 = SpinChain::sh_energy(&s3, &f);

        assert!((s1doth - f[2]).abs() < 0.0001);
        assert!((s2doth - f[0]).abs() < 0.0001);
        assert!((s3doth - f[1]).abs() < 0.0001);
    }

    #[test]
    fn initialise_energy() {
        let sc: SpinChain = SpinChain::new(None);
        let abs_diff = (sc.vars.ednsty - sc.total_energy()).abs();
        assert!(abs_diff < 0.01);
    }

}
