use serde::{Deserialize, Serialize};

#[derive(Deserialize, Serialize, Debug, Clone)]
pub enum DriveType {
    xyplane,
    uniaxial,
    none,
}

#[derive(Deserialize, Serialize, Debug, Clone)]
#[serde(default)]
/// This stores a configuration for the classical spin chain
pub struct Config {
    /// the size of the entire system and reservoir
    pub hsize: i32,
    ///is the size of the system
    pub ssize: i32,
    ///the current time of the simulation
    pub t: f64,
    ///timestep
    pub dt: f64,
    ///how many runs to perform
    pub runs: i32,
    ///trel
    pub trel: f64,
    ///threads
    pub threads: usize,
    ///relaxation time
    pub tau: f64,
    ///anisotropy parameter
    pub lambda: f64,
    ///static magnetic field [x,y,z]
    pub hfield: Vec<f64>,
    ///static magnetic field on subsystem [x,y,z]
    pub hs: Vec<f64>,
    ///variance in spin-spin coupling
    pub jvar: f64,
    ///variance in magnetic field around hfield
    pub hvar: f64,
    ///method to use for Trotter-Suzuki decomposition
    pub method: i32,
    ///target energy density for chain initialisation
    pub ednsty: f64,
    /// Logging file
    pub file: String,
    /// Stroboscopic on?
    pub strob: bool,
    /// Offset for log
    pub offset: u32,
    /// Inverse temperature for Monte Carlo
    pub beta: f64,
    pub mc_points: usize,
    pub drive: DriveType,
}

impl Default for Config {
    fn default() -> Self {
        Config {
            hsize: 512,
            ssize: 8,
            t: 256.0,
            dt: 0.02,
            runs: 2,
            threads: 2,
            trel: 0.0,
            tau: 10.0,
            lambda: 1.0,
            hfield: [0.0, 0.0, 0.0].to_vec(),
            hs: [0.0, 0.0, 0.0].to_vec(),
            jvar: 0.001,
            hvar: 0.0,
            method: 2,
            ednsty: -0.66,
            file: "log".to_string(),
            strob: false,
            beta: 2.88,
            mc_points: 1000,
            offset: 0,
            drive: DriveType::xyplane,
        }
    }
}
