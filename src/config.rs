use serde::{Deserialize, Serialize};

#[derive(Deserialize, Serialize, Debug)]
#[serde(default)]
/// This stores a configuration for the classical spin chain
pub struct Config {
    /// the size of the entire system and reservoir
    pub hsize: i32,
    ///is the size of the system
    pub ssize: i32,
    ///the current time of the simulation
    pub t: f32,
    ///timestep
    pub dt: f32,
    ///how many runs to perform
    pub runs: i32,
    ///trel
    pub trel: f32,
    ///relaxation time
    pub tau: f32,
    ///anisotropy parameter
    pub lambda: f32,
    ///static magnetic field [x,y,z]
    pub hfield: Vec<f32>,
    ///variance in spin-spin coupling
    pub jvar: f32,
    ///variance in magnetic field around hfield
    pub hvar: f32,
    ///method to use for Trotter-Suzuki decomposition
    pub method: i32,
    ///target energy density for chain initialisation
    pub ednsty: f32,
    /// Logging file
    pub file: String,
}

impl Default for Config {
    fn default() -> Self {
        Config {
            hsize: 512,
            ssize: 8,
            t: 256.0,
            dt: 0.02,
            runs: 2,
            trel: 0.0,
            tau: 10.0,
            lambda: 1.0,
            hfield: [0.0, 0.0, 0.0].to_vec(),
            jvar: 0.001,
            hvar: 0.1,
            method: 2,
            ednsty: -0.66,
            file: "log.dat".to_string(),
        }
    }
}