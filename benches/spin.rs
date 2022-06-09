///Struct for storing a 3-dimensional normalised spin
#[derive(Debug, Clone)]
pub struct Spin {
    pub dir: [f64; 3],
}

impl Spin {
    #[allow(dead_code)]
    pub fn new() -> Spin {
        Spin {
            dir: [1.0, 0.0, 0.0],
        }
    }

    #[allow(dead_code)]
    pub fn new_xyz(input: &[f64]) -> Spin {
        Spin {
            dir: [input[0], input[1], input[2]],
        }
    }

    #[allow(dead_code)]
    pub fn new_from_angles(theta: f64, phi: f64) -> Spin {
        Spin {
            dir: [theta.cos() * phi.sin(), theta.sin() * phi.sin(), phi.cos()],
        }
    }

    #[allow(dead_code)]
    /// Sets the angles of the spin
    pub fn set_angles(&mut self, theta: f64, phi: f64) {
        self.dir = [theta.cos() * phi.sin(), theta.sin() * phi.sin(), phi.cos()];
    }
    #[allow(dead_code)]
    /// Returns a vector of coordinates [x,y,z] of the spin
    pub fn xyz(&self) -> Vec<f64> {
        self.dir.to_vec()
    }

    #[allow(dead_code)]
    pub fn rotate(&mut self, field: &[f64; 3], dt: f64) {
        let fs: f64 = field.iter().map(|x| x * x).sum::<f64>().sqrt();
        let f: [f64; 3] = [field[0] / fs, field[1] / fs, field[2] / fs];
        let dts: f64 = fs * dt;

        let s: &[f64] = &self.xyz();
        let axp: f64 = s.iter().zip(f.iter()).map(|(x, y)| x * y).sum();
        let axo: Vec<f64> = vec![
            f[1] * s[2] - f[2] * s[1],
            f[2] * s[0] - f[0] * s[2],
            f[0] * s[1] - f[1] * s[0],
        ];

        let sxyz: Vec<f64> = (0..3)
            .map(|k| s[k] * dts.cos() + f[k] * axp * (1.0 - dts.cos()) - axo[k] * dts.sin())
            .collect();

        //        eprintln!("Omega . S: {}", self.xyz().iter().zip( field.iter()).map(|(x,y)| x*y).sum::<f64>());
        //        eprintln!("Omega . S': {}", sxyz.iter().zip( field.iter()).map(|(x,y)| x*y).sum::<f64>());
        let ss: f64 = sxyz.iter().map(|x| x * x).sum::<f64>().sqrt();
        self.dir = [sxyz[0] / ss, sxyz[1] / ss, sxyz[2] / ss];
    }
}

impl PartialEq<Spin> for Spin {
    fn eq(&self, other: &Spin) -> bool {
        if (self.dir[0] - other.dir[0]).abs() < 0.001
            && (self.dir[1] - other.dir[1]).abs() < 0.001
            && (self.dir[2] - other.dir[2]).abs() < 0.001
        {
            return true;
        }
        false
    }
}

/*
impl IndexMut<usize> for Vec<&Spin> {
    fn index_mut(&mut self, index: usize) -> &mut Spin {
        &self[index]
    }
}
*/

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn spin_xyz_new() {
        let s1: Spin = Spin::new();
        assert_eq!(s1.xyz(), vec![1.0, 0., 0.]);
    }

    #[test]
    fn spin_xyz_other() {
        let pi: f64 = std::f64::consts::PI;
        let s1: Spin = Spin::new_from_angles(pi / 4.0, pi / 2.0);
        let mut s: f64 = 1.0 / 2.0;
        s = s.sqrt();
        let target: Vec<f64> = vec![s, s, 0.];

        let abs_diff: f64 = s1
            .xyz()
            .iter()
            .zip(target.iter())
            .map(|(x, y)| (x - y).abs())
            .sum();

        assert!(abs_diff < 0.001);
    }

    #[test]
    fn rotation_test_xy_plane() {
        let pi: f64 = std::f64::consts::PI;
        let mut s1: Spin = Spin::new_from_angles(0.0, pi / 2.0);
        let field: [f64; 3] = [0.0, 0.0, -1.0];
        let dt: f64 = 0.01;
        let mut t: f64 = 0.0;
        while t + dt < pi / 2.0 {
            t += dt;
            s1.rotate(&field, dt);
        }
        let abs_diff: f64 = (1.0 - s1.dir[1]).abs();
        println!("t-pi/2: {}", t - pi / 2.0);
        println!("Diff: {}", abs_diff);
        assert!(abs_diff < 0.001);
    }
    #[test]
    fn rotation_test_xz_plane() {
        let pi: f64 = std::f64::consts::PI;
        let mut s1: Spin = Spin::new_from_angles(0.0, pi / 2.0);
        let field: [f64; 3] = [0.0, 1.0, 0.0];
        let dt: f64 = 0.01;
        let mut t: f64 = 0.0;
        while t + dt < pi / 2.0 {
            t += dt;
            s1.rotate(&field, dt);
        }
        let abs_diff: f64 = (1.0 - s1.dir[2]).abs();
        println!("Diff: {}", abs_diff);
        assert!(abs_diff < 0.001);
    }
    #[test]
    fn erotation_test_xz_plane() {
        let pi: f64 = std::f64::consts::PI;
        let mut s1: Spin = Spin::new_from_angles(0.0, pi / 2.0);
        let field: [f64; 3] = [0.0, 1.0, 0.0];
        let dt: f64 = 0.01;
        let mut t: f64 = 0.0;
        while t + dt < pi / 2.0 {
            t += dt;
            s1.rotate(&field, dt);
        }
        let abs_diff: f64 = (1.0 - s1.dir[2]).abs();
        assert!(abs_diff < 0.001);
    }

    #[test]
    fn spin_equality() {
        let pi: f64 = std::f64::consts::PI;
        let s1: Spin = Spin::new_from_angles(0.0, pi / 2.0);
        let s2: Spin = Spin::new_from_angles(0.0, pi / 2.0);
        let s3: Spin = Spin::new_from_angles(0.1, pi / 2.0);
        let s4: Spin = Spin::new_from_angles(0.10000001, pi / 2.0);

        assert_eq!(s1, s2);
        assert_ne!(s1, s3);
        assert_ne!(s2, s3);
        assert_eq!(s3, s4);
    }
}
