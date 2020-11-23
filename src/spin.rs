///Struct for storing a 3-dimensional normalised spin
#[derive(Debug, Clone)]
pub struct Spin {
    pub x: f32,
    pub y: f32,
    pub z: f32,
}

impl Spin {
    #[allow(dead_code)]
    pub fn new() -> Spin {
        Spin {
            x: 1.0,
            y: 0.0,
            z: 0.0,
        }
    }

    #[allow(dead_code)]
    pub fn new_xyz(input: &[f32]) -> Spin {
        Spin {
            x: input[0],
            y: input[1],
            z: input[2],
        }
    }

    #[allow(dead_code)]
    pub fn new_from_angles(theta: f32, phi: f32) -> Spin {
        Spin {
            x: theta.cos() * phi.sin(),
            y: theta.sin() * phi.sin(),
            z: phi.cos(),
        }
    }

    #[allow(dead_code)]
    /// Sets the angles of the spin
    pub fn set_angles(&mut self, theta: f32, phi: f32) {
        self.x = theta.cos() * phi.sin();
        self.y = theta.sin() * phi.sin();
        self.z = phi.cos();
    }

    /// Returns a vector of coordinates [x,y,z] of the spin
    pub fn xyz(&self) -> Vec<f32> {
        vec![self.x, self.y, self.z]
    }

    pub fn rotate(&mut self, field: &[f32], dt: f32) {
        let fs: f32 = field.iter().map(|x| x * x).sum::<f32>().sqrt();
        let f: Vec<f32> = field.iter().map(|x| x / fs).collect();
        let dts: f32 = fs * dt;

        let s: &[f32] = &self.xyz();
        let axp: f32 = s.iter().zip(f.iter()).map(|(x, y)| x * y).sum();
        let axo: Vec<f32> = vec![
            f[1] * s[2] - f[2] * s[1],
            f[2] * s[0] - f[0] * s[2],
            f[0] * s[1] - f[1] * s[0],
        ];

        let sxyz: Vec<f32> = (0..3)
            .map(|k| s[k] * dts.cos() + f[k] * axp * (1.0 - dts.cos()) - axo[k] * dts.sin())
            .collect();

//        eprintln!("Omega . S: {}", self.xyz().iter().zip( field.iter()).map(|(x,y)| x*y).sum::<f32>());
//        eprintln!("Omega . S': {}", sxyz.iter().zip( field.iter()).map(|(x,y)| x*y).sum::<f32>());
        let ss: f32 = sxyz.iter().map(|x| x * x).sum::<f32>().sqrt();
        self.x = sxyz[0] / ss;
        self.y = sxyz[1] / ss;
        self.z = sxyz[2] / ss;
    }
}

impl PartialEq<Spin> for Spin {
    fn eq(&self, other: &Spin) -> bool {
        if (self.x - other.x).abs() < 0.001
            && (self.y - other.y).abs() < 0.001
            && (self.z - other.z).abs() < 0.001
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
        let pi: f32 = std::f32::consts::PI;
        let s1: Spin = Spin::new_from_angles(pi / 4.0, pi / 2.0);
        let mut s: f32 = 1.0 / 2.0;
        s = s.sqrt();
        let target: Vec<f32> = vec![s, s, 0.];

        let abs_diff: f32 = s1
            .xyz()
            .iter()
            .zip(target.iter())
            .map(|(x, y)| (x - y).abs())
            .sum();

        assert!(abs_diff < 0.001);
    }

    #[test]
    fn rotation_test_xy_plane() {
        let pi: f32 = std::f32::consts::PI;
        let mut s1: Spin = Spin::new_from_angles(0.0, pi / 2.0);
        let field: Vec<f32> = vec![0.0, 0.0, -1.0];
        let dt: f32 = 0.01;
        let mut t: f32 = 0.0;
        while t + dt < pi / 2.0 {
            t += dt;
            s1.rotate(&field, dt);
        }
        let abs_diff: f32 = (1.0 - s1.y).abs();
        println!("t-pi/2: {}", t - pi / 2.0);
        println!("Diff: {}", abs_diff);
        assert!(abs_diff < 0.001);
    }
    #[test]
    fn rotation_test_xz_plane() {
        let pi: f32 = std::f32::consts::PI;
        let mut s1: Spin = Spin::new_from_angles(0.0, pi / 2.0);
        let field: Vec<f32> = vec![0.0, 1.0, 0.0];
        let dt: f32 = 0.01;
        let mut t: f32 = 0.0;
        while t + dt < pi / 2.0 {
            t += dt;
            s1.rotate(&field, dt);
        }
        let abs_diff: f32 = (1.0 - s1.z).abs();
        println!("Diff: {}", abs_diff);
        assert!(abs_diff < 0.001);
    }
    #[test]
    fn erotation_test_xz_plane() {
        let pi: f32 = std::f32::consts::PI;
        let mut s1: Spin = Spin::new_from_angles(0.0, pi / 2.0);
        let field: Vec<f32> = vec![0.0, 1.0, 0.0];
        let dt: f32 = 0.01;
        let mut t: f32 = 0.0;
        while t + dt < pi / 2.0 {
            t += dt;
            s1.rotate(&field, dt);
        }
        let abs_diff: f32 = (1.0 - s1.z).abs();
        assert!(abs_diff < 0.001);
    }

    #[test]
    fn spin_equality() {
        let pi: f32 = std::f32::consts::PI;
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
