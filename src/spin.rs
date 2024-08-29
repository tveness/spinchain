use std::f64::consts::PI;

///Struct for storing a 3-dimensional normalised spin

#[derive(Debug, Clone)]
pub struct Spin {
    /// [x,y,z] coordinates
    pub dir: [f64; 3],
}

impl Spin {
    /// Basic constructor for a Spin
    ///
    /// # Examples
    ///
    /// ```
    /// # use sc::spin::Spin;
    /// let s: Spin = Spin::new();
    ///
    /// assert_eq!(s.dir, [1.0,0.0,0.0]);
    /// ```
    #[allow(dead_code)]
    pub fn new() -> Spin {
        Spin {
            dir: [1.0, 0.0, 0.0],
        }
    }

    /// Constructs a spin from three Cartesian co-ordinates, and ensures normalisation of the spin
    ///
    /// # Examples
    /// ```
    /// # use sc::spin::Spin;
    /// let s_1: Spin = Spin::new_xyz(&[1.0,0.0,1.0]);
    /// let s_2: Spin = Spin::new_xyz(&[2.0,0.0,2.0]);
    /// assert_eq!(s_1,s_2);
    /// ```
    #[allow(dead_code)]
    pub fn new_xyz(input: &[f64]) -> Spin {
        let norm: f64 = input.iter().map(|x| x * x).sum::<f64>().sqrt();
        Spin {
            dir: [input[0] / norm, input[1] / norm, input[2] / norm],
        }
    }

    /// Constructs a spin from three polar angles
    ///
    /// # Examples
    /// ```
    /// # use sc::spin::Spin;
    /// # use std::f64::consts::PI;
    /// let s_xyz: Spin = Spin::new_xyz(&[1.0,0.0,1.0]);
    /// let s_angles: Spin = Spin::new_from_angles(0.0,PI/4.0);
    /// assert_eq!(s_xyz,s_angles);
    /// ```
    #[allow(dead_code)]
    pub fn new_from_angles(theta: f64, phi: f64) -> Spin {
        Spin {
            dir: [theta.cos() * phi.sin(), theta.sin() * phi.sin(), phi.cos()],
        }
    }

    /// Simple Gram-Schmidt process
    /// Subtracts the `old` component from `perp` and returns a normalised vector in this direction
    /// # Examples
    /// ```
    /// # use sc::spin::Spin;
    /// let old: [f64; 3] = [1.0,0.0,0.0];
    /// let perp: [f64; 3] = [1.0, 1.0, 0.0];
    /// let gs_v: [f64; 3] = Spin::gs(&old, &perp);
    ///
    /// assert_eq!(gs_v, [0.0,1.0,0.0]);
    /// ```
    ///
    /// ```
    /// # use sc::spin::Spin;
    /// let old: [f64; 3] = [2.0,1.0,3.0];
    /// let perp: [f64; 3] = [1.0, 1.0, 0.0];
    /// let gs_v: [f64; 3] = Spin::gs(&old, &perp);
    ///
    /// let dot_product: f64 = old.iter().zip(gs_v.iter()).map(|(x,y)| x*y).sum::<f64>();
    /// assert_eq!(dot_product,0.0);
    /// ```
    pub fn gs(old: &[f64; 3], perp: &[f64; 3]) -> [f64; 3] {
        let dot: f64 = perp[0] * old[0] + perp[1] * old[1] + perp[2] * old[2];
        let on: f64 = old.iter().map(|x| x * x).sum::<f64>();
        let sub: [f64; 3] = [
            perp[0] - dot * old[0] / on,
            perp[1] - dot * old[1] / on,
            perp[2] - dot * old[2] / on,
        ];
        let s: f64 = sub.iter().map(|x| x * x).sum::<f64>().sqrt();

        [sub[0] / s, sub[1] / s, sub[2] / s]
    }

    pub fn rot_perp(old: &[f64; 2], rot_dir: &[f64; 3], rot_angle: f64) -> [f64; 2] {
        let mut old_spin = Spin::new_from_angles(old[0], old[1]);
        //        let oo_spin = old_spin.clone();

        // Get perp dir via Gram-Schmidt (i.e. just subtract off old-component
        let perp_dir: [f64; 3] = Spin::gs(&old_spin.dir, rot_dir);

        /*
        println!(
            "Norm of perp_dir: {}",
            perp_dir[0] * perp_dir[0] + perp_dir[1] * perp_dir[1] + perp_dir[2] * perp_dir[2]
        );
        println!(
            "Perp_dir*old_spin: {}",
            perp_dir[0] * old_spin.dir[0]
                + perp_dir[1] * old_spin.dir[1]
                + perp_dir[2] * old_spin.dir[2]
        );
        */

        old_spin.rotate(&perp_dir, -rot_angle);

        /*
        let ss: f64 = oo_spin.dir[0] * old_spin.dir[0]
            + oo_spin.dir[1] * old_spin.dir[1]
            + oo_spin.dir[2] * old_spin.dir[2];
        println!(
            "acos(s.s'): {}, rot_angle: {}, delta: {}",
            ss.acos(),
            rot_angle,
            ss.acos() - rot_angle
        );
        */

        // [theta, phi]
        [old_spin.theta(), old_spin.phi()]
    }

    pub fn theta(&self) -> f64 {
        (self.dir[1] / self.dir[0]).atan()
            + match self.dir[0] {
                x if x > 0.0 => 0.0,
                _ => PI,
            }
    }

    pub fn phi(&self) -> f64 {
        self.dir[2].acos()
    }

    #[allow(dead_code)]
    /// Sets the angles of the spin
    pub fn set_angles(&mut self, theta: f64, phi: f64) {
        self.dir = [theta.cos() * phi.sin(), theta.sin() * phi.sin(), phi.cos()];
    }

    /// Returns a vector of coordinates [x,y,z] of the spin
    pub fn xyz(&self) -> Vec<f64> {
        self.dir.to_vec()
    }

    /// Rotate spin around direction `field` by an angle determined by norm of `field` times `dt`
    ///
    /// # Examples
    /// ```
    /// # use sc::spin::Spin;
    /// # use std::f64::consts::PI;
    /// let mut s: Spin = Spin::new_xyz(&[1.0,0.0,0.0]);
    /// let field: [f64;3] = [0.0,1.0,0.0];
    /// let dt: f64 = PI/2.0;
    /// s.rotate(&field,dt);
    /// assert_eq!(s.dir[2], 1.0);
    /// ```
    pub fn rotate(&mut self, field: &[f64; 3], dt: f64) {
        let fs: f64 = field.iter().map(|x| x * x).sum::<f64>().sqrt();
        let f: [f64; 3] = [field[0] / fs, field[1] / fs, field[2] / fs];
        let dts: f64 = fs * dt;
        let dc: f64 = dts.cos();
        let ds: f64 = dts.sin();

        let s: &[f64] = &self.dir;
        let axp: f64 = s.iter().zip(f.iter()).map(|(x, y)| x * y).sum();
        let axo: Vec<f64> = vec![
            f[1] * s[2] - f[2] * s[1],
            f[2] * s[0] - f[0] * s[2],
            f[0] * s[1] - f[1] * s[0],
        ];

        let sxyz: [f64; 3] = [
            (s[0] * dc + f[0] * axp * (1.0 - dc) - axo[0] * ds),
            (s[1] * dc + f[1] * axp * (1.0 - dc) - axo[1] * ds),
            (s[2] * dc + f[2] * axp * (1.0 - dc) - axo[2] * ds),
        ];
        let ss: f64 = sxyz.iter().map(|x| x * x).sum::<f64>().sqrt();

        self.dir = [sxyz[0] / ss, sxyz[1] / ss, sxyz[2] / ss]
    }
}

impl Default for Spin {
    fn default() -> Self {
        Spin::new()
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
