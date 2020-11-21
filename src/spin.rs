///Struct for storing a 3-dimensional normalised spin
#[derive(Debug, Clone)]
pub struct Spin {
    ///Azimuth angle (angle in x-y plane)
    pub theta: f32,
    ///Zenith angle (angle from z-axis)
    pub phi: f32,
}

impl Spin {
    /// Returns x-coordinate x = sin \phi cos \theta
    pub fn x(&self) -> f32 {
        self.phi.sin() * self.theta.cos()
    }
    /// Returns y-coordinate y = sin \phi sin \theta
    pub fn y(&self) -> f32 {
        self.phi.sin() * self.theta.sin()
    }
    /// Returns z-coordinate z = cos \phi
    pub fn z(&self) -> f32 {
        self.phi.cos()
    }

    /// Returns a new spin pointing along the z-axis
    pub fn new() -> Spin {
        Spin { theta: 0., phi: 0. }
    }

    pub fn new_from_angles(theta: f32, phi: f32) -> Spin {
        Spin { theta, phi }
    }

    /// Sets the angles of the spin
    pub fn set_angles(&mut self, theta: f32, phi: f32) {
        self.theta = theta;
        self.phi = phi;
    }

    /// Returns a vector of coordinates [x,y,z] of the spin
    pub fn xyz(&self) -> Vec<f32> {
        vec![self.x(), self.y(), self.z()]
    }

    pub fn rotate(&mut self, field: &[f32], dt: f32) {
        // d \phi / d t= fy cos theta - fx sin theta
        let dphi: f32 = field[1] * self.theta.cos() - field[0] * self.theta.sin();
        let dtheta: f32 = 
             field[2]
                - field[1] * self.theta.sin() * self.phi.cos() / self.phi.sin()
                - field[0] * self.theta.cos() * self.phi.cos() / self.phi.sin();

        self.phi += dt * dphi;
        self.theta += dt * dtheta;
    }
}

impl PartialEq<Spin> for Spin {
    fn eq(&self, other: &Spin) -> bool {

        if (self.theta-other.theta).abs() < 0.001 && (self.phi-other.phi).abs() < 0.001 {
            return true;
        }
        return false;
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
        assert_eq!(s1.xyz(), vec![0., 0., 1.]);
    }

    #[test]
    fn spin_xyz_other() {
        let pi: f32 = std::f32::consts::PI;
        let s1: Spin = Spin::new_from_angles(pi/4.0,pi/2.0);
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
        let mut s1: Spin = Spin::new_from_angles(0.0,pi/2.0);
        let field: Vec<f32> = vec![0.0,0.0,1.0];
        let dt: f32 = 0.01;
        let mut t: f32 = 0.0;
        while t+dt<pi/2.0 {
            t+=dt;
            s1.rotate(&field,dt);
        }
        let abs_diff: f32 = (pi/2.0 - s1.theta).abs();
        println!("t-pi/2: {}",t-pi/2.0);
        println!("Diff: {}", abs_diff);
        assert!(abs_diff < 0.001);

    }
    #[test]
    fn rotation_test_xz_plane() {
        let pi: f32 = std::f32::consts::PI;
        let mut s1: Spin = Spin::new_from_angles(0.0,pi/2.0);
        let field: Vec<f32> = vec![0.0,-1.0,0.0];
        let dt: f32 = 0.01;
        let mut t: f32 = 0.0;
        while t+dt<pi/2.0 {
            t+=dt;
            s1.rotate(&field,dt);
        }
        let abs_diff: f32 = (0.0 - s1.phi).abs();
        println!("phi: {}",s1.phi);
        println!("Diff: {}", abs_diff);
        assert!(abs_diff < 0.001);

    }

    #[test]
    fn spin_equality() {
        let pi: f32 = std::f32::consts::PI;
        let s1: Spin = Spin::new_from_angles(0.0,pi/2.0);
        let s2: Spin = Spin::new_from_angles(0.0,pi/2.0);
        let s3: Spin = Spin::new_from_angles(0.1,pi/2.0);
        let s4: Spin = Spin::new_from_angles(0.10000001,pi/2.0);

        assert_eq!(s1,s2);
        assert_ne!(s1,s3);
        assert_ne!(s2,s3);
        assert_eq!(s3,s4);
    }



}
