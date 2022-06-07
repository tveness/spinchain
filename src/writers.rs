use super::*;

pub fn write_in_bg(file: File, rx: std::sync::mpsc::Receiver<String>) {
    thread::spawn(move || {
        for received in rx {
            writeln!(&file, "{}", &received).unwrap();
        }
    });
}

pub fn profile_in_bg(
    file: File,
    rx: std::sync::mpsc::Receiver<String>,
    size: usize,
    pts: usize,
) -> thread::JoinHandle<()> {
    thread::spawn(move || {
        let mut e_vec_mu: Vec<f64> = vec![0.0; size];
        let mut mx_vec_mu: Vec<f64> = vec![0.0; size];
        let mut my_vec_mu: Vec<f64> = vec![0.0; size];
        let mut mz_vec_mu: Vec<f64> = vec![0.0; size];
        // Issues with hanging to join threads: how to fix?
        //        for received in rx {
        for _ in 0..(size * pts) {
            let received = rx.recv().unwrap();
            let v: Vec<&str> = received.split(' ').collect();
            //            println!("{:?}",v);
            let index: usize = v[0].parse::<usize>().unwrap();
            let e: f64 = v[1].parse::<f64>().unwrap();
            let sx: f64 = v[2].parse::<f64>().unwrap();
            let sy: f64 = v[3].parse::<f64>().unwrap();
            let sz: f64 = v[4].parse::<f64>().unwrap();

            e_vec_mu[index] += e;
            mx_vec_mu[index] += sx;
            my_vec_mu[index] += sy;
            mz_vec_mu[index] += sz;
        }
        println!("Now writing data");
        for i in 0..size {
            let mu_e: f64 = e_vec_mu[i] / (pts as f64);
            let mu_x: f64 = mx_vec_mu[i] / (pts as f64);
            let mu_y: f64 = my_vec_mu[i] / (pts as f64);
            let mu_z: f64 = mz_vec_mu[i] / (pts as f64);
            writeln!(&file, "{} {} {} {} {}", i, mu_e, mu_x, mu_y, mu_z).unwrap();
        }
    })
}
