use criterion::{black_box, criterion_group, criterion_main, Criterion};

mod spin;
use self::spin::Spin;

fn rotate_original(spin: Spin, field: &[f64; 3], dt: f64) -> [f64; 3] {
    let fs: f64 = field.iter().map(|x| x * x).sum::<f64>().sqrt();
    let f: [f64; 3] = [field[0] / fs, field[1] / fs, field[2] / fs];

    let dts: f64 = fs * dt;

    let s: &[f64] = &spin.xyz();
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
    [sxyz[0] / ss, sxyz[1] / ss, sxyz[2] / ss]
}

fn rotate(spin: Spin, field: &[f64; 3], dt: f64) -> [f64; 3] {
    let fs: f64 = field.iter().map(|x| x * x).sum::<f64>().sqrt();
    let f: [f64; 3] = [field[0] / fs, field[1] / fs, field[2] / fs];

    let dts: f64 = fs * dt;
    let dc: f64 = dts.cos();
    let ds: f64 = dts.sin();

    let s: &[f64] = &spin.dir;
    let axp: f64 = s.iter().zip(f.iter()).map(|(x, y)| x * y).sum();
    let axo: [f64; 3] = [
        f[1] * s[2] - f[2] * s[1],
        f[2] * s[0] - f[0] * s[2],
        f[0] * s[1] - f[1] * s[0],
    ];

    let sxyz: [f64; 3] = [
        (s[0] * dc + f[0] * axp * (1.0 - dc) - axo[0] * ds),
        (s[1] * dc + f[1] * axp * (1.0 - dc) - axo[1] * ds),
        (s[2] * dc + f[2] * axp * (1.0 - dc) - axo[2] * ds),
    ];

    sxyz
}

fn fibonacci(n: u64) -> u64 {
    (1..n + 1).product()
}

fn spin_bench(c: &mut Criterion) {
    let mut group = c.benchmark_group("Spin rotate");
    group.bench_function("rotate original", |b| {
        b.iter(|| rotate_original(Spin::new(), black_box(&[2.1, 1.2, 3.7]), 0.88))
    });
    group.bench_function("rotate", |b| {
        b.iter(|| rotate(Spin::new(), black_box(&[2.1, 1.2, 3.7]), 0.88))
    });
    group.finish();
}

criterion_group!(benches, spin_bench);
criterion_main!(benches);
