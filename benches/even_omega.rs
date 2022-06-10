use criterion::{black_box, criterion_group, criterion_main, Criterion};

mod spin;
use self::spin::Spin;
use std::f64::consts::PI;

#[inline]
fn sum_vec(a: &[f64], b: &[f64]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

#[inline]
fn ssh_sum(l: &[f64], r: &[f64], h: &[f64]) -> [f64; 3] {
    /*
    l.iter()
    .zip(r.iter().zip(h.iter()))
    .map(|(x, (y, z))| x + y - z)
    .collect()
    */
    [l[0] + r[0] - h[0], l[1] + r[1] - h[1], l[2] + r[2] - h[2]]
}

/*
#[inline]
fn fibonacci(n: u64) -> u64 {
    (1..n + 1).product()
}
*/

#[inline]
fn h_ext(t: f64) -> [f64; 3] {
    let tau: f64 = 10.0;
    let phi: f64 = 2.0 * PI * t / tau;
    [phi.cos(), phi.sin(), 0.0]
}

#[inline]
fn j_s(j: &[f64; 3], s: &Spin) -> [f64; 3] {
    [j[0] * s.dir[0], j[1] * s.dir[1], j[2] * s.dir[2]]
}

#[inline]
fn even_omega(
    spins: &[Spin],
    static_h: &[[f64; 3]],
    j_couple: &[[f64; 3]],
    delta_t: f64,
) -> Vec<[f64; 3]> {
    let mut result: Vec<[f64; 3]> = vec![];
    let ssize: usize = 10;
    let h_ext: [f64; 3] = h_ext(delta_t);
    let l: usize = spins.len() as usize / 2;

    for n in 0..l {
        // J_{2n-1} S_{2n-1} + J_{2n} S_{2n+1} - B
        let left_s: [f64; 3] = match n {
            _ if n == 0 => j_s(&j_couple[2 * l - 1], &spins[2 * l - 1]),
            x => j_s(&j_couple[2 * x - 1], &spins[2 * x - 1]),
        };

        let right_s: [f64; 3] = j_s(&j_couple[2 * n], &spins[2 * n + 1]);

        let h: [f64; 3] = match n {
            n if 2 * n < ssize as usize => sum_vec(&static_h[2 * n], &h_ext),
            _ => static_h[2 * n],
        };

        result.push(ssh_sum(&left_s, &right_s, &h));
    }

    result
}

#[inline]
fn even_omega_new(
    spins: &[Spin],
    static_h: &[[f64; 3]],
    j_couple: &[[f64; 3]],
    _delta_t: f64,
) -> Vec<[f64; 3]> {
    let mut result: Vec<[f64; 3]> = Vec::with_capacity(spins.len() as usize);
    let ssize: usize = 10;
    //    let h_ext: [f64; 3] = h_ext(delta_t);
    let l: usize = spins.len() as usize / 2;
    let _s: usize = ssize / 2;

    //n=0 loop explicit
    {
        let left_s: [f64; 3] = j_s(&j_couple[2 * l - 1], &spins[2 * l - 1]);
        let right_s: [f64; 3] = j_s(&j_couple[0], &spins[1]);

        result.push(ssh_sum(&left_s, &right_s, &static_h[0]));
    }
    //h_ext drive explicit

    for n in 1..l {
        // J_{2n-1} S_{2n-1} + J_{2n} S_{2n+1} - B
        let left_s: [f64; 3] = j_s(&j_couple[2 * n - 1], &spins[2 * n - 1]);
        let right_s: [f64; 3] = j_s(&j_couple[2 * n], &spins[2 * n + 1]);
        result.push(ssh_sum(&left_s, &right_s, &static_h[2 * n]));
    }

    result
}

fn even_omega_benchmark(c: &mut Criterion) {
    let mut group = c.benchmark_group("Even omega");
    group.sample_size(500);
    let j: [f64; 3] = [1.23, 4.56, 7.89];
    let _s: Spin = Spin::new_xyz(&[1.0, 0.0, 0.0]);
    let spins: Vec<Spin> = vec![Spin::new(); 100];
    //Was j.clone() here before, but clone_on_copy
    let j_couple: Vec<[f64; 3]> = vec![j; 200];
    let static_h: Vec<[f64; 3]> = vec![j; 200];
    //    c.bench_function("j_s ", |b| b.iter(|| j_s(black_box(&j),black_box(&s))));
    group.bench_function("original", |b| {
        b.iter(|| {
            even_omega(
                black_box(&spins),
                black_box(&static_h),
                black_box(&j_couple),
                black_box(0.02),
            )
        })
    });
    group.bench_function("optimised", |b| {
        b.iter(|| {
            even_omega_new(
                black_box(&spins),
                black_box(&static_h),
                black_box(&j_couple),
                black_box(0.02),
            )
        })
    });

    group.finish();
}

criterion_group!(benches, even_omega_benchmark);
criterion_main!(benches);
