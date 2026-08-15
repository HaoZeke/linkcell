//! Tight k-nearest loop for perf / cachegrind / roofline.
//!
//! Env: `KNN_HOT_N` (default 4096), `KNN_HOT_K` (default 4),
//! `KNN_HOT_REPS` (default 40).

use linkcell::{knearest_into, Cell};

fn env_usize(key: &str, default: usize) -> usize {
    std::env::var(key)
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(default)
}

fn lattice(n: usize, boxl: f64) -> Vec<[f64; 3]> {
    let nside = ((n as f64).cbrt().ceil() as usize).max(1);
    let a = boxl / nside as f64;
    let mut xyz = Vec::with_capacity(n);
    for i in 0..n {
        let ix = i % nside;
        let iy = (i / nside) % nside;
        let iz = i / (nside * nside);
        xyz.push([ix as f64 * a, iy as f64 * a, iz as f64 * a]);
    }
    xyz
}

fn main() {
    let n = env_usize("KNN_HOT_N", 4096);
    let k = env_usize("KNN_HOT_K", 4);
    let reps = env_usize("KNN_HOT_REPS", 40);
    let boxl = 50.0;
    let xyz = lattice(n, boxl);
    let cell = Cell::ortho(boxl, boxl, boxl).expect("box");
    let mut out = vec![-1i32; xyz.len() * k];
    let mut acc = 0i32;
    for _ in 0..reps {
        knearest_into(&xyz, &cell, k, None, Some(3.0), &mut out).expect("knearest");
        acc = acc.wrapping_add(out[0]);
    }
    std::hint::black_box(acc);
    eprintln!("n={} k={} reps={}", xyz.len(), k, reps);
}
