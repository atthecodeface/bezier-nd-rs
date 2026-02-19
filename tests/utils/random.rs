#![allow(dead_code)]

use bezier_nd::{BezierOps, Num};

use rand::prelude::*;
use rand_chacha::ChaCha8Rng;

pub fn max<F: Num>(a: F, b: F) -> F {
    if a > b {
        a
    } else {
        b
    }
}
pub fn min<F: Num>(a: F, b: F) -> F {
    if a < b {
        a
    } else {
        b
    }
}

/// Make a random number generator given a seed
pub fn make_random(seed_text: &str) -> impl Rng {
    let mut seed = [0_u8; 32];
    let seed_bytes = seed_text.as_bytes();
    for (i, b) in seed_bytes.iter().enumerate() {
        let j = i % 32;
        seed[j] = seed[j] * 13 + (seed[j] >> 4) + *b;
    }
    ChaCha8Rng::from_seed(seed)
}

/// Create a random point using a RNG
pub fn random_point<R: Rng, F: Copy + From<f32>, D: Distribution<f32>, const N: usize>(
    rng: &mut R,
    dist: &D,
) -> [F; N] {
    let mut pt = [(0.0_f32).into(); N];
    for p in pt.iter_mut() {
        *p = dist.sample(rng).into();
    }
    pt
}

/// Update a set of random points
pub fn set_random_point_array<R: Rng, F: Copy + From<f32>, D: Distribution<f32>, const N: usize>(
    rng: &mut R,
    dist: &D,
    pts: &mut [[F; N]],
) {
    for p in pts.iter_mut() {
        *p = random_point(rng, dist);
    }
}

/// Update a Bezier with a set of random points
pub fn set_random_bezier<
    R: Rng,
    F: Num + From<f32>,
    D: Distribution<f32>,
    const N: usize,
    B: BezierOps<F, [F; N]>,
>(
    rng: &mut R,
    dist: &D,
    bezier: &mut B,
) {
    let mut f = |pts| {
        set_random_point_array(rng, dist, pts);
        true
    };
    bezier.map_all_pts(&mut f);
}

/// Make an array of random points
pub fn new_random_point_array<
    R: Rng,
    F: Copy + From<f32>,
    D: Distribution<f32>,
    const M: usize,
    const N: usize,
>(
    rng: &mut R,
    distribution: &D,
) -> [[F; N]; M] {
    let mut pts = [[(0.0_f32).into(); N]; M];
    set_random_point_array(rng, distribution, &mut pts);
    pts
}

/// Make an array of random points
pub fn new_random_point_vec<R: Rng, F: Copy + From<f32>, D: Distribution<f32>, const N: usize>(
    rng: &mut R,
    distribution: &D,
    num_pts: usize,
) -> Vec<[F; N]> {
    let mut pts = vec![[(0.0_f32).into(); N]; num_pts];
    set_random_point_array(rng, distribution, &mut pts);
    pts
}
