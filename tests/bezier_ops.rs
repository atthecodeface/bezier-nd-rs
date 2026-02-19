//a Imports
use bezier_nd::Num;
use bezier_nd::{BasicBezier, BezierND};
mod utils;
use geo_nd::vector;
use rand::prelude::*;
use utils::{assert_near_equal, assert_near_equal_scale, test_beziers_approx_eq};

fn test_bezier_ops<
    F: Num + From<f32>,
    const D: usize,
    const NP: usize,
    B: BasicBezier<F, [F; D]> + std::fmt::Debug,
>(
    bezier: &B,
    pts: [[F; D]; NP],
) {
    eprintln!("Sum points of b");
    let b_pts_total = bezier
        .control_points()
        .iter()
        .fold([F::ZERO; D], |acc, p| vector::add(acc, p, F::ONE));

    let mut b2 = bezier.clone();
    b2.scale(2.0_f32.into());

    eprintln!("Sum points of 2*b");
    let b2_pts_total = b2
        .control_points()
        .iter()
        .fold([F::ZERO; D], |acc, p| vector::add(acc, p, F::ONE));
    assert_near_equal_scale(&b_pts_total, &b2_pts_total, 2.0_f32);

    eprintln!("b3 = b2 + b; b1 = b3 - b2");
    let mut b3 = b2.clone();
    b3.add(bezier);
    let mut b1 = b3.clone();
    b1.sub(&b2);

    eprintln!("Sum b3");
    let b3_pts_total = b3
        .control_points()
        .iter()
        .fold([F::ZERO; D], |acc, p| vector::add(acc, p, F::ONE));
    assert_near_equal_scale(&b_pts_total, &b3_pts_total, 3.0_f32);

    eprintln!("Sum b1");
    let b1_pts_total = b1
        .control_points()
        .iter()
        .fold([F::ZERO; D], |acc, p| vector::add(acc, p, F::ONE));
    assert_near_equal_scale(&b_pts_total, &b1_pts_total, 1.0_f32);

    eprintln!("Compare b and b1");
    test_beziers_approx_eq(bezier, &b1);

    eprintln!("b -= b1, compare with zero");
    b1.map_pts(&|i, p| vector::sub(pts[i], p, F::ONE));
    for p in b1.control_points() {
        assert_near_equal(p, &[F::ZERO; D]);
    }
}

fn test_ops<
    F: Num + From<f32>,
    R: Rng,
    D: Distribution<f32>,
    const N: usize,
    const NP: usize,
    B: BasicBezier<F, [F; N]> + std::fmt::Debug,
    M: Fn([[F; N]; NP]) -> B,
>(
    rng: &mut R,
    distribution: &D,

    create_bezier: M,
) {
    for _ in 0..10 {
        let pts: [[F; N]; NP] = utils::new_random_point_array(rng, distribution);
        let bezier = create_bezier(pts);
        test_bezier_ops(&bezier, pts);
    }
}

#[test]
fn test_all_bezier_ops() {
    let mut rng = utils::make_random("test_all_bezier_ops_seed");
    let distribution = rand::distr::Uniform::new(-10.0_f32, 10.0).unwrap();

    eprintln!(
        "*************************************************\nTesting 'BezierND<>' type degree 6"
    );
    test_ops::<f32, _, _, 3, 2, _, _>(&mut rng, &distribution, |pts| {
        BezierND::<_, 8, _>::new(&pts)
    });

    eprintln!("*************************************************\nTesting 'Vec' type");
    test_ops::<f32, _, _, 5, 7, _, _>(&mut rng, &distribution, |pts| {
        let b: Vec<_> = pts.into();
        b
    });

    eprintln!("*************************************************\nTesting 'Array'[2] type");
    test_ops::<f32, _, _, 4, 2, _, _>(&mut rng, &distribution, |pts| pts);

    eprintln!("*************************************************\nTesting 'Array'[3] type");
    test_ops::<f32, _, _, 4, 3, _, _>(&mut rng, &distribution, |pts| pts);

    eprintln!("*************************************************\nTesting 'Array'[4] type");
    test_ops::<f32, _, _, 3, 4, _, _>(&mut rng, &distribution, |pts| pts);
}
