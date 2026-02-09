//a Imports
use bezier_nd::{BasicBezier, Bezier, BezierND};
use bezier_nd::{Float, Num};
mod utils;
use geo_nd::vector;
use rand::prelude::*;
use utils::{assert_near_equal, assert_near_equal_scale, test_beziers_approx_eq};

fn test_bezier_ops<
    F: Num,
    const D: usize,
    const NP: usize,
    B: BasicBezier<F, [F; D]> + std::fmt::Debug,
>(
    bezier: &B,
    pts: [[F; D]; NP],
) {
    fn sum_pts<'a, F: Num, const D: usize>(
        total: &'a mut [F; D],
    ) -> impl FnMut(usize, &[F; D]) + 'a {
        *total = [F::ZERO; D];
        move |_n, p| {
            *total = vector::add(*total, p, F::ONE);
        }
    }

    eprintln!("Sum points of b");
    let mut b_pts_total = [F::ZERO; D];
    bezier.for_each_control_point(&mut sum_pts(&mut b_pts_total));

    let mut b2 = bezier.clone();
    b2.scale(2.0_f32.into());

    eprintln!("Sum points of 2*b");
    let mut b2_pts_total = [F::ZERO; D];
    b2.for_each_control_point(&mut sum_pts(&mut b2_pts_total));
    assert_near_equal_scale(&b_pts_total, &b2_pts_total, 2.0_f32);

    eprintln!("b3 = b2 + b; b1 = b3 - b2");
    let mut b3 = b2.clone();
    b3.add(bezier);
    let mut b1 = b3.clone();
    b1.sub(&b2);

    eprintln!("Sum b3");
    let mut b3_pts_total = [F::ZERO; D];
    b3.for_each_control_point(&mut sum_pts(&mut b3_pts_total));
    assert_near_equal_scale(&b_pts_total, &b3_pts_total, 3.0_f32);

    eprintln!("Sum b1");
    let mut b1_pts_total = [F::ZERO; D];
    b1.for_each_control_point(&mut sum_pts(&mut b1_pts_total));
    assert_near_equal_scale(&b_pts_total, &b1_pts_total, 1.0_f32);

    eprintln!("Compare b and b1");
    test_beziers_approx_eq(bezier, &b1);

    eprintln!("b -= b1, compare with zero");
    b1.map_pts(&|i, p| vector::sub(pts[i], p, F::ONE));
    b1.for_each_control_point(&mut |_, p| {
        assert_near_equal(p, &[F::ZERO; D]);
    });
}

fn test_ops<
    F: Float,
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

    eprintln!("*************************************************\nTesting 'Bezier<>' type (from 'try_into()), line");
    test_ops::<f32, _, _, 3, 2, _, _>(&mut rng, &distribution, |pts| {
        let b: Bezier<_, _> = pts.as_ref().try_into().unwrap();
        b
    });

    eprintln!("*************************************************\nTesting 'Bezier<>' type again (from 'into'), quadratic");
    test_ops::<f32, _, _, 2, 3, _, _>(&mut rng, &distribution, |pts| {
        let b: Bezier<_, _> = pts.into();
        b
    });

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
