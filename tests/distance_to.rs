//a Imports
use bezier_nd::{
    Approximation, BasicBezier, Bezier, BezierDistance, BezierEval, BezierMinMax, BezierND,
    BezierSplit,
};
use bezier_nd::Float;
mod utils;
use geo_nd::vector;
use rand::prelude::*;

fn test_distance_to<
    F: Float,
    const D: usize,
    B: BezierEval<F, [F; D]>
        + BezierDistance<F, [F; D]>
        + BezierMinMax<F>
        + BezierSplit
        + Clone
        + std::fmt::Debug,
>(
    bezier: &B,
    test_pts: &[[F; D]],
) {
    let mut bbox_min = [F::ZERO; D];
    let mut bbox_max = [F::ZERO; D];
    bezier.for_each_control_point(&mut |i, pt| {
        if i == 0 {
            vector::set_bbox(&[*pt], &mut bbox_min, &mut bbox_max);
        } else {
            vector::update_bbox(&[*pt], &mut bbox_min, &mut bbox_max);
        }
    });
    for n in 0..D {
        if let Some((_t, c)) = bezier.t_coord_at_min_max(true, n) {
            bbox_max[n] = c;
        }
        if let Some((_t, c)) = bezier.t_coord_at_min_max(false, n) {
            bbox_min[n] = c;
        }
    }
    let eps: F = 1E-3_f32.into();
    for t in utils::float_iter(1000) {
        let pt = bezier.point_at(t);
        let within_bbox = pt
            .iter()
            .zip(bbox_min.iter().zip(bbox_max.iter()))
            .all(|(c, (min, max))| (*c >= *min - eps) && (*c <= *max + eps));
        assert!(
            within_bbox,
            "at t {t}: Point {pt:?} not within bbox {bbox_min:?} {bbox_max:?}"
        );
    }
    vector::expand_bbox(1.0_f32.into(), &mut bbox_min, &mut bbox_max);
    const NUM_APPROX_PTS: usize = 100;
    let num_approx_pts: F = (NUM_APPROX_PTS as f32).into();
    let _approx_t_spacing = F::ONE / num_approx_pts;
    let approx = Approximation::of_regular_t(bezier, NUM_APPROX_PTS, None);

    for t in test_pts {
        let mut pt = [F::ZERO; D];
        for i in 0..D {
            pt[i] = t[i] * bbox_max[i] - (t[i] - 1.0_f32.into()) * bbox_min[i];
        }

        eprintln!("Bezier {bezier:?} pt:{pt:?}");
        let (approx_t, min_dist_sq_to_approx) =
            approx
                .iter_t_pts()
                .fold((F::ZERO, 1E6_f32.into()), |acc, (t, p)| {
                    let d_sq = vector::distance_sq(&p, &pt);
                    if d_sq < acc.1 {
                        (t, d_sq)
                    } else {
                        acc
                    }
                });

        let est = bezier.est_min_distance_sq_to(&pt);
        eprintln!("Pt {pt:?} is approx t/dist_sq {approx_t}/{min_dist_sq_to_approx} from approx [pt {:?}], and estimated at {est} from Bezier {bezier:?}", bezier.point_at(approx_t));
        assert!(est <= min_dist_sq_to_approx * (1.05_f32.into()), "Estimate of distance {est} must be less than the actual distance (to approximation it is {min_dist_sq_to_approx}");
        if let Some((_t, d_sq)) = bezier.t_dsq_closest_to_pt(&pt) {
            // eprintln!("Pt {pt:?} closest on Bezier at {t}/{d_sq} [={:?}], approx closest at {approx_t}/{min_dist_sq_to_approx}, from Bezier {bezier:?}", bezier.point_at(t));
            assert!(d_sq <= min_dist_sq_to_approx + (1E-3_f32).into(), "Actual distance {d_sq} must be less than the closest approximation point of {min_dist_sq_to_approx}");

            // The following check is not correct if the point is equidistant to *two* points on the Bezier
            // assert!(t >= approx_t - approx_t_spacing && t <= approx_t + approx_t_spacing,
            //        "Expected closest point on bezier (supposedly at {t}) to have a `t` value within {approx_t_spacing} of the closest `t` on the Approximation  of {approx_t}");
        }
    }
}

fn test_bezier<
    F: Float,
    R: Rng,
    D: Distribution<f32>,
    const N: usize,
    const FD: usize,
    B: BasicBezier<F, [F; FD]> + BezierDistance<F, [F; FD]> + BezierSplit + Clone + std::fmt::Debug,
    M: Fn([[F; FD]; N]) -> B,
>(
    rng: &mut R,
    distribution: &D,

    create_bezier: M,
) {
    for _ in 0..10 {
        let pts: [[F; FD]; N] = utils::new_random_point_array(rng, distribution);
        let test_pts: [[F; FD]; 100] =
            utils::new_random_point_array(rng, &rand::distr::Uniform::new(0.0_f32, 1.0).unwrap());
        let bezier = create_bezier(pts);
        test_distance_to::<_, _, _>(&bezier, &test_pts);
    }
}
#[test]
fn distance_to() {
    let mut rng = utils::make_random("test_cubics_seed");
    let distribution = rand::distr::Uniform::new(-10.0_f32, 10.0).unwrap();

    eprintln!("*************************************************\nTesting 'Bezier<>' 2D type line (from 'try_into())");
    test_bezier::<f32, _, _, 2, 2, _, _>(&mut rng, &distribution, |pts| {
        let b: Bezier<_, _> = pts.as_ref().try_into().unwrap();
        b
    });

    eprintln!("*************************************************\nTesting 'Bezier<>' 2D type again quad (from 'into')");
    test_bezier::<f32, _, _, 3, 2, _, _>(&mut rng, &distribution, |pts| {
        let b: Bezier<_, _> = pts.into();
        b
    });

    eprintln!("*************************************************\nTesting 'Bezier<>' 2D type again cubic (from 'into')");
    test_bezier::<f32, _, _, 4, 2, _, _>(&mut rng, &distribution, |pts| {
        let b: Bezier<_, _> = pts.into();
        b
    });

    eprintln!("*************************************************\nTesting 'Bezier<>' 4D type line (from 'try_into())");
    test_bezier::<f32, _, _, 2, 4, _, _>(&mut rng, &distribution, |pts| {
        let b: Bezier<_, _> = pts.as_ref().try_into().unwrap();
        b
    });

    eprintln!("*************************************************\nTesting 'Bezier<>' 4D type again quad (from 'into')");
    test_bezier::<f32, _, _, 3, 4, _, _>(&mut rng, &distribution, |pts| {
        let b: Bezier<_, _> = pts.into();
        b
    });

    eprintln!("*************************************************\nTesting 'Bezier<>' 4D type again cubic (from 'into')");
    test_bezier::<f32, _, _, 4, 4, _, _>(&mut rng, &distribution, |pts| {
        let b: Bezier<_, _> = pts.into();
        b
    });

    eprintln!("*************************************************\nTesting 'BezierND<>' type line");
    test_bezier::<f32, _, _, 2, 3, _, _>(&mut rng, &distribution, |pts| {
        BezierND::<_, 6, _>::new(&pts)
    });

    eprintln!("*************************************************\nTesting 'BezierND<>' type quad");
    test_bezier::<f32, _, _, 3, 1, _, _>(&mut rng, &distribution, |pts| {
        BezierND::<_, 6, _>::new(&pts)
    });

    eprintln!("*************************************************\nTesting 'BezierND<>' type cubic");
    test_bezier::<f32, _, _, 4, 2, _, _>(&mut rng, &distribution, |pts| {
        BezierND::<_, 6, _>::new(&pts)
    });

    eprintln!(
        "*************************************************\nTesting 'BezierND<>' type degree 6"
    );
    test_bezier::<f32, _, _, 6, 3, _, _>(&mut rng, &distribution, |pts| {
        BezierND::<_, 8, _>::new(&pts)
    });

    eprintln!("*************************************************\nTesting 'Vec' of line type");
    test_bezier::<f32, _, _, 2, 2, _, _>(&mut rng, &distribution, |pts| {
        let b: Vec<_> = pts.into();
        b
    });

    eprintln!("*************************************************\nTesting 'Vec' of quad type");
    test_bezier::<f32, _, _, 3, 4, _, _>(&mut rng, &distribution, |pts| {
        let b: Vec<_> = pts.into();
        b
    });

    eprintln!("*************************************************\nTesting 'Vec' of cubic type");
    test_bezier::<f32, _, _, 4, 2, _, _>(&mut rng, &distribution, |pts| {
        let b: Vec<_> = pts.into();
        b
    });

    eprintln!("*************************************************\nTesting 'Vec' of degree 5 type");
    test_bezier::<f32, _, _, 6, 2, _, _>(&mut rng, &distribution, |pts| {
        let b: Vec<_> = pts.into();
        b
    });

    eprintln!("*************************************************\nTesting 'Array' of line type");
    test_bezier::<f32, _, _, 2, 2, _, _>(&mut rng, &distribution, |pts| pts);

    eprintln!("*************************************************\nTesting 'Array' of quad type");
    test_bezier::<f32, _, _, 3, 2, _, _>(&mut rng, &distribution, |pts| pts);

    eprintln!("*************************************************\nTesting 'Array' of cubic type");
    test_bezier::<f32, _, _, 4, 2, _, _>(&mut rng, &distribution, |pts| pts);
}
