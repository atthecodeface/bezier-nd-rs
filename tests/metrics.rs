use bezier_nd::BezierEval;
use bezier_nd::{metrics, BezierBuilder, BezierND, BezierOps};
mod utils;
use geo_nd::vector;

fn check_metrics_of_beziers<
    const D: usize,
    B: BezierEval<f32, [f32; D]> + BezierOps<f32, [f32; D]> + Clone + std::fmt::Debug,
>(
    b0: &B,
    b1: &B,
) {
    let mut max_d = 0.0;
    let mut max_t = 0.0;
    for t in utils::float_iter(101) {
        let p0 = b0.point_at(t);
        let p1 = b1.point_at(t);
        let dp = (p1[0] - p0[0]).abs();
        if dp > max_d {
            max_d = dp;
            max_t = t;
        }
    }
    eprintln!(
        "at t {max_t} dp {max_d} b0 {:?} b1 {:?}",
        b0.point_at(max_t),
        b1.point_at(max_t)
    );

    let dl_sq = metrics::dl_sq_est(b0, b1, 10000).unwrap();
    let dm_sq = metrics::dm_sq_est(b0, b1, 10000).unwrap();
    let dc_sq = metrics::dc_sq(b0, b1).unwrap();
    let df_sq = metrics::df_sq(b0, b1).unwrap();
    eprintln!("b0-b1 l2 {dl_sq} dm {dm_sq} dc {dc_sq} df {df_sq}",);

    assert!(dl_sq <= dm_sq * 1.1, "dl_sq < dm_sq (with margin)");
    assert!(dm_sq <= dc_sq * 1.1, "dm_sq < dc_sq (with margin)");
    assert!(dc_sq <= df_sq * 1.01, "dc_sq < df_sq (with margin)");
    assert!(
        df_sq <= dc_sq * (b0.num_control_points() as f32) * 1.01,
        "df_sq < N*dc_sq (with margin)"
    );

    let df_sq = metrics::df_sq_from_line(b0);
    let dc_sq = b0.closeness_sq_to_line();
    assert!(dc_sq <= df_sq * 1.01, "dc_sq < df_sq (with margin)");
    assert!(
        df_sq <= dc_sq * (b0.num_control_points() as f32) * 1.01,
        "df_sq < N*dc_sq (with margin)"
    );

    let mut b0_relative_to_line = b0.clone();
    let (l0, l1) = b0.endpoints();
    let l0 = *l0;
    let l1 = *l1;
    b0_relative_to_line.map_pts(&|n, pt| {
        let t = (n as f32) / (b0.degree() as f32);
        vector::sum_scaled(&[*pt, l0, l1], &[1.0, t - 1.0, -t])
    });
    eprintln!("B0 {b0:?} relative {b0_relative_to_line:?}");
    let f_sq = metrics::f_sq(&b0_relative_to_line);
    let c_sq = metrics::c_sq(&b0_relative_to_line);

    utils::approx_eq(
        f_sq,
        df_sq,
        1E-4,
        "f_sq and df_sq should match if bezier has endpoints at origin",
    );
    utils::approx_eq(
        c_sq,
        dc_sq,
        1E-4,
        "c_sq and dc_sq should match if bezier has endpoints at origin",
    );
}

#[test]
fn metrics_are_ordered() {
    let mut b0 = BezierBuilder::default();
    b0.add_point_at(0.0_f32, [0.0_f32]);
    b0.add_point_at(0.25, [-0.8]);
    b0.add_point_at(0.5, [0.2]);
    b0.add_point_at(0.75, [0.8]);
    b0.add_point_at(1.0, [0.0]);

    let mut b1 = BezierBuilder::default();
    b1.add_point_at(0.0_f32, [1.0_f32]);
    b1.add_point_at(0.25, [0.2]);
    b1.add_point_at(0.5, [0.8]);
    b1.add_point_at(0.75, [0.2]);
    b1.add_point_at(1.0, [1.0]);

    eprintln!("{b0:?}");
    eprintln!("{b1:?}");

    let b0 = b0.construct::<BezierND<_, 10, _>>().unwrap();
    let b1: BezierND<_, 10, _> = b1.construct().unwrap();

    eprintln!("{b0}");
    eprintln!("{b1}");
    check_metrics_of_beziers(&b0, &b1);

    let seed = "banana4";
    let mut rng = utils::make_random(seed);
    let distribution = rand::distr::Uniform::new(-10.0_f32, 10.0).unwrap();
    for _ in 0..10 {
        let b0: [[f32; 4]; 6] = utils::new_random_point_array(&mut rng, &distribution);
        let b1: [[f32; 4]; 6] = utils::new_random_point_array(&mut rng, &distribution);
        let b0 = BezierND::<_, 6, _>::new(&b0);
        let b1 = BezierND::<_, 6, _>::new(&b1);
        check_metrics_of_beziers(&b0, &b1);
    }
}
