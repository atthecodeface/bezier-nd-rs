use bezier_nd::BezierEval;
use bezier_nd::{BezierBuilder, BezierND};
mod utils;

fn check_metrics_of_beziers<const N: usize, const D: usize>(
    b0: &BezierND<f32, N, D>,
    b1: &BezierND<f32, N, D>,
) {
    let mut max_d = 0.0;
    let mut max_t = 0.0;
    for t in utils::float_iter(0., 1., 101) {
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

    let l2_sq = b0.metric_l2_est(&b1, 10000);
    let dm_sq = b0.metric_dm_sq_est(&b1, 10000);
    let dc_sq = b0.metric_dc_sq(&b1);
    let df_sq = b0.metric_df_sq(&b1);
    eprintln!("b0-b1 l2 {l2_sq} dm {dm_sq} dc {dc_sq} df {df_sq}",);

    assert!(l2_sq <= dm_sq * 1.1, "L_SQ < dm_sq (with margin)");
    assert!(dm_sq <= dc_sq * 1.1, "dm_sq < dc_sq (with margin)");
    assert!(dc_sq <= df_sq * 1.01, "dc_sq < df_sq (with margin)");
    assert!(
        df_sq <= dc_sq * (N as f32) * 1.01,
        "df_sq < N*dc_sq (with margin)"
    );
}

#[test]
fn l2_is_not_halo() {
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

    let b0 = BezierND::<_, 10, _>::of_builder(b0).unwrap();
    let b1 = BezierND::<_, 10, _>::of_builder(b1).unwrap();

    eprintln!("{b0}");
    eprintln!("{b1}");
    check_metrics_of_beziers(&b0, &b1);

    let seed = "banana4";
    let mut rng = utils::make_random(seed);
    let distribution = rand::distr::Uniform::new(-10.0_f32, 10.0).unwrap();
    for _ in 0..1000 {
        let b0: [[f32; 4]; 6] = utils::new_random_point_array(&mut rng, &distribution);
        let b1: [[f32; 4]; 6] = utils::new_random_point_array(&mut rng, &distribution);
        let b0 = BezierND::<_, 6, _>::new(&b0);
        let b1 = BezierND::<_, 6, _>::new(&b1);
        check_metrics_of_beziers(&b0, &b1);
    }
}
