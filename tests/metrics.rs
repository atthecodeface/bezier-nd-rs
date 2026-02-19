use bezier_nd::BezierEval;
use bezier_nd::{metrics, BezierBuilder, BezierMetric, BezierND, BezierOps, Num};
mod utils;
use geo_nd::vector;
use rand::prelude::*;

fn check_metrics_of_beziers<
    F: Num + From<f32>,
    const D: usize,
    B: BezierEval<F, [F; D]> + BezierOps<F, [F; D]> + Clone + std::fmt::Debug,
>(
    b0: &B,
    b1: &B,
) {
    let mut max_d = F::ZERO;
    let mut max_t = F::ZERO;
    let mut total_d = F::ZERO;
    for t in utils::float_iter(1001) {
        let p0 = b0.point_at(t);
        let p1 = b1.point_at(t);
        let dp = vector::distance_sq(&p0, &p1);
        total_d += dp / 1001.0_f32.into();
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

    if b0.num_control_points() > 0 {
        assert!(
            b0.metric_from(Some(&[]), BezierMetric::SumControlSquared)
                .is_none(),
            "Metric to different sized bezier should be None"
        );
        assert!(
            b0.metric_from(Some(&[]), BezierMetric::MaxControlSquared)
                .is_none(),
            "Metric to different sized bezier should be None"
        );
        assert!(
            b0.metric_from(Some(&[]), BezierMetric::MaxDistanceSquared(10000))
                .is_none(),
            "Metric to different sized bezier should be None"
        );
        assert!(
            b0.metric_from(Some(&[]), BezierMetric::MaxDistanceSquared(10000))
                .is_none(),
            "Metric to different sized bezier should be None"
        );
    }

    let dl_sq = b0
        .metric_from(
            Some(b1.control_points()),
            BezierMetric::SumDistanceSquared(10000),
        )
        .unwrap();
    let dm_sq = b0
        .metric_from(
            Some(b1.control_points()),
            BezierMetric::MaxDistanceSquared(10000),
        )
        .unwrap();
    let dc_sq = b0
        .metric_from(Some(b1.control_points()), BezierMetric::MaxControlSquared)
        .unwrap();
    let df_sq = b0
        .metric_from(Some(b1.control_points()), BezierMetric::SumControlSquared)
        .unwrap();
    eprintln!("b0-b1 l2 {dl_sq} dm {dm_sq} dc {dc_sq} df {df_sq}",);

    let margin_1pc = 1.01_f32.into();
    let margin_5pc_l = 0.95_f32.into();
    let margin_5pc_h = 1.05_f32.into();
    let margin_10pc = 1.10_f32.into();
    assert!(
        dl_sq  >= total_d* margin_5pc_l && dl_sq  <= total_d*margin_5pc_h,
        "Sum distance squared metric {dl_sq} should be the same as the test calculated {total_d} within margin"
    );
    assert!(
        dm_sq >= max_d *margin_5pc_l && dm_sq  <= max_d* margin_5pc_h,
        "Max distance squared metric {dm_sq} should be the same as the test calculated {max_d} within margin"
    );
    assert!(dl_sq <= dm_sq * margin_10pc, "dl_sq < dm_sq (with margin)");
    assert!(dm_sq <= dc_sq * margin_10pc, "dm_sq < dc_sq (with margin)");
    assert!(dc_sq <= df_sq * margin_1pc, "dc_sq < df_sq (with margin)");
    assert!(
        df_sq <= dc_sq * (b0.num_control_points() as f32).into() * margin_1pc,
        "df_sq < N*dc_sq (with margin)"
    );

    let dl_sq = b0
        .metric_from(None, BezierMetric::SumDistanceSquared(10000))
        .unwrap();
    let dm_sq = b0
        .metric_from(None, BezierMetric::MaxDistanceSquared(10000))
        .unwrap();
    let df_sq = b0
        .metric_from(None, BezierMetric::SumControlSquared)
        .unwrap();
    let dc_sq = b0.dc_sq_from_line();
    assert!(dc_sq <= df_sq * margin_1pc, "dc_sq < df_sq (with margin)");
    assert!(
        df_sq <= dc_sq * (b0.num_control_points() as f32).into() * margin_1pc,
        "df_sq < N*dc_sq (with margin)"
    );

    if b0.num_control_points() < 2 {
        return;
    }
    let mut b0_relative_to_line = b0.clone();
    let (l0, l1) = b0.endpoints();
    b0_relative_to_line.map_pts(&|n, pt| {
        let t: F = ((n as f32) / (b0.degree() as f32)).into();
        eprintln!("{n} {t}");
        vector::sum_scaled(&[*pt, l0, l1], &[F::ONE, t - F::ONE, -t])
    });
    eprintln!("B0 {b0:?} relative {b0_relative_to_line:?}");

    let l_sq = b0_relative_to_line
        .metric_from(None, BezierMetric::SumDistanceSquared(10000))
        .unwrap();
    let m_sq = b0_relative_to_line
        .metric_from(None, BezierMetric::MaxDistanceSquared(10000))
        .unwrap();
    let f_sq = b0_relative_to_line
        .metric_from(None, BezierMetric::SumControlSquared)
        .unwrap();
    let c_sq = b0_relative_to_line
        .metric_from(None, BezierMetric::MaxControlSquared)
        .unwrap();

    if dl_sq.is_unreliable_divisor() {
        assert!(l_sq.is_unreliable_divisor());
        assert!(m_sq.is_unreliable_divisor());
    } else {
        utils::approx_eq(
            l_sq / dl_sq,
            F::ONE,
            1.001,
            "l_sq and dl_sq should match if bezier has endpoints at origin",
        );

        utils::approx_eq(
            m_sq / dm_sq,
            F::ONE,
            1.001,
            "m_sq and dm_sq should match if bezier has endpoints at origin",
        );
    }

    if df_sq.is_unreliable_divisor() {
        assert!(f_sq.is_unreliable_divisor());
        assert!(c_sq.is_unreliable_divisor());
    } else {
        utils::approx_eq(
            f_sq / df_sq,
            F::ONE,
            1.001,
            "f_sq and df_sq should match if bezier has endpoints at origin",
        );

        utils::approx_eq(
            c_sq / dc_sq,
            F::ONE,
            1.001,
            "c_sq and dc_sq should match if bezier has endpoints at origin",
        );
    }
}

fn check_metrics_of_beziers_test<
    R: Rng,
    D: Distribution<f32>,
    F: Num + From<f32>,
    const N: usize,
    B: BezierEval<F, [F; N]> + BezierOps<F, [F; N]> + Clone + std::fmt::Debug,
>(
    rng: &mut R,
    distribution: &D,
    b0: &mut B,
    b1: &mut B,
) {
    utils::set_random_bezier(rng, distribution, b0);
    utils::set_random_bezier(rng, distribution, b1);
    check_metrics_of_beziers::<F, N, B>(&b0, &b1);
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
        let mut b0 = BezierND::<_, 6, _>::new(&[[0.0_f32; 2]; 3]);
        let mut b1 = BezierND::<_, 6, _>::new(&[[0.0_f32; 2]; 3]);
        check_metrics_of_beziers_test(&mut rng, &distribution, &mut b0, &mut b1);

        let mut b0 = BezierND::<_, 6, _>::new(&[[0.0_f32; 2]; 6]);
        let mut b1 = BezierND::<_, 6, _>::new(&[[0.0_f32; 2]; 6]);
        check_metrics_of_beziers_test(&mut rng, &distribution, &mut b0, &mut b1);

        let mut b0 = [[0.0_f32; 3]; 2];
        let mut b1 = [[0.0_f32; 3]; 2];
        check_metrics_of_beziers_test(&mut rng, &distribution, &mut b0, &mut b1);

        let mut b0 = [[0.0_f32; 3]; 3];
        let mut b1 = [[0.0_f32; 3]; 3];
        check_metrics_of_beziers_test(&mut rng, &distribution, &mut b0, &mut b1);

        let mut b0 = [[0.0_f32; 3]; 4];
        let mut b1 = [[0.0_f32; 3]; 4];
        check_metrics_of_beziers_test(&mut rng, &distribution, &mut b0, &mut b1);

        let mut b0 = vec![[0.0_f32; 2]; 2];
        let mut b1 = vec![[0.0_f32; 2]; 2];
        check_metrics_of_beziers_test(&mut rng, &distribution, &mut b0, &mut b1);

        let mut b0 = vec![[0.0_f32; 2]; 3];
        let mut b1 = vec![[0.0_f32; 2]; 3];
        check_metrics_of_beziers_test(&mut rng, &distribution, &mut b0, &mut b1);

        let mut b0 = vec![[0.0_f32; 2]; 4];
        let mut b1 = vec![[0.0_f32; 2]; 4];
        check_metrics_of_beziers_test(&mut rng, &distribution, &mut b0, &mut b1);

        let mut b0 = vec![[0.0_f32; 2]; 7];
        let mut b1 = vec![[0.0_f32; 2]; 7];
        check_metrics_of_beziers_test(&mut rng, &distribution, &mut b0, &mut b1);

        let mut b0 = vec![[0.0_f32; 2]; 0];
        let mut b1 = vec![[0.0_f32; 2]; 0];
        check_metrics_of_beziers_test(&mut rng, &distribution, &mut b0, &mut b1);

        let mut b0 = vec![[0.0_f32; 2]; 1];
        let mut b1 = vec![[0.0_f32; 2]; 1];
        check_metrics_of_beziers_test(&mut rng, &distribution, &mut b0, &mut b1);
    }
}
