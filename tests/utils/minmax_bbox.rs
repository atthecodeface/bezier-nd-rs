use bezier_nd::{bernstein_fns, BezierEval, BezierSplit};

use bezier_nd::Num;
use geo_nd::{matrix, vector};

pub fn assert_min_max_coords<
    F: Num + From<f32>,
    const D: usize,
    B: BezierEval<F, [F; D]> + std::fmt::Debug,
>(
    bezier: &B,
) {
    eprintln!("Min/max of {bezier:?}");
    let pts = super::BezierPtSet::of_point_at(bezier, 1000);
    let bbox = pts.bbox();
    eprintln!("Bbox of BezierPtSet of point_at(0<=t<=1) = {bbox:?}");
    for d in 0..D {
        let (omin, omax) = bezier.t_coords_at_min_max(d, false, false);
        assert!(omin.is_none());
        assert!(omax.is_none());
        let (omin, omax) = bezier.t_coords_at_min_max(d, true, false);
        assert!(omin.is_some());
        assert!(omax.is_none());
        let (omin, omax) = bezier.t_coords_at_min_max(d, false, true);
        assert!(omin.is_none());
        assert!(omax.is_some());
        let (omin, omax) = bezier.t_coords_at_min_max(d, true, true);
        let (t_min_d, f_min_d) = omin.unwrap();
        let (t_max_d, f_max_d) = omax.unwrap();
        eprintln!(
            "BezierMinMax for coord {d} has min {f_min_d} <> max {f_max_d} at {t_min_d} {t_max_d}"
        );
        super::approx_eq(f_min_d, bezier.point_at(t_min_d)[d], 1E-6, "Pt at t_min_d");
        super::approx_eq(f_max_d, bezier.point_at(t_max_d)[d], 1E-6, "Pt at t_min_d");
        super::approx_eq(f_min_d, bbox.0[d], 1E-4, "Bbox min coord d");
        super::approx_eq(f_max_d, bbox.1[d], 1E-4, "Bbox min coord d");
    }
}
