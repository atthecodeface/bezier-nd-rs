use crate::polynomial::{PolyNewtonRaphson, Polynomial};
use crate::Num;
use geo_nd::vector;

pub fn min_tc<F: Num, T: Copy>((ta, a): (T, F), (tb, b): (T, F)) -> (T, F) {
    if a < b {
        (ta, a)
    } else {
        (tb, b)
    }
}

pub fn max_tc<F: Num, T: Copy>((ta, a): (T, F), (tb, b): (T, F)) -> (T, F) {
    if a > b {
        (ta, a)
    } else {
        (tb, b)
    }
}

pub fn opt_min_and_max_tc<F: Num, T: Copy>(
    give_min: bool,
    give_max: bool,
    ta_a: (T, F),
    tb_b: (T, F),
    opt_tc_c: Option<(T, F)>,
) -> (Option<(T, F)>, Option<(T, F)>) {
    if let Some(tc_c) = opt_tc_c {
        (
            give_min.then(|| min_tc(tc_c, min_tc(ta_a, tb_b))),
            give_max.then(|| max_tc(tc_c, max_tc(ta_a, tb_b))),
        )
    } else {
        (
            give_min.then(|| min_tc(ta_a, tb_b)),
            give_max.then(|| max_tc(ta_a, tb_b)),
        )
    }
}
/// Iterate in 'n' steps from t0 to t1 inclusive
pub fn float_iter_between<N: Num>(t0: N, t1: N, n: usize) -> impl Iterator<Item = N> {
    let step = {
        match n {
            0 => N::ZERO,
            1 => N::ZERO,
            n => (t1 - t0) / N::from_usize(n - 1).unwrap(),
        }
    };
    (0..n).map(move |i: usize| t0 + step * (N::from_usize(i).unwrap()))
}

/// Iterate in 'n' steps from 0 to 1 inclusive
pub fn float_iter<F: Num>(n: usize) -> impl Iterator<Item = F> {
    float_iter_between(F::ZERO, F::ONE, n)
}

/// Generate a linear mix of an equivalent straight Bezier from a parameter iterator
pub fn linear_iter<'a, F: Num, const D: usize, I: Iterator<Item = F> + 'a>(
    b: &'a [[F; D]],
    iter: I,
) -> impl Iterator<Item = [F; D]> + 'a {
    let p0 = b.first().unwrap();
    let p1 = b.last().unwrap();
    iter.map(|t| vector::mix(p0, p1, t))
}

/// Calculate the projection of the point onto the line, the length squared of a line,
/// and if the line is approximately zero length
pub fn relative_to_line<F: Num, const D: usize>(
    pt: &[F; D],
    l0: &[F; D],
    l1: &[F; D],
) -> (F, F, bool) {
    // Make everything relative to line.0
    //
    // so point is P and line is L (i.e. on the line = k*L for some 0 <= k <= 1)
    //
    // Note P = kL + N for some k and normal to the line N (N.L = 0)
    //
    // P.L = k(L.L) + N.L = k|L|^2
    let l1_m_l0 = vector::sub(*l1, l0, F::ONE);
    // len_l1_m_l0_sq = |L|^2
    let len_l1_m_l0_sq = vector::length_sq(&l1_m_l0);

    let pt_m_l0 = vector::sub(*pt, l0, F::ONE);

    if len_l1_m_l0_sq.is_unreliable_divisor() {
        return (F::ZERO, len_l1_m_l0_sq, false);
    }

    // pt_along_line = P . L = k|L|^2
    let pt_along_line = vector::dot(&pt_m_l0, &l1_m_l0);
    (pt_along_line, len_l1_m_l0_sq, true)
}

/// Calculate the distance squared to a line segment l0<>l1
pub fn distance_sq_to_line_segment<F: Num, const D: usize>(
    pt: &[F; D],
    l0: &[F; D],
    l1: &[F; D],
) -> F {
    let (t_times_len_sq, line_len_sq, valid) = crate::utils::relative_to_line(pt, l0, l1);
    if valid && t_times_len_sq >= F::ZERO && t_times_len_sq <= line_len_sq {
        let t = t_times_len_sq / line_len_sq;
        let u = F::ONE - t;
        vector::distance_sq(pt, &vector::sum_scaled(&[*l0, *l1], &[u, t]))
    } else if valid && t_times_len_sq >= line_len_sq {
        vector::distance_sq(pt, l1)
    } else {
        vector::distance_sq(pt, l0)
    }
}

/// Calculate Barycentric coordinates for a point projected onto a triangle,
/// or None if the triangle is degenerate (effectively a line)
///
/// The Barycentric coordinates for the projection of the point P onto the triangle
/// (at Q) are (k0, k1, k2) with k0+k1+k2=1, where Q=k0.P0 + k1.P1 + k2.P2
///
/// We can evaluate k1 and kw by making everything relative to P0 (denoted below with ', so Q' = Q-P0)
///
/// Q' = k1.P1' + k2.P2'; P'.P1' = Q'.P1' and P'.P2'=Q'.P2' (since Q is the projection onto O.P1'.P2')
///
/// Hence P'.P1' = Q'.P1' = k1.|P1'|^2 + k2.(P2'.P1'); P'.P2' = Q'.P2' = k2.|P2'|^2 + k1.(P2'.P1')
///
/// k1 =  ( |P2'|^2.(P'.P1') - (P2'.P1')(P'.P2') ) / (|P1'|^2.|P1'|^2 - (P1'.P2')^2)
/// k2 =  ( |P1'|^2.(P'.P2') - (P2'.P1')(P'.P1') ) / (|P1'|^2.|P1'|^2 - (P1'.P2')^2)
///
/// If |P1'|^2.|P1'|^2 ~= (P1'.P2')^2 then the triangle is degenerate
///
/// k0 = 1-k1-k2
///
/// A point is within or on the edge of the triangle if 0<=k0,k1,k2<=1
///
/// Only k1 and k2 are returned, and None if the triangle is degenerate
pub fn barycentric_coordinates<F: Num, const D: usize>(
    pt: &[F; D],
    p0: &[F; D],
    p1: &[F; D],
    p2: &[F; D],
) -> Option<(F, F)> {
    let p1 = vector::sub(*p1, p0, F::ONE);
    let p2 = vector::sub(*p2, p0, F::ONE);
    let pt = vector::sub(*pt, p0, F::ONE);
    let p_dot_p1 = vector::dot(&pt, &p1);
    let p_dot_p2 = vector::dot(&pt, &p2);
    let p1_dot_p2 = vector::dot(&p1, &p2);
    let p1_sq = vector::length_sq(&p1);
    let p2_sq = vector::length_sq(&p2);

    // Note that det is always positive
    let det = p1_sq * p2_sq - p1_dot_p2 * p1_dot_p2;
    if det.is_unreliable_divisor() {
        None
    } else {
        let k1 = (p2_sq * p_dot_p1 - p1_dot_p2 * p_dot_p2) / det;
        let k2 = (p1_sq * p_dot_p2 - p1_dot_p2 * p_dot_p1) / det;
        Some((k1, k2))
    }
}

/// Calculate how far from straight a cubic is using the distance squared from the
/// control point of a linear Bezier with the same endpoints
///
/// Effectively this is an estimate of the maximum squared size of Cubic - Linear Bezier
/// with the same endpoints
///
/// Cubic  = u^3.P0 + 3u^2t.P1      + 3ut^2.P2       + t^3.P3
/// Linear = u^3.P0 + u^2t.(2P0-P3) + ut^2.(P3-2P0)  + t^3.P3
///
/// Difference = u^2t.(3P1-2P0+P3) + 3ut^2.(3P2+2P0-P3)
/// Difference = ut.((1-t)(3P1-2P0+P3) + t(3P2+2P0-P3)) = ut.(u.A + t.B)
/// Difference sqaured = (ut)^2 . |u.A + t.B|^2
/// |uA + tB| (with u=1-t, 0<=t<=1) has a maximum value of max(|A|, |B|)
/// Hence difference squared <= (ut)^2 . max(A^2, B^2) <= max(A^2, B^2)
pub fn straightness_sq_of_cubic<F: Num, const D: usize>(cubic: &[[F; D]; 4]) -> F {
    let m_third = F::frac(-1, 3);
    let m_twothird = F::frac(-2, 3);
    let dv_0 = vector::sum_scaled(cubic, &[m_twothird, F::ONE, F::ZERO, m_third]);
    let dc2_0 = vector::length_sq(&dv_0);
    let dv_1 = vector::sum_scaled(cubic, &[m_third, F::ZERO, F::ONE, m_twothird]);
    let dc2_1 = vector::length_sq(&dv_1);
    if dc2_0 < dc2_1 {
        dc2_1
    } else {
        dc2_0
    }
}

pub fn find_real_roots_quad_num<F: Num>(poly: &[F]) -> (Option<F>, Option<F>) {
    assert!(
        poly.len() >= 3,
        "Root of a quadratic polynomial requires at least three coefficients (and [3..] should be zero)"
    );
    if poly[2].is_unreliable_divisor() {
        (crate::polynomial::find_real_roots_linear(poly), None)
    } else {
        let two_a = poly[2] + poly[2];
        let b = poly[1];
        let c = poly[0];
        let disc = b * b - F::of_i32(2) * two_a * c;
        if disc < F::zero() {
            (None, None)
        } else {
            let disc_sq = disc.sqrt_est();
            (Some((-b + disc_sq) / two_a), Some((-b - disc_sq) / two_a))
        }
    }
}

pub fn find_root_cubic_num<F: Num>(mut poly: [F; 4]) -> (Option<F>, Option<F>, Option<F>) {
    if poly[3].is_unreliable_divisor() {
        let (root_a, root_b) = find_real_roots_quad_num(&poly);
        (root_a, root_b, None)
    } else {
        poly[0] = poly[0] / poly[3];
        poly[1] = poly[1] / poly[3];
        poly[2] = poly[2] / poly[3];
        poly[3] = F::ONE;

        // If a NR fails at some point from t=-1 then probably at that t (a) the polynomial has
        // a small absolute value (b) d/dt of the polynomial is near zero (i.e. t is a root for that);
        //
        // If it fails again at a different point from t=2 then presuambly we have found the *two*
        // turning points of the cubic, which is unfortunate; but we can try the midpoint of the two
        // as that *must* be closer to the root with d/dt being larger
        //
        // If it fails again at the same point from t=2
        let mut t0 = F::ZERO;
        let mut p_t0 = F::of_i32(i32::MAX);
        for t in [0_i32, 1, 2, -1] {
            let (root_est, _root_dt) =
                poly.find_root_nr_with_err(F::of_i32(t), F::frac(1, 1_000_000), 20);
            let p_t = poly.calc(root_est).nabs();
            if p_t < p_t0 {
                t0 = root_est;
                p_t0 = p_t;
            }
        }
        if p_t0.nabs() > F::frac(1, 1_000) {
            // If it is not actually a root, i.e poly(t) > 0.001
            panic!("Failed to find root for polynomial {poly:?}");
            return (None, None, None);
        }

        eprintln!("Poly at t0 {t0} - {p_t0}");
        let (b, c) = {
            if t0.is_unreliable_divisor() {
                // Divide by x
                (poly[2], poly[1])
            } else {
                // Divide by (x-t0)
                let b = poly[2] + t0; // == (c - d_dt_sq_poly[1]) / t0;
                let c = poly[1] + t0 * b; // == d_dt_sq_poly[0] / t0;
                (b, c)
            }
        };
        let (opt_t1, opt_t2) = find_real_roots_quad_num(&[c, b, F::ONE]);
        (Some(t0), opt_t1, opt_t2)
    }
}

/// Estimate a square root given a Num
///
/// This uses Newton-Raphson on a scaled version of `sq`; it will always return a value
/// below or equal to the actual square root if 'min' is true, else a value larger or equal to
///
/// A value of N=5 will yield a square root accurate to one part in a 1E6.
pub fn sqrt_est<F: Num, const N: usize>(sq: F, min: bool) -> F {
    if sq < F::frac(1, 10_000) {
        if min {
            F::ZERO
        } else {
            sq
        }
    } else {
        if sq.is_unreliable_divisor() {
            return F::ZERO;
        }
        let scale = F::of_usize(4);
        let scale_r = F::frac(1, 4);
        let scale_sq = F::of_usize(16); // scale^2
        let scale_sq_r = F::frac(1, 16); // 1/scale^2
        let mut est = F::of_usize(2); // scale/2.0
        let half: F = F::frac(1, 2);

        let mut tsq = sq;
        while tsq < F::ONE {
            est *= scale_r;
            tsq *= scale_sq;
        }
        while tsq > scale_sq {
            est *= scale;
            tsq *= scale_sq_r;
        }
        for _ in 0..N {
            est = half * (est + sq / est);
        }
        let est2 = sq / est;
        if min == (est < est2) {
            est
        } else {
            est2
        }
    }
}

pub fn est_d_m_c_from_dsq_m_dcsq<F: Num>(d_sq: F, dc_sq: F) -> F {
    // Get min estimate for d_sq (i.e. d_est = d_sq.sqrt() - eps)
    let d_est = sqrt_est::<_, 4>(d_sq, true);
    // Get max estimate for dc_sq (i.e. dc_est = dc_sq.sqrt() + eps)
    let dc_est = sqrt_est::<_, 4>(dc_sq, false);

    if d_est > dc_est {
        F::ZERO
    } else {
        d_est - dc_est
    }
}

/// Estimate a cube root given a Num
///
/// This uses Newton-Raphson on a scaled version of `cb`
///
/// A value of N=3 will yield a square root accurate to one part in a 1E4.
///
/// A value of N=5 will yield a square root accurate to one part in a 1E7.
pub fn cbrt_est<F: Num, const N: usize>(cb: F) -> F {
    if cb.is_unreliable_divisor() {
        return F::ZERO;
    }
    let scale = F::of_usize(4);
    let scale_r = F::frac(1, 4);
    let scale_cb = F::of_usize(64); // scale^3
    let scale_cb_r = F::frac(1, 64); // 1/scale^3
    let mut est = F::of_usize(2); // scale/2.0

    let mut tcb = {
        if cb < F::ZERO {
            est = -est;
            -cb
        } else {
            cb
        }
    };
    while tcb < F::ONE {
        est *= scale_r;
        tcb *= scale_cb;
    }
    while tcb > scale_cb {
        est *= scale;
        tcb *= scale_cb_r;
    }
    let two_cb = cb + cb;
    for _ in 0..N {
        let est_cb = est * est * est;
        est *= (est_cb + two_cb) / (est_cb + est_cb + cb);
    }
    est
}

#[test]
fn test_sqrt() {
    fn test_x<F: Num + From<f32>, const N: usize>(x: f32, accuracy: f32) {
        let x: F = x.into();
        let x_sq = x * x;
        let est_min = sqrt_est::<_, N>(x_sq, true);
        let est_max = sqrt_est::<_, N>(x_sq, false);
        if x < 1E-6_f32.into() {
            assert!(est_min < 1E-6_f32.into());
            assert!(est_max < 1E-6_f32.into());
            return;
        }
        assert!(
            est_min <= x,
            "Minimum est of sqrt(x_sq) {est_min}<>{est_max} must be <= x {x}"
        );
        assert!(
            est_max >= x,
            "Maximum est of sqrt(x_sq) {est_min}<>{est_max} must be >= x {x}"
        );
        assert!(
            (x - est_min) / x < accuracy.into(),
            "Result {est_min}<>{est_max} should be within {accuracy} of x {x}"
        );
        assert!(
            (est_max - x) / x < accuracy.into(),
            "Result {est_min}<>{est_max} should be within {accuracy} of x {x}"
        );
    }
    for x in [0.0_f32, 1., 2., 3., 4., 5., 6., 7.] {
        for y in [1.0_f32, 1E2, 1E3, 1E4, 1E5, 1E-1, 1E-2] {
            test_x::<f32, 3>(x * y, 1E-3);
            test_x::<f32, 5>(x * y, 1E-7);
        }
    }
}

#[test]
fn test_cbrt() {
    fn test_x<F: Num + From<f32>, const N: usize>(x: f32, accuracy: f32) {
        let x: F = x.into();
        let x_cb = x * x * x;
        let est = cbrt_est::<_, N>(x_cb);
        if x > (-1E-6_f32).into() && x < 1E-6_f32.into() {
            assert!(est < 1E-2_f32.into());
            assert!(est > (-1E-2_f32).into());
            return;
        }
        assert!(
            (x - est) / x < accuracy.into(),
            "Result {est} should be within {accuracy} of x {x}"
        );
    }
    for x in [0.0_f32, 1., 2., 3., 4., 5., 6., 7.] {
        for y in [1.0_f32, 1E2, 1E3, 1E4, 1E5, 1E-1, 1E-2] {
            test_x::<f32, 3>(x * y, 1E-3);
            test_x::<f32, 3>(x * y, 1E-3);
            test_x::<f32, 4>(-x * y, 1E-6);
            test_x::<f32, 4>(-x * y, 1E-6);
        }
    }
}
