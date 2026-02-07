use crate::Num;
use geo_nd::vector;

pub fn max<F: Num>(a: F, b: F) -> F {
    if a < b {
        b
    } else {
        a
    }
}

pub fn min_or_max<F: Num>(use_max: bool, ta: F, a: F, tb: F, b: F) -> (F, F) {
    if use_max == (a < b) {
        (tb, b)
    } else {
        (ta, a)
    }
}

/// Iterate in 'n' steps from t0 to t1 inclusive
pub fn float_iter_between<N: Num>(t0: N, t1: N, n: usize) -> impl Iterator<Item = N> {
    assert!(
        n >= 2,
        "Float iterator must have at least two steps for begin and end"
    );
    let r = t1 - t0;
    assert!(r >= N::zero(), "Float range must be positive");
    (0..n).map(move |i: usize| t0 + r * N::from_usize(i).unwrap() / N::from_usize(n - 1).unwrap())
}

/// Iterate in 'n' steps from 0 to 1 inclusive
pub fn float_iter<F: Num>(n: usize) -> impl Iterator<Item = F> {
    float_iter_between(F::ZERO, F::ONE, n)
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
pub fn straightness_sq_of_cubic<F: crate::Num, const D: usize>(cubic: &[[F; D]; 4]) -> F {
    let m_third = (-0.33333333_f32).into();
    let m_twothird = (-0.666_666_7_f32).into();
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
