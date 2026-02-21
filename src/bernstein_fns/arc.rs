use crate::{BezierEval, Num};
use geo_nd::vector;

fn lambda_of_k_d<F: Num>(k: F, d: F) -> F {
    // There are numerous versions of calculating
    // the lambda for the arc from the angle of the arc
    //
    // For a 90 degree arc the *best* values is 0.2652165 apparently
    //
    // One equation that provides this  is
    //   lambda = four_thirds * radius * (d/k - one);
    //
    // Another is:
    //   theta = (k/d).asin() / F::int(4);
    //   lambda = four_thirds * theta.tan();
    //
    // This table is captures the values for this second
    //
    // Actually attempting a better approximation leads to the following for (k/d)^2 -> lambda
    //
    // 0.0011099165 0.009397572
    // 0.004424813 0.04415609
    // 0.008196682 0.06044403
    // 0.012195113 0.07385845
    // 0.019999988 0.09472578
    // 0.038461603 0.13201918
    // 0.100000046 0.21637033
    // 0.100000046 0.21637033
    // 0.137931 0.2567711
    // 0.20000009 0.314736
    // 0.3076923 0.40363038
    // 0.5 0.5519717
    // 0.6923078 0.71254206
    // 0.8000001 0.822074
    // 0.862069 0.89976513
    // 0.8999999 0.9571549
    // 0.9615385 1.0864261
    // 0.98 1.1479391
    // 0.99180335 1.2072284
    // 0.9955752 1.2359663
    // 0.9988901 1.2764238
    //
    // With a quintic polynomial of coeffs (x^0 + x^1 +... + x^5):
    // 3.1603235091816735e-002
    // 2.7950542994656820e+000
    // -1.1486743224812313e+001
    // 2.8975368657401102e+001
    // -3.2845222512637491e+001
    // 1.3779429574112177e+001
    //
    // Or for r^2/d^2 -> lambda
    // 0.9988901 0.009397572
    // 0.9955752 0.04415609
    // 0.99180335 0.06044403
    // 0.9878049 0.07385845
    // 0.98 0.09472578
    // 0.9615384 0.13201918
    // 0.9 0.21637033
    // 0.9 0.21637033
    // 0.862069 0.2567711
    // 0.79999995 0.314736
    // 0.6923077 0.40363038
    // 0.5 0.5519717
    // 0.3076923 0.71254206
    // 0.20000002 0.822074
    // 0.13793105 0.89976513
    // 0.10000005 0.9571549
    // 0.038461536 1.0864261
    // 0.02000001 1.1479391
    // 0.008196682 1.2072284
    // 0.0044247806 1.2359663
    // 0.0011098981 1.2764238
    //
    // 1.2494900596889080e+000
    // -4.2639321404424191e+000
    // 1.6162330324360198e+001
    // -3.5388797293367219e+001
    // 3.6051953254575963e+001
    // -1.3779440945693199e+001

    let k_d = k / d;
    // let four_thirds  = F::frac(4,3);
    // let  theta = (k/d).asin() / F::int(4);
    // let  lambda = four_thirds * theta.tan();
    // lambda
    let k_d = k_d * k_d;
    let a0 = F::frac(31_603_236, 1_000_000_000);
    let a1 = F::frac(2_795_0542, 10_000_000);
    let a2 = F::frac(-11_486_743, 1_000_000);
    let a3 = F::frac(28_975_368, 1_000_000);
    let a4 = F::frac(-32_845_222, 1_000_000);
    let a5 = F::frac(13_779_429, 1_000_000);
    a0 + a1 * k_d
        + a2 * k_d * k_d
        + a3 * k_d * k_d * k_d
        + a4 * k_d * k_d * k_d * k_d
        + a5 * k_d * k_d * k_d * k_d * k_d
}

/// Create a Cubic Bezier that approximates closely a circular arc
///
/// The arc has a center C, a radius R, and is of an angle (should be <= PI/2).
///
/// The arc sweeps through points a distance R from C, in a circle
/// using a pair of the planar unit vectors in the vector space for the
/// points.
///
/// The arc will be between an angle A1 and A2, where A2-A1 == angle, and A1==rotate
///
pub fn arc<F: Num + num_traits::Float, const D: usize>(
    angle: F,
    radius: F,
    center: &[F; D],
    unit: &[F; D],
    normal: &[F; D],
    rotate: F,
) -> [[F; D]; 4] {
    let two = F::of_usize(2);
    let half_angle = angle / two;
    let s = half_angle.sin();
    let lambda = radius * lambda_of_k_d(s, F::one());

    let d0a = rotate;
    let (d0s, d0c) = d0a.sin_cos();
    let d1a = rotate + angle;
    let (d1s, d1c) = d1a.sin_cos();

    let mut p0 = [F::zero(); D];
    let mut p1 = [F::zero(); D];
    let mut c0 = [F::zero(); D];
    let mut c1 = [F::zero(); D];
    for i in 0..D {
        p0[i] = center[i] + unit[i] * (d0c * radius) + normal[i] * (d0s * radius);
        p1[i] = center[i] + unit[i] * (d1c * radius) + normal[i] * (d1s * radius);

        c0[i] = p0[i] - unit[i] * (d0s * lambda) + normal[i] * (d0c * lambda);
        c1[i] = p1[i] + unit[i] * (d1s * lambda) - normal[i] * (d1c * lambda);
    }
    [p0, c0, c1, p1]
}

/// Create a Cubic Bezier that is a circular arc focused on the corner point,
/// with v0 and v1 are vectors IN to the point (P)
///
/// As it is a circular arc we have a kite P, P+k.v0, C, P+k.v1, where
///
/// ```text
/// |P+k.v0 - C| = |P+k.v1 - C| = r; |P-C| = d (i.e. side lengths are r, r, k, k)
/// ```
///
/// with two corners being right-angles. (and d is the length of
/// the kite diagonal opposite these right-angles).
///
/// The kite is formed from two d, r, k right-angled triangles; it
/// has two other angles, alpha and 2*theta, (alpha = angle
/// between v0 and v1). Hence alpha = 180 - 2*theta, theta = 90-(alpha/2)
///
/// Hence d^2 = r^2 + k^2; r/d = cos(theta), k/d=sin(theta)
///
/// We know cos(alpha) = v0.v1 (assuming unit vectors).
///
/// ```text
/// cos(alpha) = cos(180-2*theta)
///            = -cos(2*theta)
///            = -(2cos^2(theta) - 1)
///            = 1 - 2cos^2(theta)
///
/// cos^2(theta) = (1 - cos(alpha)) / 2 = r^2/d^2
///
/// sin^2(theta) = (1 + cos(alpha)) / 2
///
/// => d^2 = 2*r^2  / (1 - cos(alpha))
/// ```
///
/// Hence also k^2, and hence d and k.
///
/// Then we require an arc given the angle of the arc is 2*theta
pub fn of_round_corner<F: Num, const D: usize>(
    corner: &[F; D],
    v0: &[F; D],
    v1: &[F; D],
    radius: F,
) -> [[F; D]; 4] {
    let nearly_one = F::frac(999_999, 1_000_000);
    let one = F::ONE;
    let two = F::of_usize(2);
    let v0_l = vector::length_sq(v0).sqrt_est();
    let v1_l = vector::length_sq(v1).sqrt_est();
    let v0 = if v0_l.is_unreliable_divisor() {
        *v0
    } else {
        vector::reduce(*v0, v0_l)
    };
    let v1 = if v1_l.is_unreliable_divisor() {
        *v1
    } else {
        vector::reduce(*v1, v1_l)
    };
    let cos_alpha = vector::dot(&v0, &v1);
    if cos_alpha >= nearly_one || cos_alpha < -nearly_one {
        // v0 and v1 point in the same direction
        let mut p0 = [F::zero(); D];
        let mut p1 = [F::zero(); D];
        for i in 0..D {
            p0[i] = corner[i] - radius * v0[i];
            p1[i] = corner[i] - radius * v1[i];
        }
        let c0 = vector::sum_scaled(&[p0, *corner], &[F::frac(1, 3), F::frac(2, 3)]);
        let c1 = vector::sum_scaled(&[p1, *corner], &[F::frac(1, 3), F::frac(2, 3)]);
        [p0, c0, c1, p1]
    } else {
        let r2 = radius * radius;
        let d2 = two * r2 / (one - cos_alpha);
        let k2 = d2 - r2;
        let d = d2.sqrt_est();
        let k = k2.sqrt_est();

        let lambda = radius * lambda_of_k_d(k, d);

        let mut p0 = [F::zero(); D];
        let mut p1 = [F::zero(); D];
        let mut c0 = [F::zero(); D];
        let mut c1 = [F::zero(); D];
        for i in 0..D {
            p0[i] = corner[i] - k * v0[i];
            p1[i] = corner[i] - k * v1[i];
            c0[i] = p0[i] + lambda * v0[i];
            c1[i] = p1[i] + lambda * v1[i];
        }
        [p0, c0, c1, p1]
    }
}

/// Find the center and radius of a Bezier if it is assumed to
/// be a circular arc
///
/// To calculate the center of the circle given point p0 and tangent t0 and
/// point p1 and tangent t1:
///
/// ```text
/// |p0-c|^2 = |p1-c|^2 = r^2
/// (p0-c) . t0 = 0 (t0 is tangent to circle, p0-c is radius)
/// (p1-c) . t1 = 0 (t1 is tangent to circle, p1-c is radius)
/// ```
///
/// Consider c = k0.t0 + k1.t1 for some k0 and k1
///
/// ```text
///    (p0-c) . t0 = 0
///         0 = (p0 - k0.t0 - k1.t1).t0
///     p0.t0 = k0(t0.t0) + k1(t1.t0). [1]
/// &   p1.t1 = k0(t1.t0) + k1(t1.t1)  [2] (similarly)
/// => [2] * (t1.t0)
///  (t1.t0) * (p1.t1)         = k0(t1.t0)^2 + k1(t1.t1)(t1.t0) [3]
/// => [1]*(t1.t1) - [3]
///  (p0.t0)(t1.t1) - (t1.t0)(p1.t1) = k0(t0.t0)(t1.t1) + k1(t1.t0)(t1.t1) - k0.(t1.t0)^2 - k1(t1.t0)(t1.t1)
///  p0.t0*(t1.t1) - (t1.t0) * (p1.t1) = k0 ((t1.t1)(t0.t0) - (t1.t0)^2)
/// => k0 = (p0.t0 * |t1|^2 - p1.t1 * t1.t0) / ( |t0|^2.|t1|^2 - (t1.t0)^2)
///  & k1 = (p1.t1 * |t0|^2 - p0.t0 * t1.t0) / ( |t0|^2.|t1|^2 - (t1.t0)^2)
/// ```
pub fn center_radius_of_bezier_arc<
    F: Num,
    const D: usize,
    B: BezierEval<F, [F; D]> + std::fmt::Debug,
>(
    bezier: &B,
) -> ([F; D], F) {
    let p0 = bezier.endpoints().0;
    let p1 = bezier.endpoints().1;
    let (_sc0, t0) = bezier.derivative_at(F::ZERO);
    let (_sc1, t1) = bezier.derivative_at(F::ONE);
    let t0_l2 = vector::length_sq(&t0);
    let t1_l2 = vector::length_sq(&t1);
    let t1_d_t0 = vector::dot(&t1, &t0);
    let p0_d_t0 = vector::dot(&p0, &t0);
    let p1_d_t1 = vector::dot(&p1, &t1);

    let det = t1_l2 * t0_l2 - t1_d_t0 * t1_d_t0;
    let k0 = (p0_d_t0 * t1_l2 - p1_d_t1 * t1_d_t0) / det;
    let k1 = (p1_d_t1 * t0_l2 - p0_d_t0 * t1_d_t0) / det;

    let c = vector::sum_scaled(&[t0, t1], &[k0, k1]);
    let r_sq_est_0 = vector::distance_sq(&c, &p0);
    let r_sq_est_1 = vector::distance_sq(&c, &p1);
    let r_sq_est = (r_sq_est_0 + r_sq_est_1) * F::frac(1, 2);
    let r = r_sq_est.sqrt_est();
    (c, r)
}
