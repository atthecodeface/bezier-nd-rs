//a Imports
use crate::Float;
use geo_nd::vector;

use crate::BezierEval;

//fi lambda_of_k_d
fn lambda_of_k_d<F: Float>(k: F, d: F) -> F {
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
    let a0: F = 3.160_323_6e-2_f32.into();
    let a1: F = 2.795_054_2_f32.into();
    let a2: F = (-1.148_674_3e1_f32).into();
    let a3: F = 2.897_536_8e1_f32.into();
    let a4: F = (-3.284_522_2e1_f32).into();
    let a5: F = 1.377_942_9e1_f32.into();
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
pub fn arc<F: Float, const D: usize>(
    angle: F,
    radius: F,
    center: &[F; D],
    unit: &[F; D],
    normal: &[F; D],
    rotate: F,
) -> [[F; D]; 4] {
    let two = (2.0_f32).into();
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
pub fn of_round_corner<F: Float, const D: usize>(
    corner: &[F; D],
    v0: &[F; D],
    v1: &[F; D],
    radius: F,
) -> [[F; D]; 4] {
    let nearly_one = (0.999_999_f32).into();
    let one = F::one();
    let two: F = (2.0_f32).into();
    let v0 = vector::normalize(*v0);
    let v1 = vector::normalize(*v1);
    let cos_alpha = vector::dot(&v0, &v1);
    if cos_alpha >= nearly_one || cos_alpha < -nearly_one {
        // v0 and v1 point in the same direction
        let mut p0 = [F::zero(); D];
        let mut p1 = [F::zero(); D];
        for i in 0..D {
            p0[i] = corner[i] - radius * v0[i];
            p1[i] = corner[i] - radius * v1[i];
        }
        let c0 = vector::sum_scaled(
            &[p0, *corner],
            &[0.3333333_f32.into(), 0.66666667_f32.into()],
        );
        let c1 = vector::sum_scaled(
            &[p1, *corner],
            &[0.3333333_f32.into(), 0.66666667_f32.into()],
        );
        [p0, c0, c1, p1]
    } else {
        let r2 = radius * radius;
        let d2 = two * r2 / (one - cos_alpha);
        let k2 = d2 - r2;
        let d = d2.sqrt();
        let k = k2.sqrt();

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
/// what is the center of the circle
/// given point p0 and unit tangent t0
/// and point p1 and unit tangent t1
///
/// ```text
/// |p0-c|^2 = |p1-c|^2 = r^2
/// (p0-c) . t0 = 0
/// (p1-c) . t1 = 0
/// ```
///
/// Consider c = k0.t0 + k1.t1
///
/// (given t0.t0 == 1 and t1.t1==1)
///
/// ```text
/// (p0-c) . t0 = (p0 - k0.t0 - k1.t1).t0 = 0
///       p0.t0 = k0 + k1(t1.t0)
/// ```
/// similarly
/// ```text
///       p1.t1 = k1 + k0(t1.t0)
/// ```
///
/// hence
/// ```text
///  (t1.t0) * (p1.t1)         = k0.(t1.t0)^2 + k1(t1.t0)
///  p0.t0 - (t1.t0) * (p1.t1) = k0 ( 1 - (t1.t0)^2)
///  k0 = (p0.t0 - p1.t1 * t1.t0) / ( 1 - (t1.t0)^2)
///  k1 = (p1.t1 - p0.t0 * t1.t0) / ( 1 - (t1.t0)^2)
/// ```
pub fn center_radius_of_bezier_arc<F: Float, const D: usize, B: BezierEval<F, [F; D]>>(
    bezier: &B,
) -> ([F; D], F) {
    let zero = F::zero();
    let one = F::one();
    let p0 = bezier.point_at(zero);
    let p1 = bezier.point_at(one);
    let (_sc0, t0) = bezier.derivative_at(zero);
    let (_sc1, t1) = bezier.derivative_at(one);
    let t0 = vector::normalize(t0);
    let t1 = vector::normalize(t1);
    let t1_d_t0 = vector::dot(&t1, &t0);
    let p0_d_t0 = vector::dot(&p0, &t0);
    let p1_d_t1 = vector::dot(&p1, &t1);
    let k0 = (p0_d_t0 - p1_d_t1 * t1_d_t0) / (one - t1_d_t0 * t1_d_t0);
    let k1 = (p1_d_t1 - p0_d_t0 * t1_d_t0) / (one - t1_d_t0 * t1_d_t0);

    let mut c = [F::zero(); D];
    for i in 0..D {
        c[i] = t0[i] * k0 + t1[i] * k1;
    }

    let r = (vector::distance(&c, &p0) + vector::distance(&c, &p1)) / (2.0_f32).into();
    (c, r)
}
