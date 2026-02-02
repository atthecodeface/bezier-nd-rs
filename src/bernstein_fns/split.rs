use crate::Num;
use geo_nd::vector;

/// Use de Casteljau's algorithm to split a Bernstein Bezier control
/// points set into two other Bernstein Bezier control point sets
/// at a given parameter t
///
/// The first Bezier returned has parameter t0 where 0<=t0<=1 maps to 0<=t*t0<=t
///
/// The second Bezier returned has parameter t1 where 0<=t1<=1 maps to t<=t+(1-t)*t1<=1
///
/// This destroys the provided points
pub fn into_two_at_de_cast<F: Num, const D: usize>(
    pts: &mut [[F; D]],
    t: F,
    b0: &mut [[F; D]],
    b1: &mut [[F; D]],
) {
    // eprintln!("Split {pts:?} at {t}");
    // Beta[0][i] = pts[i]
    // for j = 1..=n
    //  for i = 0..=n-j
    //   Beta[j][i] = (1-t)*Beta[j-1][i] + t*Beta[j-1][i+1]
    let n = pts.len();
    assert!(b0.len()>= n, "Splitting Bezier with {n} control points requires target Bezier 0 to have the same number of control points");
    assert!(b1.len()>= n, "Splitting Bezier with {n} control points requires target Bezier 1 to have the same number of control points");
    b0[0] = pts[0];
    b1[n - 1] = pts[n - 1];
    let u = F::one() - t;
    for j in 1..n {
        // For j=1, p[0] = u*p[0]+t*p[1], p[1] = u*p[1]+t*p[2], ... p[n-2] = u*p[n-2]+t*p[n-1]
        // For j=n-1, p[0] = u*p[0] + t*p[1]
        for i in 0..(n - j) {
            pts[i] = vector::add(vector::scale(pts[i], u), &pts[i + 1], t);
        }
        b0[j] = pts[0];
        b1[n - 1 - j] = pts[n - 1 - j];
    }
    // eprintln!("Split into {b0:?} {b1:?}");
}
