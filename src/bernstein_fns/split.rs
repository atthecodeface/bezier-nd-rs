use crate::Num;
use geo_nd::vector;

/// Use de Casteljau's algorithm to split a Bernstein Bezier control
/// points set into two other Bernstein Bezier control point sets
/// at a given parameter t. The Bezier 'prior to t' is returned in `first`,
/// and `pts` is updated to be `beyond t`
///
/// The first Bezier returned has parameter t0 where 0<=t0<=1 maps to 0<=t*t0<=t
///
/// This destroys the provided points
pub fn into_two_at_de_cast<F: Num, const D: usize>(pts: &mut [[F; D]], t: F, first: &mut [[F; D]]) {
    // eprintln!("Split {pts:?} at {t}");
    // Beta[0][i] = pts[i]
    // for j = 1..=n
    //  for i = 0..=n-j
    //   Beta[j][i] = (1-t)*Beta[j-1][i] + t*Beta[j-1][i+1]
    let n = pts.len();
    assert!(first.len()>= n, "Splitting Bezier with {n} control points requires target Bezier 0 to have the same number of control points");
    first[0] = pts[0];
    let u = F::one() - t;
    for j in 1..n {
        // For j=1, p[0] = u*p[0]+t*p[1], p[1] = u*p[1]+t*p[2], ... p[n-2] = u*p[n-2]+t*p[n-1]
        // For j=n-1, p[0] = u*p[0] + t*p[1]
        for i in 0..(n - j) {
            pts[i] = vector::add(vector::scale(pts[i], u), &pts[i + 1], t);
        }
        first[j] = pts[0];
    }
    // eprintln!("Split into {b0:?} {b1:?}");
}

/// Use de Casteljau's algorithm to determine a Bernstein Bezier control
/// points set for the Bezier between t and 1.0
///
/// The Bezier is updated such that it has parameter t1 where 0<=t1<=1 maps to t<=t+(1-t)*t1<=1
pub fn bezier_to_de_cast<F: Num, const D: usize>(pts: &mut [[F; D]], t: F) {
    let n = pts.len();
    let u = F::one() - t;
    for j in 1..n {
        // For j=1, p[0] = u*p[0]+t*p[1], p[1] = u*p[1]+t*p[2], ... p[n-2] = u*p[n-2]+t*p[n-1]
        // For j=n-1, p[0] = u*p[0] + t*p[1] = u^(n-1)*p[0] +
        for i in (j..n).rev() {
            pts[i] = vector::add(vector::scale(pts[i - 1], u), &pts[i], t);
        }
    }
}

/// Use de Casteljau's algorithm to determine a Bernstein Bezier control
/// points set for the Bezier between t and 1.0
///
/// The Bezier is updated such that it has parameter t1 where 0<=t1<=1 maps to t<=t+(1-t)*t1<=1
pub fn bezier_from_de_cast<F: Num, const D: usize>(pts: &mut [[F; D]], t: F) {
    let n = pts.len();
    let u = F::one() - t;
    for j in 1..n {
        // For j=1, p[0] = u*p[0]+t*p[1], p[1] = u*p[1]+t*p[2], ... p[n-2] = u*p[n-2]+t*p[n-1]
        // For j=n-1, p[0] = u*p[0] + t*p[1] = u^(n-1)*p[0] +
        for i in 0..(n - j) {
            pts[i] = vector::add(vector::scale(pts[i], u), &pts[i + 1], t);
        }
    }
}
