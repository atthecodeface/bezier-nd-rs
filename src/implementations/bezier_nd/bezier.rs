use crate::{BezierSplit, Num};
use geo_nd::vector;

/// This type stores a Bernstein Bezier of up to N control points each of dimension D
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct Bezier<F, const N: usize, const D: usize>
where
    F: Num,
{
    /// Degree, which is one less than number of valid control points
    ///
    /// i.e. for a cubic this is three
    ///
    /// A Bezier of degree 0 is a fixed point
    pub(super) degree: usize,
    /// Control points - 0..=degree are valid
    pub(super) pts: [[F; D]; N],
}

//ti Default for Bezier
impl<F, const N: usize, const D: usize> std::default::Default for Bezier<F, N, D>
where
    F: Num,
{
    fn default() -> Self {
        let pts = [[F::ZERO; D]; N];
        Self { degree: 0, pts }
    }
}

//ti Display for Bezier
impl<F, const N: usize, const D: usize> std::fmt::Display for Bezier<F, N, D>
where
    F: Num,
{
    //mp fmt - format a `Bezier` for display
    /// Display the `Bezier' as sets of points
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "[")?;
        let mut needs_comma = false;
        for p in self.pts.iter().take(self.degree + 1) {
            if needs_comma {
                write!(f, ", ")?;
            }
            needs_comma = true;
            vector::fmt(f, p)?;
        }
        write!(f, "]")
    }

    //zz All done
}

//ip Bezier iterators
impl<F, const N: usize, const D: usize> Bezier<F, N, D>
where
    F: Num,
{
    /// Find the maximum degree this type can encode
    pub fn max_degree() -> usize {
        N - 1
    }

    //mp Reduce-and-split iterator
    /// Apply a (new_degree+1) by (degree+1) matrix to the points to generate a new Bezier
    /// of a new degree
    pub fn reduce_and_split_iter<'a>(
        &'a self,
        reduce_matrix: &'a [F],
        elev_reduce_matrix: &'a [F],
        reduce_degree: usize,
        max_dc_sq: F,
    ) -> BezierReduceIter<'a, F, N, D> {
        BezierReduceIter {
            reduce_matrix,
            elev_reduce_matrix,
            reduce_degree,
            max_dc_sq,
            stack: vec![],
        }
    }
}

//tp BezierReduceIter
/// A type that provided an iterator implementaion for splitting a Bezier
/// into reduced Beziers given a maximum 'dc' metric
#[derive(Clone, PartialEq, Debug)]
pub struct BezierReduceIter<'a, F, const N: usize, const D: usize>
where
    F: Num,
{
    reduce_matrix: &'a [F],
    elev_reduce_matrix: &'a [F],
    reduce_degree: usize,
    max_dc_sq: F,
    stack: Vec<(usize, Bezier<F, N, D>)>,
}

//ip Iterator for BezierReduceIter
impl<'a, F, const N: usize, const D: usize> std::iter::Iterator for BezierReduceIter<'a, F, N, D>
where
    F: Num,
{
    type Item = (usize, Bezier<F, N, D>);
    fn next(&mut self) -> Option<(usize, Bezier<F, N, D>)> {
        while let Some((n, b)) = self.stack.pop() {
            let dc2 = b.dc2_of_ele_red(self.elev_reduce_matrix);
            if dc2 > self.max_dc_sq {
                let (b0, b1) = b.split();
                self.stack.push((n + 1, b1));
                self.stack.push((n + 1, b0));
            } else {
                return Some((n, b.apply_matrix(self.reduce_matrix, self.reduce_degree)));
            }
        }
        None
    }
}
