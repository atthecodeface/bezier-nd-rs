use crate::BezierSplit;
use crate::Num;

/// An iterator with Item = (F, F, Bezier), allowing for splitting dynamically during
/// iteration
#[derive(Clone)]
pub enum BezierSplitTIter<B: Clone, F: Num> {
    /// An iterator that *splits* the Bezier into halves recursively
    Split {
        /// A stack of future beziers to examine
        /// The top of the stack is p0->p1; below that is p1->p2, etc
        /// These beziers may need to be split to achieve some critertion
        stack: Vec<(F, F, B)>,
    },
    /// An iterator over a Bezier using a range of `t` values
    TRange {
        /// The Bezier that is being iterated over
        bezier: B,
        /// Initial value of `t`
        t0: F,
        /// The amount to increase `t` by for each iteration
        t_scale: F,
        /// The step of iteration (from 0 to num_steps)
        step: usize,
        /// The number of iteration steps
        num_steps: usize,
    },
}

//pi BezierLineIter
impl<B, F> BezierSplitTIter<B, F>
where
    B: Clone,
    F: Num,
{
    /// Create a new Bezier line iterator for a given Bezier and
    /// straightness
    ///
    /// This clones the Bezier.
    pub fn split_bezier(bezier: &B) -> Self {
        let stack = vec![(F::ZERO, F::ONE, bezier.clone())];
        Self::Split { stack }
    }

    /// Create a new Bezier line iterator for a given Bezier and
    /// straightness
    ///
    /// This clones the Bezier.
    pub fn uniform(bezier: &B, t0: F, t1: F, num_steps: usize) -> Self {
        let bezier = bezier.clone();
        let t_scale = {
            if num_steps <= 1 {
                F::ONE
            } else {
                (t1 - t0) / (F::from_usize(num_steps).unwrap() - F::ONE)
            }
        };
        let step = 0;
        Self::TRange {
            bezier,
            t0,
            t_scale,
            step,
            num_steps,
        }
    }

    /// Splits and adds both Bezier to the iterator; the first half
    /// of the split Bezier will be returned next by the iterator
    pub fn add_split(&mut self, (t0, t1, bezier): &(F, F, B))
    where
        B: BezierSplit<F>,
    {
        match self {
            Self::Split { ref mut stack } => {
                let tm = (*t0 + *t1) * F::frac(1, 2);
                let (b0, b1) = bezier.split();
                stack.push((tm, *t1, b1));
                stack.push((*t0, tm, b0));
            }
            _ => {
                panic!("Cannot split a TRange");
            }
        }
    }
}

//ip Iterator for BezierLineIter
impl<B, F> std::iter::Iterator for BezierSplitTIter<B, F>
where
    B: Clone,
    F: Num,
{
    /// Item is a pair of points that make a straight line
    type Item = (F, F, B);

    /// next - return the next Bezier
    ///
    /// It pops the first Bezier from the stack; a client may then
    /// consume this Bezier, or it might choose to add_split the Bezier, or whatever.
    fn next(&mut self) -> Option<Self::Item> {
        match self {
            Self::Split { stack } => stack.pop(),
            Self::TRange {
                bezier,
                t0,
                t_scale,
                ref mut step,
                num_steps,
            } => {
                if *step >= *num_steps {
                    None
                } else {
                    let t = *t0 + F::from_usize(*step).unwrap() * *t_scale;
                    *step += 1;
                    Some((t, t + *t_scale, bezier.clone()))
                }
            }
        }
    }
}
