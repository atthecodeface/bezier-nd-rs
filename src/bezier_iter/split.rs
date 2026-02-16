use crate::BezierSplit;
use crate::Num;

/// An iterator with Item = (F, F, Bezier), allowing for splitting dynamically during
/// iteration
#[derive(Clone)]
pub enum BezierSplitTIter<B: Clone, F: Num> {
    Split {
        /// A stack of future beziers to examine
        /// The top of the stack is p0->p1; below that is p1->p2, etc
        /// These beziers may need to be split to achieve some critertion
        stack: Vec<(F, F, B)>,
    },
    TRange {
        bezier: B,
        /// Range of t
        t0: F,
        t_scale: F,
        step: F,
        num_steps: F,
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
        let num_steps = F::from_usize(num_steps).unwrap();
        let t_scale = {
            if num_steps <= F::ONE {
                F::ONE
            } else {
                (t1 - t0) / (num_steps - F::ONE)
            }
        };
        let step = F::ZERO;
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
                let tm = (*t0 + *t1) / (2.0_f32).into();
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
                    let t = *t0 + *step * *t_scale;
                    *step += F::ONE;
                    Some((t, t + *t_scale, bezier.clone()))
                }
            }
        }
    }
}
