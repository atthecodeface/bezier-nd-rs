//!  A [Bezier] is an implementation of an arbitrary bezier up to a maximum degree
//! given in the type.
//!
//! The control points are stored to be used with Bernstein basis polynomials
//!
//! ``` ignore
//!    b[i,n](t) = n!/(i!(n-i)!) . u^i . t^(n-i) where u=1-t, n=Bezier degree, 0<=i <=n
//! ```
//! A point at parameter t (0 <= t <= 1) is then Sum(b[i,n](t).pts[i])
//!
//! If n were two then this is u^2.pts[0] + 2.u.t.pts[1] + t^2.pts[2], which can be rewritten as:
//! ``` ignore
//!   (1-2t+t^2).pts[0] + 2.(1-t).t.pts[1] + t^2.pts[2]
//!   = t^2.(pts[0] - 2.pts[1] + pts[2]) + 2t.(pts[1] - pts[0]) + pts[0]
//! ```
//!
//! i.e. a set of different point values where a monomial basis (1, t, t^2, ...) is used
//!
//! To convert from Bernstein basis to monomial basis for degree 4 (five points) the pts would map
//! using a 5x5 matrix:
//!
//! ``` ignore
//!     |  1   0   0   0   0 |
//!     | -4   4   0   0   0 |
//!     |  6 -12   6   0   0 |
//!     | -4  12 -12   4   0 |
//!     |  1  -4   6  -4   1 |
//!```
//!
//! Actual value Mij for degree N is (-1)^(j-i) . (n C j) . (j c i) for i<=j, 0 otherwise (by inspection...)
//!
//! Note that this matrix is guaranteed to be invertible, hence all Bernstein representations are equivalent
//! to monomial representations, and vice versa
//!
//! Note that P[N] is only used in conjunction with monomial t^N, and hence if this were to be used
//! to derive points etc then it would be unstable near t=1. However, it has the potentially useful
//! property (for t<1) that 1 > t > t^2 > t^3 > 0.
//!
//! Another possible basis would be for the basis 1, (t-T), (t-T)^2, etc for a fixed T (effectively
//! monomial is this with T=0). This is termed the Taylor basis vectors in https://arxiv.org/pdf/2201.07834
//! This has the property that it is locally (at t=0) approximatable (as one can ignore high powers of (t-T)).
//!
//! Metrics ('distance' between two Beziers)
//!
//! Various metrics (real number) can be used for 'measuring' the difference between two Beziers. Note that a metric
//! is a special kind of measurement (such as the distance between two points), where:
//!
//! * The metric is always >= 0
//!
//! * If the metric is zero then the two items are the same
//!
//! * The metric between A and B is the same as that between B and A
//!
//! * A metric between A and B, plus the metric between B and C, must be less than or equal to that between A and C
//!
//! One metric (maximum parameterwise distance, call it dM) would be the maximum distance between two Beziers
//! for the same parameter t (i.e. max[t] |B(t)-C(t)|).
//!
//! One metric is the square-root of the sum of the distance-squred between the control points of two Beziers
//! (this is termed here dF)
//!
//! Another is the maximum of the difference betwen the control points of two Beziers (termed dC). This
//! has a neat property that both Beziers with a metric will lie within (1) the first Bezier plus a
//! ball (of D dimensions) of radius dC, and (2) within the second Bezier plus a similar ball. Another way to put
//! this is that the maximum distance between any two points on the two Beziers for the same parameter.
//!
//! Note that dF is guaranteed to be >= dC (since it is only *one* control point), but df <= sqrt(N)*dC as there are
//! N points.
//!
//! As more points are added to a Bezier, N increases, and dF can increase; dM, howwver, is independent of N.
//!
//! # Least-squares reduction of degree N+1 to degree N
//!
//! The reduction matrix for degree n+1 to n is given by
//!
//! RL2 = E.transpose() * (E*E.transpose()).inverse()
//!
//! Note that E * RL2 = E * E.transpose * (E*E.transpose()).inverse() = X * X.inverse() = identity
//!
//! # Basis reduction using points on the Bezier
//!
//! The polynomial stored in the structure used the Bernstein polynomials as a basis.
//!
//! That is, point(t) = BernsteinMatrix(t) * Polynomials
//!
//! If we select 'N' different values of t, ti, with (0<=ti<=1), then we find 'N' points pi
//! where pi = Polynomials * B(ti)
//! i.e. [pi] = Polynomials * [B(ti)]
//! or Polynomials = [pi] * [B(ti]].inverse()
//!
//! So if we want to determine Polynomials we can do so from N points pi at N values ti. These are points on the
//! Bezier.
//!
//! To reduce from degree N+1 to N we will require N+1 points (as a degree N curve requires N+1 points). We can
//! use *any* N+1 ti in the range 0<=ti<=1, as long as they are all different values.
//!
//! We must in some sense find pi for these N+1 points given these ti for order N+1, and then map them back to
//! the Bernstein basis for order N.
//!
//! We will get
//!
//!   Polynomials[N] = [pi] * [B[N](ti)].inverse()
//!   Polynomials[N] = Polynomials[N+1] * [B[N+1](ti)] * B[N](ti)].inverse()
//!
//! i.e. we can use a reduction matrix
//!
//!   R[t0..tN] = [B[N+1](ti)] * [B[N](ti)].inverse()
//!
//! Note that elevating-by-one of B[N](t) produces B[N+1](t), so
//!   E[N] * R[t0..tN] = E[N] * [B[N+1](ti)] * [B[N](ti)].inverse()
//!                    = [E[N] * B[N+1](ti)] * [B[N](ti)].inverse()
//!                    = [B[N](ti)] * [B[N](ti)].inverse()
//!                    = Identity
//!
//! i.e. the reduction chosen has the property that elevate(reduce) is the identity, which we need for a reduction

mod accessors;
mod bezier;
mod constructors;
mod evaluation;
mod manipulation;
mod metrics;

pub mod bezier_fns;

pub use bezier::{Bezier, BezierReduceIter};
