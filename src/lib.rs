//a Documentation
#![warn(missing_docs)]
/*!

TO DO:

* New building of Bezier
* Reduce to quadratics, reduce to straight lines, ele_red matrix minus identity creation and use for dc metrics

# Bezier curve library

This library provides traits for implementing Bezier curves, such as
evaluation at parameter values, derivatives, splitting into lines, and
reduction/elevation of degree.

Implementations are provided for arrays of vectors, where a vector is an array of D dimensions of floats (such as f32 or f64).

The Beziers are stored using standard control points, and hence use Bernstein polynomials of a parameter 't' to find trace the locus of the Bezier.

# Fundamental numeric types

Traits are used throughout that define the numeric types that Bezier curves use as parameters or for distances etc.
This is fundamentally to provide support for 'f32' and 'f64' as generic floats, but it (in general) also allows for
rational number types - provided they are Copy.

Two traits are used. Firstly there is Num, which is effectively a collection of num_traits, with From<32>, PartialOrd, Copy and Display.
This trait has a blanket implementation for all types that support the required traits. This 'Num' trait is all that is required for the Bezier traits.

Secondly there is Float, which is Num plus the
addition of num_traits::Float - which provides sqrt() and cbrt(). These latter are required for some root finding, or minimum distance
determination. Hence some Bezier *implementations*, such as for [[F;D];4] for cubic Bezier curves, require F to support not just Num but Float, as the
mathematics for the analytical functions require sqrt() or cbrt().

The fundamental types 'f32' and 'f64' implement Num and Float.

# Bezier curve types

The simplest Bezier types are just arrays of points, where points are D-dimensional arrays of a type F.
An array of two points can be used as a linear Bezier; an array of three points can be used as a quadratic Bezier; and array of four points can be used
as a cubic Bezier.

A legacy generic type Bezier is provided which can represent any one of linear, quadratic or cubic Bezier.

A BezierN type is provided that has a generic usize 'N' indicating the *maximum* degree the individual curve can have, but the
type also includes a degree that the curve actually has. The type stores its control points in an array [[F;D];N], and so
every Bezier curve the type can describes requires the same memory size.

# Bezier traits

Bezier curves have many uses, and there are many ways in which they may be manipulated and used. A variety of traits
are provided by this library to encompass those use cases. Some example applications are:

* Font outlines (quadratic or potentially cubic Beziers)
* Drawing programs (generally cubic or lower, but potentially arbitrary degree)
* Path description for postiion of objects (such as robotic joints), with arbitrary degreee Beziers

The broad category of functionality required is:

* Evaluation of points and derivatives
* Arc length calculation
* Curve simplification (reduction to lines, lower Beziers) with some tolerance
* Curve splitting
* Minimum distance from a curve, and closest point (point projection)
* Reparametrization of a curve
* Curvature
* Intersection with a line
* Self-intersectoin
* Curve building (from points and derivatives)
* Bounding box

The library's traits expect Bezier's to use a parameter 't' of some type 'F' which must be similar to a float (e.g. f32, f64) and points that are clonable.
Normally the trait implementation will be for D-dimensional points with coordinates of type 'F'.

## BezierEval

The simplest trait provided is [BezierEval], which is a dyn-compatible trait that allows a type to provide:

* Access to the degree and control points of a Bezier
* Evaluation of points and first derivatives at a parameter value 't'
* Determination of how close to a straight line a Bezier is
* Determination of the worst case distance between points on the Bezier and a Bezier reduced simply to a quadratic or cubic Bezier with the same endpoints

This is implemented for:

* [[F; D]; 2] as a linear Bezier, with parameter of type F and points of type [F;D]
* [[F; D]; 3] as a quadratic Bezier, with parameter of type F and points of type [F;D]
* [[F; D]; 4] as a cubic Bezier, with parameter of type F and points of type [F;D]

## BezierSplit

This trait provides a method to split a Bezier curve into two of the same degree; it requires the type
to be Sized, and is not dyn compatible as the split method must return two Beziers.

Splitting a Bezier curve into two (curves A and B, at the point with t=0.5) yields two Beziers that provide precisely the
same points; the first for P(t)=Pa(2*t) 0<=t<=0.5, and the second for P(t)=Pb(2*t-0.5) 0.5<=t<=1.

If a type provides this trait then a BezierSplitIter can be created, which is an iterator that can provide
for recursive splitting of the Bezier; in a loop such as

```ignore
let mut split = bezier.split_iter();
while let Some(b) = split.next() {
  if some_criterion {
    // Argh! Bezier 'b' is not refined enough
    //
    // Recurse using b split into two Bezier curves at t=0.5!
    split.add_split(b);
  } else {
    ...
  }
}
```

An example of this provides for splitting a Bezier curve into an iterator of connected line segments, if the
Bezier supports Clone, BezierSplit and BezierEval, such as:

```ignore
fn draw_bezier(bezier:&Bez, tolerance_sq:f32) {
  for (p0, p1) in BezierLineIter::new(&bezier, tolerance_sq) {
    // Next line in approximation of bezier, connected to the previous line (if any)
    // is p0 -> p1
    draw_line(p0, p1);
  }
}
```

## BezierReduce

This trait is again not dyn-compatible; it provides methods for reducing a Bezier from one degree to
a lower degree, and potentially down to a cubic or quadratic. The latter methods expect the reduced Bezier
to have the same endpoints (although this is not an absolute guarantee, if supported it enables the
use of recursive splitting of a high degree Bezier into connected cubic or quadratic Bezier curves).

The trait also has methods for determining how close the Bezier is to its reductions (be they by one degree,
or to quadratic or cubic). With these methods, an algorithm can recursively split a Bezier into smaller Beziers (if BezierSplit
is provided for the type), until a smaller Bezier is close enough to a cubic (or quadratic, etc) so that
the initial Bezier ends up being split into a list of connected cubic (or quadratic etc) Bezier curves.

A Bezier that supports this trait is normally expected to support reduction; however, it can be useful to implement this for
linear Beziers, where a reduction makes no sense, and so there is a 'can_reduce' method here that should be used to gate
the reduction when required.

## BezierDistance

This trait provides methods that can estimate and potentially can accurately provide the (minimum) distance (squared) between
a point P (of the type the Bezier supports) and the Bezier.

A quadratic Bezier is expected to provide an accurate distance to a point; a linear Bezier (a line segment) *must* provide
an accurate distance. Higher degree Bezier curves need not (i.e. can return None from such methods).

If the Bezier supports BezierSplit and BezierReduce, then the distance of a point from the Bezier can be determined
to a tolerance - by splitting/reducing the curve to the tolerance into quadratics or lines, and then determining the
minimum distance from the point to the quadratics or lines.

The estimation method in this trait should be relatively light-weight; the aim is that a set of Bezier curves
(even of high degree) that support BezierSplit can be iterated over to find the distance between a point P
and the curves (and which curve, and where on the curve). This allows extending the distance-to-Bezier algorithm, given
a current distance to a closest Bezier, to test if *another* Bezier is potentially closer (i.e. its estimation
method returns a value less than the closest current distance) to see if it should be split (or if it is quadratic or linear
then a precise determination used).

The estimation method must therefore always provide a result that is less than or equal to the *actual* distance from
the point to the Bezier.

## BoxedBezier (requires BezierEval)

This trait is explicitly dyn-compatible, and it provides for a mechanism for splitting and reducing Beziers
(akin to BezierSplit and/or BezierReduce) - into Boxed dyn BoxedBezier. If a type does not support reduce
but does support split then it can implement this trait; the methods all return Option of a Boxed object.

## Legacy

This library provides linear, quadratic and cubic Bezier curves, using
generic types for the scalar that they are made up of; this can be
either f32 or f64 (anything that supports the geo_nd::Float trait).

Points in the Bezier are [Float; D] for a D-dimension Bezier curve.

## Bezier types

A linear Bezier has two points, p0 and p1, and provides points
along the line as:

```text
  p(t,u=1-t) = u*p0 + t*p1
```

A quadratic Bezier has three points, p0, c and p1, and provides
points along the curve as:

```text
p(t,u=1-t) = u^2.p0 + 2.t.u.c + t^2.p1
```

or, viewing it is a linear Bezier between two linear Beziers:

```text
p(t) = u(u.p0 + t.c) + t(u.c + t.p1)
```

Cubic Beziers go one step further, being effectively linear Bezier
curves of two quadratic Beziers, but explicitly a point at parameter t
is:


```text
p(t,u=1-t) = u^3.p0 + 3.t.u^2.c0 + 3.u.t^2.c1 + t^3.p1
```

## Overview of Bezier curve type

The Bezier type supports construction and interrogation of the type
instance.  The type utilizes a `num_pts` field, which is 2 for linear
Beziers, 3 for quadratic, and 4 for cubic.  The type includes an array
of *four* control points for *every* Bezier - the first two points are
the endpoints of the Bezier (p0 and p1, in notation). The third
control point is not used for linear Beziers, and the fourth control
point is only used for cubic Beziers.

The Bezier can be interrogated, to determine any point along the curve
given the parameter `t` (from 0 to 1), as can its tangent (which is a
simple linear combination of the control points - the derivative of
the Bezier function with respect to `t`).

Any portion of a Bezier curve between two values of `t` is also a
Bezier curve of the same type; hence a method is supplied to permit
sectioning of a Bezier into another Bezier.

## Straightness and Bisection

The Bezier curve can be bisected in to two Bezier curves of the same
type (such as cubic); the two halves of a Bezier will be closer
approximations to a straight line, and hence to render a Bezier it is
common to bisect a Bezier recursively until the elements are 'straight
enough'. One can consider a cubic Bezier curve, for example being
mapped to a list of linear Beziers which join end-to-end, which
approximate the original Bezier to within a 'straightness' measure.

This concept can be used not only to render a Bezier as straight line
segments, but also to approximate the length of a Bezier, or to find a
value for the parameter `t` for a particular distance along the
Bezier.

Methods are provided to generate iterators that will yield straight
line segments or points, which approximate a Bezier curve to a
particular straightness; further methods are supplied to calculate the
length of a Bezier, or to find the value of `t` for s distance along
the Bezier.

The actual straightness measure used for a Bezier curve in this
library is the sum of the square of the distances of the intermediate
control points from the straight line between the endpoints.

For example, if the units of the Bezier coordinates are pixels then a
straightness of 1.0 will produce straight line segments that never
exceed 1 pixel away from the original Bezier.

## Circular arcs and rounding

Bezier curves cannot express circular arcs precisely; they are
different mathematical beasts. However, a circular arc may be
approximated quite closely with a cubic Bezier curve.

It can be useful, when using two-dimensional graphics, to generate
rounded rectangles (or other polygons). For these case the corner
points of the rectangles are known, and the directions of the lines
into the corner are known, and it is useful to generate a Bezier that
represents the arc which can replace this corner. This is another
circular arc, but described by the corner, approach vectors, and
rounding radius.

# Examples

An example using 2-dimensional vectors of f32

```
type Bezier = bezier_nd::Bezier<f32, 2>;
use bezier_nd::{BezierEval, BezierSplit};

let dx = [1.,0.];
let dy = [0.,1.];
let line = Bezier::line( &dx, &dy );
assert_eq!( line.degree(), 1);
assert!( line.is_straight(0.));

// A 30-degree arc of radius 3; has length of PI/2,
// it is fairly close to being a straight line, subjectively...
let arc = Bezier::arc( 30.0f32.to_radians(), 3.0, &[3.,4.].into(), &dx, &dy, 0.);
assert_eq!( arc.degree(), 3);
assert!(!arc.is_straight(0.));

// Breaking the arc with a large value of straightness yields only 3 points:
assert_eq!( arc.as_points(0.01).count(), 3);
// This is 1.5662874
println!( "Arc length when straightened to '0.1' is {}", arc.length(0.1) );

// Breaking the arc with a small value of straightness yields 17 points!
assert_eq!( arc.as_points(0.00001).count(), 9);
// This is 1.5707505 - a lot closer to PI/2 = 1.57079632679
println!( "Arc length when straightened to '0.1' is {}", arc.length(0.1) );

for (a,b) in arc.as_lines(0.05) {
    println!("Line from {a:?} to {b:?}\n");
}

let q = Bezier::quadratic( &[2.,6.].into(),
                           &[3.5,8.].into(),
                           &[4.,7.].into());
assert_eq!( q.control_point(0)[0], 2., "Start point X of Bezier" );
assert_eq!( q.control_point(0)[1], 6., "Start point Y of Bezier" );
assert_eq!( q.control_point(1)[0], 4., "End point X of Bezier" );
assert_eq!( q.control_point(1)[1], 7., "End point Y of Bezier" );
```

An example using 3D vectors of f64

```
type Bezier = bezier_nd::Bezier<f64, 3>;

let c = Bezier::cubic( &[1.,0.,0.],
                           &[2.5,0.,-1.],
                           &[0.,2.5,1.],
                           &[0.,1.,0.]);
// This is just over 3.283
println!( "3D cubic length when straightened to '0.1' is {}", c.length(0.1) );
// But this is just over 3.29
println!( "3D cubic length when straightened to '0.01' is {}", c.length(0.01) );

```

# Circular arc algorithm

There are many papers and articles on constructing Bezier curves that
approximate circular arcs; from an academic perspective these provide
analytical approximation generation. This library is more practical:
it provides an implementation that is simple, fast, and works for all
angles. Effectively it is table-based; in fact, it uses a quartic
polynomial of the square of the sine of the angle of the arc, with the
polynomial matching a table generated from experimental data.

The approach taken was to start with an analytical value for the
coefficient (`lambda`) applied to the corner direction vectors to
generate the two control points, and to determine the
mean-square-error this generated, for a number of points along the
Bezier curve, in the distance from the point to the centre of the arc
compared to the radius of the arc. The value of 'lambda' was adjusted
in small steps, to reach a minimum mean-square-error for the angle of
the arc.

Numerous tables of values of best-lambda compared to different
functions of the angle (such as sin(angle), sin(angle)^2, etc) were
produced, along with reciprocal tables, and a best low-order
polynomial match was sought which would yield a very high correlation
to the table.

This polynomial is encoded in the Bezier function - currently not
particularly well optimized, though. The accuracy of the Bezier curves
produced is compared in regression testing to be within 0.1% (using
the mean-square-error of a number of points along the Bezier).

This is not looking for the worst excursion; no visual comparisons
have been done; the purpose is to provide something that is clean,
works across the range, and is efficient. There is no perfect solution
to this problem!

# References

A Primer on Bézier Curves ("Pomax")

https://pomax.github.io/bezierinfo/

Adaptive Bezier Degree Reduction and Splitting
for Computationally Efficient Motion Planning (Omur Arslan and Aron Tiemessen):

https://arxiv.org/pdf/2201.07834

!*/

/*a Imports
*/
pub mod bernstein;
mod builder;
pub(crate) mod constants;
pub(crate) mod utils;

mod traits;
pub use traits::{Float, Num};

mod curve;
mod distance;
pub(crate) mod polynomial;

mod bezier_iter;
mod farray_cubic;
mod farray_line;
mod farray_quadratic;

pub mod bignum;

/*a Exports
*/
pub use bezier_iter::{BezierLineIter, BezierPointIter, BezierQuadIter, BezierSplitIter};
pub use traits::{BezierDistance, BezierEval, BezierReduce, BezierSplit, BoxedBezier};

pub use builder::BezierBuilder;
pub use constants::*;
pub use curve::Bezier;
pub use distance::BezierDistance2D;
pub use polynomial::{PolyFindRoots, PolyNewtonRaphson, Polynomial};

pub use bernstein::Bezier as BezierND;
