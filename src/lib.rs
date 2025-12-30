//a Documentation
#![warn(missing_docs)]
/*!

# Bezier curve library

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
assert_eq!( arc.as_points(0.1).count(), 3);
// This is 1.5662874
println!( "Arc length when straightened to '0.1' is {}", arc.length(0.1) );

// Breaking the arc with a small value of straightness yields 17 points!
assert_eq!( arc.as_points(0.001).count(), 17);
// This is 1.5707505 - a lot closer to PI/2 = 1.57079632679
println!( "Arc length when straightened to '0.1' is {}", arc.length(0.1) );

for (a,b) in arc.as_lines(0.05) {
    println!("Line from {a:?} to {b:?}\n");
}

let q = Bezier::quadratic( &[2.,6.].into(),
                           &[3.5,8.].into(),
                           &[4.,7.].into());
assert_eq!( q.borrow_pt(0)[0], 2., "Start point X of Bezier" );
assert_eq!( q.borrow_pt(0)[1], 6., "Start point Y of Bezier" );
assert_eq!( q.borrow_pt(1)[0], 4., "End point X of Bezier" );
assert_eq!( q.borrow_pt(1)[1], 7., "End point Y of Bezier" );
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

!*/

/*a Imports
*/
mod curve;
mod line;
mod point;

/*a Exports
*/
pub use self::line::BezierLineIter;
pub use self::point::BezierPointIter;
pub use curve::Bezier;
