//a Documentation
#![warn(missing_docs)]
/*!

# TODO

Adder BuilderError to BezierBuilder, to provide feedback on why a build fails

Add tests for:
* t_closest
* reduce
* elevate_by(>1)
* BezierOps
* Polynomial (set_multiply, set_divide, find cubic real root, tests for [F])
* section_length with t0<0, t1>1, t1<t0
* Approximation t_of_distance beyond length
* Approximation derivative_at of t<0 or t>1
* Approximation into iterator
* Appromixation for_each_control_point

* add 'BezierCubic' trait? - support arc
* add 'BezierQuadratic' trait?
* Amalgamate BezierSplit and BezierSection
* Move BezierOps, Distance, MinMax, etc (all dyn compatible) to BezierEval If Approxmiation can do them
* Add c_sq, f_sq, df_sq_from_line to BezierEval

# Bezier curve library

This library provides traits for implementing Bezier curves, such as
evaluation at parameter values, derivatives, splitting into lines, and
reduction/elevation of degree.

Implementations are provided for arrays of points], where a point is an array of D dimensions of floats (such as f32 or f64), and `Vec` of points.

The Beziers are stored using standard control points, and hence use Bernstein polynomials of a parameter 't' to find trace the locus of the Bezier.

# Fundamental numeric types

Traits are used throughout that define the numeric types that Bezier curves use as parameters or for distances etc.
This is fundamentally to provide support for 'f32' and 'f64' as generic floats, but it (in general) also allows for
rational number types - provided they are Copy.

Two traits are used. Firstly there is [Num], which is effectively a collection of num_traits, with From<32>, PartialOrd, Copy and Display.
This trait has a blanket implementation for all types that support the required traits. This [Num] trait is all that is required for the Bezier traits.

Secondly there is [Float], which is Num plus the
addition of [num_traits::Float] - which provides sqrt() and cbrt(). These latter are required for some root finding, or minimum distance
determination. Hence some Bezier *implementations*, such as for `[[F;D];4]` for cubic Bezier curves, require F to support not just Num but Float, as the
mathematics for the analytical functions require sqrt() or cbrt().

The fundamental types 'f32' and 'f64' implement Num and Float.

# Bezier curve types

The simplest Bezier types are just arrays of points, where points are D-dimensional arrays of a type F.
An array of two points can be used as a linear Bezier; an array of three points can be used as a quadratic Bezier; and array of four points can be used
as a cubic Bezier.

A legacy generic type Bezier is provided which can represent any one of linear, quadratic or cubic Bezier.

A BezierN type is provided that has a generic usize 'N' indicating the *maximum* degree the individual curve can have, but the
type also includes a degree that the curve actually has. The type stores its control points in an array `[[F;D];N]`, and so
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
* Curve simplification (reduction to lines, lower Beziers) within some tolerance
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

## [BezierEval]

The simplest trait provided is [BezierEval], which is a dyn-compatible trait that allows a type to provide:

* Access to the degree and control points of a Bezier
* Evaluation of points and first derivatives at a parameter value 't'
* Determination of how close to a straight line a Bezier is
* Determination of the worst case distance between points on the Bezier and a Bezier reduced simply to a quadratic or cubic Bezier with the same endpoints

This is implemented for:

* `[[F; D]; 2]` as a linear Bezier, with parameter of type F and points of type `[F;D]`
* `[[F; D]; 3]` as a quadratic Bezier, with parameter of type F and points of type `[F;D]`
* `[[F; D]; 4]` as a cubic Bezier, with parameter of type F and points of type `[F;D]`
* `BezierND<F, N, D>` for a Bezier up to degree N with parameter of type F and points of type `[F;D]`
* `Vec<[F;D]>` as an aribtrary degree Bezier with parameter of type F and points of type `[F;D]`

## [BezierSplit]

This trait provides a method to split a Bezier curve into two of the same degree; it requires the type
to be [Sized], and is not dyn compatible as the split method must return two Beziers.

Splitting a Bezier curve into two (curves A and B, at the point with t=0.5) yields two Beziers that provide precisely the
same points; the first for `P(t)=Pa(2*t)` (0<=t<=0.5), and the second for `P(t)=Pb(2*t-0.5)` (0.5<=t<=1).

If a type provides this trait (and Clone) then a [BezierSplitIter] can be created, which is an iterator that can provide
for recursive splitting of the Bezier; in a loop such as

```ignore
let mut split = BezierSplitIter::new(&bezier);
while let Some(b) = split.next() {
  if some_criterion {
    // Argh! Bezier 'b' is not refined enough
    //
    // Recurse using b split into two Bezier curves at t=0.5!
    split.add_split(b);
  } else {
    // Yay! Use b for the operation - it is refined enough
    // for the purpose
    ...
  }
}
```

The [BezierIntoIterator] trait enhances types that implement [BezierSplit], to provide methods for iterating over points
or lines to within a degree of accuracy. For example, a Bezier curve can be drawn with a desired accuracy with:

```ignore
fn draw_bezier<B:BezierEval<f32, [f32; 2]> + BezierSplit + Clone>(bezier:&B, tolerance_sq:f32) {
  for (p0, p1) in bezier.as_points(tolerance_sq) {
    // Next line in approximation of bezier, connected to the previous line (if any)
    // is p0 -> p1
    draw_line(p0, p1);
  }
}
```

## [BezierElevate]

This trait allows a Bezier curve to be elevated by either one degree, or many.

An elevated Bezier curve is a curve of a higher degree that provides precisely the same points.

## [BezierReduce] (requres [BezierEval])

This trait provides methods for reducing a Bezier from degree 'N' to
a degree 'N-1', and potentially down to a cubic or quadratic.

Bezier curves do not have a standard mechanism for reduction of degree, as the process is lossy.
These reduction mehods *expect* the reduced Bezier
to have the same endpoints, as this enables the
use of recursive splitting of a high degree Bezier into connected Beziers of the reduced degree.

The trait also has methods for determining how close the Bezier is to its one-degree reduction; the [BezierEval]
trait includes methods for how close the Bezier is to its reduction to lines, quadratic Bezier, and cubic Bezier.

With these methods, an algorithm can recursively split a Bezier that also supports [BezierSplit] into shorter Beziers,
until the shorter Bezier is close enough to a cubic (or quadratic, etc) so that it can be replaced with the reduction.
In this way the initial Bezier ends up being split into a list of connected cubic (or quadratic etc) Bezier curves.

```ignore
fn split_into_cubics<B>(bezier: &B, tolerance_sq: f32)
where
    B: BezierSplit + BezierReduce<F, P> + Clone,
{
    let mut result: Vec<_> = vec![];
    let mut split = BezierSplitIter::new(&bezier);
    while let Some(b) = split.next() {
        if b.closeness_sq_to_cubic() > tolerance_sq {
            // Bezier 'b' is not refined enough
            //
            // Recurse using b split into two Bezier curves at t=0.5!
            split.add_split(b);
        } else {
            // Assuming we *know* that the reduction is possible (i.e. does not return None)
            result.push(b.reduced_to_cubic().unwrap());
        }
    }
}
```

A Bezier that supports this trait is normally expected to support reduction; however, it can be useful to implement this for
linear Beziers, where a reduction makes no sense, and so there is a 'can_reduce' method here that should be used to gate
the reduction when required.

## [BezierDistance]

This trait provides two methods: one to estimate the distance squred of a point from the Bezier, erring on the side of underestimation;
the other to potentially accurately provide the minimum distance squared between a point and the Bezier and the parameter 't' of the point
that is this distance.

A quadratic Bezier is expected to provide an accurate distance to a point; a linear Bezier (a line segment) *must* provide
an accurate distance. Higher degree Bezier curves need not (i.e. can return None from such methods).

If a high degree Bezier supports [BezierSplit] and [BezierReduce], then the distance of a point from the Bezier can be determined
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

## [BezierConstruct] (requires Sized)

The [BezierBuilder] provides a means to describe key points and derivatives of Bezier curves, from which
an actual curve can be derived. Implementing this trait permits the construction of a Bezier from a builder.

## [BasicBezier] (requires BezierEval, BezierSplit, BezierSection, BezierOps)

## [BezierOps]

Affine operations on Bezier curves - addition and subtraction of one curve to/from another, scaling of a curve, and
mapping of the control points of the curve through a function

## [BoxedBezier] (requires BezierEval and alloc)

This trait is explicitly dyn-compatible, and it provides for a mechanism for splitting and reducing Beziers
(akin to BezierSplit and/or BezierReduce) - into Boxed dyn BoxedBezier. If a type does not support reduce
but does support split then it can implement this trait; the methods all return Option of a Boxed object.

# Types that support Bezier traits

The simplest Bezier curve types supported by the library are arrays of points consisting of 'D' numbers (types that support Num, such as f32, f64,
and even isize) of length 2 (linear Bezier) to 4 (cubic Bezier).

## `[[F;D]; 2]` - linear Bezier using float arrays, with point `[F;D]`

Simple arrays of twp points, using the two points as the endpoints of the Bezier.

This supports [BezierEval], [BezierSplit], [BezierElevate], [BezierReduce], [BezierDistance], [BezierConstruct], [BoxedBezier].

## `[[F;D]; 3]` - quadratic Bezier using float arrays, with point `[F;D]`

Simple arrays of three points for a quadratic Bezier, with the end points being `[0]` and `[2]`.

This supports [BezierEval], [BezierSplit], [BezierElevate], [BezierReduce], [BezierDistance], [BezierConstruct], [BoxedBezier].

## `[[F;D]; 4]` - cubic Bezier using float arrays, with point `[F;D]`

Simple arrays of four points for a cubic Bezier, with the end points being `[0]` and `[3]`.

This supports [BezierEval], [BezierSplit], [BezierReduce], [BezierDistance], [BezierConstruct], [BoxedBezier].

For [BezierDistance] this does not provide an accurate 'closest point' estimation - i.e. t_dsq_closest_to_pt returns None.

## [BezierND<F, N, D>] - 'Copy' Bezier type of maximum degree N with point `[F;D]`

A structure that contains an array of N points (of `[F;D]`), and a degree

This provides a single type that can support arbitrary (up to a maximum N) Bezier degree.

This supports [BezierEval], [BezierSplit]

It SHOULD also support [BezierReduce], [BezierDistance], [BoxedBezier].

## `Vec<[F;D]>` - arbitrary Bezier with point `[F;D]`

The Bezier traits are provided for `Vec<[F;D]>` to provide for arbitrary degree Bezier curves,
when memory allocation is supported and when the application only requires `Clone` not copy.

This supports [BezierEval], [BezierSplit]

It should support [BezierReduce], [BezierDistance], [BoxedBezier].

For [BezierDistance] this does not provide an accurate 'closest point' estimation - i.e. t_dsq_closest_to_pt returns None.

## Note no implementation for `F` as a point

An implementation of (e.g.) [BezierEval] for `F:Num` as the point is not possible, as we require
this for `[F; D]` and yet in the future one might implement `Num` for arrays, and then there would be
two implementations and Rust would error

# Bezier building

It can be useful to construct Bezier curves, but not simply from control points. For example,
if some points are known to be on the Bezier, then a Bezier curve that passes through
those points may be required. Of course, there are infinitely many Bezier curves that do this,
but provided with `N` points *and* associated values for parameter `t`, there is a single Bezier
curve of degree `N-1` that fits.

Hence there is a [BezierBuilder] constructor, which can be provided with points and associated
values of `t`; from this a Bezier can be built.

Sometimes, though further contraints are required - such as the derivative of the Bezier at a parameter is a certain value.
The builder provides for specifying the `nth` derviative at a parameter `t` is a certain value. This forms an additional constratint.

A builder that has `N` constraints defines a Bezier of degree `N-1`; note that a Bezier of degree `N-1` has at `N-1`
non-zero derivatives, and so if an 'nth derivative at t' constraint is provided, there must be at least 'n' other constraints.

Bezier types may provide construction from a [BezierBuilder].

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
assert_eq!( q.endpoints().0[0], 2., "Start point X of Bezier" );
assert_eq!( q.endpoints().0[1], 6., "Start point Y of Bezier" );
assert_eq!( q.endpoints().1[0], 4., "End point X of Bezier" );
assert_eq!( q.endpoints().1[1], 7., "End point Y of Bezier" );
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

# Bezier curve background

A linear Bezier has two points, p0 and p1, and provides points
along the line as:

```text
  p(t,u=1-t) = u*p0 + t*p1
```

A quadratic Bezier has three control points, p0, c and p1, and provides
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

## Other mathematical representations (including *Monomial*)

Bezier curves with `n` control points as defined in the previous section
are described by:

```text
p(t) = Sum( (1-t)^(n-i)*t^i*n!/(i!*(n-i)!) * p[i] )
```

The coefficients of the control points are known as Bernstein basis polynomials:
from here on they will be written as `B[i,n](t)`:

```text
B[i,n](t) = (1-t)^(n-i) * t^i * n!/(i!*(n-i)!)
```

The equation for p(t) can be expanded out into a polynomial in t; each coefficient of `t^i` would
be some linear combination of the control points. Then the point p(t) would be

```text
p(t) = Sum( t^i * m[i] )
```

where the `m[i]` are those linear combinations as required. A Bezier can be represented
by these control points (described as *monomial* curves in one reference) instead of the
standard Bernstein control points.

## Bernstein and Monomial equivalence

Clearly one can convert (as described above) from a set of `N` Bernstein polynomial control points
to a set of monomial control points; it is also possible to convert in the other direction.

One way to help visualize why this is possible is to look at the coefficient of `t^0` in the (expanded)
Bernstein polynomial basis - this comes only from the first control point coefficient `(1-t)^0`, as the
other coefficients all contain `t^i` where i is at least one. Hence the *monomial* control point for `i=0`
can be used directly to find the first Bernstein polynomial control point.

The second control point coefficient (for `i=1`) is similarly only dependent on the first two Bernstein
polynomial control points, and hence (given the first one has been calculated) the second can now be calculated.

Hence Bernstein and monomial representations are equivalent - they can both express the same curves.

## Advantages of Bernstein representation

One of the interesting things about Bernstein polynomials `B[i,n](t)`, for the `i`th polynomial
in a degree `n` curve, at parameter 't', is that the all of them for `0<=i<n` and 0<=t<=1 are between 0 and 1,
and the sum of them for all i will be 1. What this means is that, for 0<=t<=1, every point
on a Bezier (using Bernstein representation, the standard approach) with 0<=t<=1 is a linear combination
of the control points with coefficients between 0 and 1, and coefficients all summing to 1: hence all the points
on the curve must lie within the bounding box of the control points (to go outside the bounding box requires negative
coefficients or coefficent sum greater than 1).

## Legacy

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

# Construct a Bezier curve from points and `t` values

A Bezier curve can be constructed purely from a given set of points,
each with an associated `t` value (where all the `t` values are different).

How can this be done?

Consider the set of points `P[i]` at the values `t[i]`; ideally we would have some function `f[i](t)` that would be 1 at `t[i]`, and 0 at every other
`t[j]` (i.e. when i!=j) - then we could just add up `P[i] * f[i](t)` - this is guaranteed to have the value of `P[i]` at each `t[i]`. What are the
(N) functions `f[i]`? Well, a function that is zero at `t[j]` when i!=j (call it `g[i](t)`) is:

`g[i](t) = (t-t[0]) . (t-t[1]) . (t-t[2]) ... (t-t[i-1]) . (t-t[i+1] ... (t-t[N-1])`

If that is multiplied out, it is clearly a polynmial of degree `N-1`. However, we are looking for `f[i](t)` that is 1 at `t[i]` and 0 at the other t;
this function (g) is `g[i](t[i])` at `t[i]` and 0 at the other t. So, it would seem that:

`f[i](t) = g[i](t) / g[i](t[i])`

Note that `g[i](t[i])` is a constant, and so `f[i](t)` is also a polynomial of degree `N-1`. Adding up all the polynomials `P[i]*f[i](t)` yields `P(t)`, which will also
be a polynomial of degree `N-1`; indeed, it is a *monomial* representation of the Bezier curve that goes through the points `P[i]` at associated parameter values `t[i]`, as required.

Since this provides a way to construct a *monomial* representation of a Bezier curve, and they are equivalent to standard (Bernstein) Bezier curves, it must
be possible to construct a standard Bezier curve from N points each for a given `t` value.

It is worth noting that the proof here that a Bezier can be found that goes through the points as required is a construction proof - i.e. it constructs a Bezier wih N control points;
This is the degrees-of-freedom that such a Bezier curve of degree `N-1` has, and hence this constructs the *only* Bezier of degree `N-1` that goes through the points as required.

# Elevation of Bezier curves

An interesting property of a Bezier curve of degree `N` is that it can be represented precisely by a Bezier curve of degree `N+1`.

The proof of this is obvious if the Bezier curve were to be expressed using *monomial* control points (rather than the standard Bernstein
polynomial control points) - because the *monomial* control points for degree `N+1` are the same as for degree `N` with the addition of one more
(for the new degree `N+1` at the origin). As noted above, *monomial* and standard Bernstein control point representations are equivalent -
hence a Bernstein (or standard) Bezier curve *can* be raised by a degree by some means to another Bezier which describes the same points for
the same parameter values of t.

This property can be useful for combining Bezier curves - for example for finding the *difference* between two curves.

# Addition/Subtraction etc of two Bezier curves

Bezier curves are mathematically a polynomial in its parameter 't' with
coefficents that are the control points. As such, they can be scaled, summed,
differenced, etc by applying the operations to the control points. This is perhaps
most meaningful when the polynomials, and hence the Beziers, are of the same
degree (for binary operations such as addition and subtraction).

If two Beziers are of different degree, though, it is always possible to elevate the one of
lower degree to the same as the other, without impacting the Bezier curve - hence binary operations
can always be performed effectively on Bezier curves of *any degree*.

Furthermore, if two Beziers of the same degree have the same endpoints, then the difference between
the two curves will be a Bezier curve of the same degree (since it is a polynomial of the same degree) whose
endoints are both at the origin.

# Maximum distance of a Bezier curve from the origin

Assuming that a Bezier curve is defined using Bernstein polynomials - the standard way that
the control points - and as noted above - every point on the curve with 0<=t<=1 lies within
the bounding box of the control points.

Another way to express this is that no point on the curve can be further away than from the origin than the furthest control point.

If the Bezier curve happens to have both endpoints at the origin (as it might if it were the difference
between two Bezier curves with the *same* endpoints) then its maximum excursion from the origin
is less than or equal to the distance from the origin of the furthest control point.

# Metrics (measures) of Bezier curve differences

A `metric` is a mathematical term that has some rigor; it is a measurement (real number) which has the properties that:

* it is non-negative
* it is zero only for identical Bezier curves (in this case - i.e. the describes the same points for the same 't')
* it is symmetric - if `A->B` has metric M then `B->A` has metric M
* it satisfies the 'triangle inequality' - i.e. `A->C` <= `A->B` + `B->C`

A good example is the Euclidean measure of distance between two points in a Euclidean space.

Many metrics of difference might be a metric of the *difference* between the two curves - which itself is as Bezier.
Thus these metrics are a measure that can be applied to a single Bezier. This is similar to the Euclidean distance between
two points, which is the same as the length of the vector between them (the difference between vectors from the origin to those points).

(To find the difference between two Bezier curves of different
degree, elevate the first Bezier (raising its degreee) to the same as the other first.)

Below the standard representations for a Bernstein Bezier `Sum[i](P[i]*B[i,n](t))` is used, i.e. control points at `P[i]`.

## L2-norm

One obvious metric is the total squared distance between the origin and the Bezier curves for every point `t`, 0<=t<=1.
This is the integral (with respect to `t` for 0<=t<=1) of `|Sum[i](P[i]*B[i,n](t))|^2`.

```text
|Sum[i](P[i]*B[i,n](t))|^2
 = ( Sum[i](P[i]*B[i,n](t)) ) . ( Sum[j](P[j]*B[j,n](t)) )
 = ( Sum[i,j](P[i].P[j} * B[i,n](t) * B[j,n](t) )
```

So the integral would come:

```text
Integral(t,0,1) [ |Sum[i](P[i]*B[i,n](t))|^2 ]
 = Integral(t,0,1) [ Sum[i,j]( P[i].P[j} * B[i,n](t) * B[j,n](t)) ]
 = Sum[i,j]( P[i].P[j} * Integral(t,0,1) [ B[i,n](t) * B[j,n](t) ] )
```

i.e. the L2 norm is a weighted sum of the dot products of the control points.

## Maximum Excursion From the Origin (l_sq)

Another obvious metric is the square of the *maximum* excursion from the origin of a Bezier, i.e.

```text
Max(Sum[i](P[i]*B[i,n](t))|^2)
 = Max[t][ Sum[i,j]( P[i].P[j} * B[i,n](t) * B[j,n](t) ) ]
```

This cannot be determined analytically, but it is guaranteed to be larger than the L2 norm.

In [BezierND] this is referred to as the `m_sq` metric.

Note that m_sq is the maximum squared distance of the Bezier from the origin; for a Bezier that
is the difference between two Bezier curves `A` and `B`, this metric (`dm_sq` for m-of-the-difference) bounds
the maximum difference between the two curves - i.e. every point on Bezier `A` is within
`dm_sq` of the equiavalent (same t) point on Bezier `B`.

## Total Control Point Length Squared

The L2-norm, which is a weighted sum of the dot products of the control point pairs, inspires
a weighted sum with 1 for `i=j`, and 0 otherwise - i.e. the sum of the length squared of all
of the control points for the Bezier - `Sum[i](P[i].P[i])`.

In [BezierND] this is referred to as the `f_sq` metric.

It is worth noting that the `f_sq` metric of an *elevated* Bezier (which describes the same
points) can be greater than that of the original Bezier (i.e. elevation may increase the `df_sq` metric)

## Maximum Control Point Length Squared

It was noted above that the maximum excursion from the origin is guaranteed to be no more than
the distance of the furthest control point from the origin, so the metric `Max(P[i]^2)` may be useful.
In [BezierND] this is referred to as the `dc` metric.

It is worth noting that the 'c_sq' metric is related to the `f_sq` metric somewhat, by `f_sq` being the *sum* of N values, and `c_sq` being the
*maximum* of those N values. Also, the `c_sq` metric of an *elevated* Bezier (which describes the same
points) is guaranteed to be no more than that of the original Bezier (i.e. elevation does not increase the `dc_sq` metric)

## Relationships between the metrics

The metrics above are ordered by magnitude in the following manner:

`l_sq <= m_sq <= c_sq <= f_sq <= N*c_sq`

Note that `c_sq` and `f_sq` can be calculated analytically; `m_sq` is the metric that might
be most useful, as it bounds the point-to-point difference between two Bezier curves, but it has
no reasonable analytical form.

# Reduction of Bezier curves

A question might be posed: 'If a Bezier curve of degree `N` can b elevated to one of degree `N+1` without loss, can it be reduced
without loss too (for any arbitrary curve)?'

Well, clearly the answer is 'no' for *monomial* representations - as to do so would mean throwing away the control point
associated with `t^(N+1)`, which would lead to different points everywhere on the curve at t!=0 unless that control point
itself was the origin.

Since *monomial* representations are equivalent to standard Bezier curves, this
means that standard Bezier curves cannot be reduced without loss in most cases.

So how does one approximate a Bezier of degree 10 with cubic Beziers? Or even with line segments - which are simply Bezier curves of degree 1?

## Least squares reduction (degree N to M)

A Bezier curve can be reduced such that the total (for all 0<=t<=1) squared distance between the reduced Bezier and the original Bezier is a minimum.
Each of the `M` control points is a linear combination of the original `N` control points, and is usually best described by an M by N matrix. For every
combination of M and N, with M<N, there is a standard matrix.

The least squares reduction, though, is frequently not what is required; usually the *endpoints* need to be kept unchanged, so that a Bezier that has
been split up into subsections that *can* be sensibly reduced - and which join up - remain after reduction as a Bezier curves that join up. It would be
fairly silly to split up for drwaing a Bezier curve into segments that can be represented by lines, but then to choose lines that do not join up end-to-end.

## Reduction keeping endpoints - equispaced points

It is possible to reduce a Bezier and preserve the endpoints.

Perhaps the simplest method is to use the 'construct a Bezier from N points each at a different `t` value' described above: to reduce a Bezier to degree M
will require M+1 points, which can be taken at equal spacing along the original Bezier curve (e.g. at t=0, t=1/M, t=2/M, up to t=1).

Given that the new Bezier curve is guaranteed to go through the points given, which includes the endpoints (for t=0 and t=1), the new Bezier curve
will pass through the endpoints.

It is worth noting that the same algorithm can be used to elevate a Bezier curve (with M being larger than the original Bezier degree of N) - this is equivalent to the standard
Bezier elevation process (it must be, as both processes produce a Bezier curve that goes through the same `M+1` points, and there is only one such curve of degree M that does)

Of course, the points used do not have to be equispaced; they could be chosen to be more dense towards the endpoints, or however is desired.

## Error in reduction

The reduction of a Bezier curve (`B[N](t)`) of degree `N` to degree `M` (`B[N](t)`) (M<N) is lossy; that is, the two
curves do not nescessarily describe the same set of points.

# Testing

There is a significant test suite for this library; code coverage is analyzed using `cargo tarpaulin --out Html`

## Polynomials

## Beziers

## Metrics

## Approxmiation

## Building

# References

A Primer on Bézier Curves ("Pomax")

<https://pomax.github.io/bezierinfo/>

Adaptive Bezier Degree Reduction and Splitting
for Computationally Efficient Motion Planning (Omur Arslan and Aron Tiemessen):

<https://arxiv.org/pdf/2201.07834>

!*/

/*a Imports
*/
mod builder;
pub(crate) mod constants;
pub(crate) mod utils;

mod traits;
pub use traits::{Float, Num};

mod approximation;
mod distance;
pub(crate) mod polynomial;

mod bezier_iter;

mod implementations;

/// Useful functions for generating Bernstein coefficients and operating on Berstein Bezier curves
pub mod bernstein_fns;

pub mod bignum;

pub mod metrics;

/*a Exports
*/
pub use approximation::Approximation;
pub use bezier_iter::{
    BezierLineIter, BezierLineTIter, BezierPointIter, BezierPointTIter, BezierQuadIter,
    BezierSplitIter, BezierSplitTIter,
};
pub use traits::{
    BasicBezier, BezierConstruct, BezierDistance, BezierElevate, BezierEval, BezierIntoIterator,
    BezierMinMax, BezierOps, BezierReduce, BezierSection, BezierSplit, BoxedBezier,
};

pub use builder::BezierBuilder;
pub use constants::*;
pub use distance::BezierDistance2D;
pub use implementations::Bezier;
pub use polynomial::{PolyFindRoots, PolyNewtonRaphson, Polynomial};

pub use implementations::bezier_nd::BezierND;
