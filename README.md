# bezier-nd

A Bezier curve class supporting linear, quadratic and cubic Bezier curves,
using an arbitrary point class.

Example uses would be for 2-dimensional Bezier curves whose
coordinates are `[f64; 2]`, or for 3-dimensional Bezier curves using
coordinates of `[f32; 3]`.

The Bezier curve supports bisection, and then splitting into straight
lines within a given `straightness` bound; iterators are provided to
automatically trace a Bezier as lines or points within such a
straightness, for rendering puroses.

The Bezier type also supports rounding of corners and circular arc
generation, utilizing a very accurate function for any angle of
rounding up to 180 degrees, derived from a curve-fit from experimental
data, rather than an explicit mathematical function for the control
point generation (the standard analytical approach).

This crate is in beta; it is used in a small number of applications,
and the functionality is mature; the API is stable, but may be enhanced.

## Usage

Add this to your `Cargo.toml`:

```toml
[dependencies]
bezier-nd = "0.1.4"
```

## Releases

Release notes are available in [RELEASES.md](RELEASES.md).

## License

Licensed under either of

 * [Apache License, Version 2.0](http://www.apache.org/licenses/LICENSE-2.0)
 * [MIT license](http://opensource.org/licenses/MIT)

at your option.

### Contribution

Unless you explicitly state otherwise, any contribution intentionally submitted
for inclusion in the work by you, as defined in the Apache-2.0 license, shall be
dual licensed as above, without any additional terms or conditions.
