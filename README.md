# geo-nd

A collection of vector, matrix and quaternion types and traits for
Rust, particularly aimed at 2D, 3D and OpenGL/Vulkan applications,
where the elements use arrays and slices of floats.

The usage model in 3D graphics is through ownership of large arrays of
floats by the application, encompassing sets of vertices and matrices,
with the library enabling in-place operations on those vertices and
matrices.

The library takes advantage of the developiong
[core_simd](https://rust-lang.github.io/stdsimd/core_simd/) crate to
provide architecture-specific implementations, without attempting to
support architectural SIMD implementations within this crate.

This crate is in alpha; it is used in a small number of applications,
and the functionality is mature, but the API may undergo some changes
in the near future (through Q3 2021) to ensure high performance OpenGL
and Vulkan operation while maintaining simplicity of operation for
other applications.

## Usage

Add this to your `Cargo.toml`:

```toml
[dependencies]
geo-nd = "0.1"
```

## Features

*NOTE* Until *core_simd* is part of std or on crates.io the following does not work - the dependency in the geo-nd/Cargo.toml for core_simd must be commented out on crates.io; so a locally modified clone of this crate is required until that point.

This crate can be used with the upcoming [core_simd](https://rust-lang.github.io/stdsimd/core_simd/) crate by enabling
the default `simd` feature. Use this in `Cargo.toml`:

```toml
[dependencies.geo-nd]
version = "0.1"
features = ["simd"]
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
