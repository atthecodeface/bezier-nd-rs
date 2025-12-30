# Release 0.6.1 (2025-11-25)

- Removed requirement for FArray / Vector for types that are needed
  for the remaining methods

# Release 0.6.0 (2025-11-22)

- Changed the API to use [Float; D] rather than types that support
  geo_nd traits explicitly

- This means that users should no longer use type like
  `geo_nd::FArray<f32,2>` to create a Bezier, but should use [f32; 2]
  instead; and points returned by iterators etc would be [f32; 2] not
  the `geo_nd::FArray` type

- Users no longer need to depend on geo_nd for most Bezier methods
  (only if the Bezier type they use is exposed as a generic might this
  be an issue)

- `bezier_between`, `arc`, `of_round_corner`, and
  `center_radius_of_bezier_arc` do require some use of `geo_nd` FArray
  and Vector, and these are currently public re-exported by this
  library; this use (and export) will go away in the future.

**Contributors**: @atthecodeface

# Release 0.5.1 (2025-11-12)

- Migrated to geo-nd version 0.6; API should be the same

**Contributors**: @atthecodeface

# Release 0.5.0 (2023-02-18)

- Migrated to geo-nd version 0.5; firming up the API as it is now in use

**Contributors**: @atthecodeface

# Release 0.1.4 (2021-06-23)

- Added documentation for the library and some simple code examples

**Contributors**: @atthecodeface

# Release 0.1.3 (2021-06-23)

- Bezier arcs and rounded corners now use
   a better approach using a matched polynomial
   to derive 'lambda'; this is now tested in
   regression to a good degree of accuracy.

**Contributors**: @atthecodeface

# Release 0.1.0 (2021-06-22)

- Publishing on crates.io for the first time

**Contributors**: @atthecodeface

