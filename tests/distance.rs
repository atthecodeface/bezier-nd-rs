//a Imports
use bezier_nd::Bezier;
use geo_nd::FArray;

// Test the distance-from-line function for points that are on the line
//
// The result for all of these points should be zero
#[test]
fn test_line_distance_on_line() {
    let line = ([0., 0.], [1.0, 0.0]);
    assert_eq!(
        Bezier::<f32, 2>::pt_distance_sq_from(&line.0, &line),
        0.0,
        "End of a line has zero distance"
    );
    assert_eq!(
        Bezier::<f32, 2>::pt_distance_sq_from(&line.1, &line),
        0.0,
        "End of a line has zero distance"
    );
    assert_eq!(
        Bezier::<f32, 2>::pt_distance_sq_from(&[0.5, 0.], &line),
        0.0,
        "Middle of a line has zero distance"
    );
    assert_eq!(
        Bezier::<f32, 2>::pt_distance_sq_from(&[0.1, 0.], &line),
        0.0,
        "Point along a line has zero distance"
    );

    // Second line - diagonal with y=-x
    let line = ([2., -2.], [4.0, -4.0]);
    assert_eq!(
        Bezier::<f32, 2>::pt_distance_sq_from(&line.0, &line),
        0.0,
        "End of a line has zero distance"
    );
    assert_eq!(
        Bezier::<f32, 2>::pt_distance_sq_from(&line.1, &line),
        0.0,
        "End of a line has zero distance"
    );
    assert_eq!(
        Bezier::<f32, 2>::pt_distance_sq_from(&[3., -3.], &line),
        0.0,
        "Middle of a line has zero distance"
    );
    assert_eq!(
        Bezier::<f32, 2>::pt_distance_sq_from(&[2.25, -2.25], &line),
        0.0,
        "Point along a line has zero distance"
    );
}

// Test the distance-from-line function for points that are beyond the ends of the line
//
// The result for all of these points should be the distance to the closest endpoint
#[test]
fn test_line_distance_beyond_ends_of_line() {
    let line = ([0., 0.], [1.0, 0.0]);
    assert_eq!(
        Bezier::<f32, 2>::pt_distance_sq_from(&[-1., 0.], &line),
        1.0,
        "Beyond the end of the line, distance from endpoint"
    );
    assert_eq!(
        Bezier::<f32, 2>::pt_distance_sq_from(&[2., 0.], &line),
        1.0,
        "Beyond the end of the line, distance from endpoint"
    );
    assert_eq!(
        Bezier::<f32, 2>::pt_distance_sq_from(&[-2., 0.], &line),
        4.0,
        "Beyond the end of the line, distance from endpoint"
    );
    assert_eq!(
        Bezier::<f32, 2>::pt_distance_sq_from(&[3., 0.], &line),
        4.0,
        "Beyond the end of the line, distance from endpoint"
    );

    let line = ([2., -2.], [4.0, -4.0]);
    assert_eq!(
        Bezier::<f32, 2>::pt_distance_sq_from(&[0., 0.], &line),
        8.0,
        "Beyond the end of the line, distance from endpoint"
    );
    assert_eq!(
        Bezier::<f32, 2>::pt_distance_sq_from(&[1., -1.], &line),
        2.0,
        "Beyond the end of the line, distance from endpoint"
    );
    assert_eq!(
        Bezier::<f32, 2>::pt_distance_sq_from(&[5., -5.], &line),
        2.0,
        "Beyond the end of the line, distance from endpoint"
    );
    assert_eq!(
        Bezier::<f32, 2>::pt_distance_sq_from(&[6., -6.], &line),
        8.0,
        "Beyond the end of the line, distance from endpoint"
    );
}

// Test the distance-from-line function for points that are off the line
//
// The result for all of these points should be the distance to the point projected to the line *or* to the closest endpoint
#[test]
fn test_line_distance_not_on_line() {
    let line = ([0., 0.], [1.0, 0.0]);
    for (i, (pt, d)) in [
        ([-1., 1.], 2.),
        ([1., -1.], 1.),
        ([-1., -1.], 2.),
        ([1., 1.], 1.),
        ([0.5, 0.5], 0.25),
        ([0.5, 1.0], 1.0),
        ([0.5, -0.5], 0.25),
        ([0.5, -1.0], 1.0),
        ([0.1, 0.5], 0.25),
        ([0.2, 1.0], 1.0),
        ([0.3, -0.5], 0.25),
        ([0.4, -1.0], 1.0),
    ]
    .iter()
    .enumerate()
    {
        assert_eq!(
            Bezier::<f32, 2>::pt_distance_sq_from(pt, &line),
            *d,
            "Point {i} {pt:?} mismatched distance from line"
        );
    }
}

#[test]
fn test_distance_bezier_line() {
    let p0: FArray<f32, 2> = [0., 0.].into();
    let p1: FArray<f32, 2> = [10., 0.].into();
    let b = Bezier::line(&p0, &p1);
    assert!(
        b.min_distance_sq_from(&p0) == 0.0,
        "Distance to points is zero"
    );
    assert!(
        b.min_distance_sq_from(&p1) == 0.0,
        "Distance to points is zero"
    );
    assert_eq!(
        b.min_distance_sq_from(&[5.0, 0.0]),
        0.0,
        "Distance to point on the line is zero"
    );
}

// Test the minimum distance squared for a quad Bezier with control point on both sides of the line
#[test]
fn test_distance_bezier_quad() {
    let p0: FArray<f32, 2> = [0., 0.].into();
    let p1: FArray<f32, 2> = [4., 0.].into();
    let c: FArray<f32, 2> = [2., 2.].into();
    let b = Bezier::quadratic(&p0, &c, &p1);
    assert!(
        b.min_distance_sq_from(&p0) == 0.0,
        "Distance to points is zero"
    );
    assert!(
        b.min_distance_sq_from(&p1) == 0.0,
        "Distance to points is zero"
    );
    assert_eq!(
        b.min_distance_sq_from(&c),
        0.0,
        "Distance to control point is zero"
    );
    for (i, pt) in [
        ([1., 1.]),
        ([1., 0.5]),
        ([2., 0.]),
        ([2., 1.]),
        ([2., 1.5]),
        ([3., 1.]),
    ]
    .iter()
    .enumerate()
    {
        assert_eq!(
            b.min_distance_sq_from(pt),
            0.0,
            "Distance to point {i} {pt:?} in the triangle is zero"
        );
    }

    for (i, (pt, d)) in [
        ([0., 2.], 2.),
        ([2., 4.], 4.),
        ([4., 2.], 2.),
        ([2., -2.], 4.),
    ]
    .iter()
    .enumerate()
    {
        assert_eq!(
            b.min_distance_sq_from(pt),
            *d,
            "Distance to point {i} {pt:?} outside the triangle"
        );
    }
}
