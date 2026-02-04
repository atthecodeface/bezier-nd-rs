mod utils;
use bezier_nd::{BezierConstruct, BezierEval, BezierND};
use utils::assert_near_equal;

use bezier_nd::BezierBuilder;

#[test]
fn build_line() {
    let p0 = [1., 2.];
    let p1 = [3., 0.];
    let p2 = [5., -2.];
    let dp = [1., 1.];
    let p3 = [2., 3.];

    let mut builder = BezierBuilder::default();
    builder.add_point_at(0.0_f32, p0);
    builder.add_point_at(0.5_f32, p1);

    let b = BezierND::<_, 10, _>::of_builder(&builder).unwrap();
    assert_near_equal(&b.point_at(0.0), &p0);
    assert_near_equal(&b.point_at(0.5), &p1);

    let b = <[[f32; 2]; 2]>::of_builder(&builder).unwrap();
    utils::vec_eq(&b[0], &p0);
    utils::vec_eq(&b[1], &p2);

    let b = <[[f32; 2]; 3]>::of_builder(&builder);
    assert!(b.is_err(), "Bezier is underconstrained for a quadratic");

    let b = <Vec<[f32; 2]>>::of_builder(&builder).unwrap();
    utils::vec_eq(&b[0], &p0);
    utils::vec_eq(&b[1], &p2);

    let mut builder = BezierBuilder::default();
    builder.add_point_at(0.0_f32, p0);
    builder.add_derivative_at(0.0_f32, 1, dp);
    let b = <Vec<[f32; 2]>>::of_builder(&builder).unwrap();
    utils::vec_eq(&b[0], &p0);
    utils::vec_eq(&b[1], &p3);
}

#[test]
fn build_quad() {
    let p0 = [1., 2.];
    let p1 = [3., 0.];
    let p2 = [10., 2.];
    let dp = [1., 1.];

    let mut builder = BezierBuilder::default();
    builder.add_point_at(0.0_f32, p0);
    builder.add_point_at(0.5_f32, p1);
    builder.add_point_at(1.0_f32, p2);

    let b = BezierND::<_, 10, _>::of_builder(&builder).unwrap();
    assert_near_equal(&b.point_at(0.0), &p0);
    assert_near_equal(&b.point_at(0.5), &p1);
    assert_near_equal(&b.point_at(1.0), &p2);

    let b = <[[f32; 2]; 2]>::of_builder(&builder);
    assert!(b.is_err(), "Bezier is overconstrained for a line");

    let b = <[[f32; 2]; 3]>::of_builder(&builder).unwrap();
    assert_near_equal(&b.point_at(0.0), &p0);
    assert_near_equal(&b.point_at(0.5), &p1);
    assert_near_equal(&b.point_at(1.0), &p2);

    let b = <[[f32; 2]; 4]>::of_builder(&builder);
    assert!(b.is_err(), "Bezier is underconstrained for a cubic");

    let b = <Vec<[f32; 2]>>::of_builder(&builder).unwrap();
    utils::vec_eq(&b[0], &p0);
    utils::vec_eq(&b[2], &p2);

    let mut builder = BezierBuilder::default();
    builder.add_point_at(0.0_f32, p0);
    builder.add_derivative_at(0.0_f32, 1, dp);
    builder.add_point_at(1.0_f32, p2);

    let b = <Vec<[f32; 2]>>::of_builder(&builder).unwrap();
    utils::vec_eq(&b[0], &p0);
    utils::vec_eq(&b[2], &p2);
    assert_near_equal(&b.point_at(0.0), &p0);
    assert_eq!(b.derivative_at(0.0).0, 2.0);
    let scale = b.derivative_at(0.0).0;
    assert_near_equal(&b.derivative_at(0.0).1, &[dp[0] * scale, dp[1] * scale]);
    assert_near_equal(&b.point_at(1.0), &p2);
}
