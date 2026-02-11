//a Imports
mod utils;
use bezier_nd::Bezier;
use utils::test_beziers_approx_eq;

#[test]
fn bisect() {
    let mut rng = utils::make_random("test_bisect_seed");
    let distribution = rand::distr::Uniform::new(-10.0_f32, 10.0).unwrap();

    let pts: [[f32; 2]; 2] = utils::new_random_point_array(&mut rng, &distribution);
    let b: Bezier<_, _> = pts.into();
    test_beziers_approx_eq(&b, &pts);
    test_beziers_approx_eq(&pts, &pts);
    test_beziers_approx_eq(&pts.to_vec(), &pts);

    let pts: [[f32; 2]; 3] = utils::new_random_point_array(&mut rng, &distribution);
    let b: Bezier<_, _> = pts.into();
    test_beziers_approx_eq(&b, &pts);
    test_beziers_approx_eq(&pts, &pts);
    test_beziers_approx_eq(&pts.to_vec(), &pts);

    let pts: [[f32; 2]; 4] = utils::new_random_point_array(&mut rng, &distribution);
    let b: Bezier<_, _> = pts.into();
    test_beziers_approx_eq(&b, &pts);
    test_beziers_approx_eq(&pts, &pts);
    test_beziers_approx_eq(&pts.to_vec(), &pts);
}
