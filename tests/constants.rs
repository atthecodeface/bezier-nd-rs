mod utils;
use bezier_nd::{bernstein_fns, bignum::RationalN};

#[test]
fn reduce_min_l2() {
    for degree in 1..7 {
        let n2 = degree + 2;
        let n1 = degree + 1;

        // r is n1 by n2
        // e is n2 by n1
        eprintln!("Reduce degree {n1} to {degree}");
        let r = bernstein_fns::elevate_reduce_matrix::reduce_l2_min_by_one_matrix::<RationalN<4>>(
            degree,
        );
        utils::display::eprintln_rational_matrix(&r, degree + 2);

        // eprintln!("{degree}: elevate");
        let (scale, mut e) = bernstein_fns::elevate_reduce_matrix::generate_elevate_by_one_matrix::<
            RationalN<4>,
        >(degree);
        for e in e.iter_mut() {
            *e /= scale;
        }
        // utils::display::eprintln_matrix(&e, n1);

        // Reduction * Elevation * Point == Point
        //
        // r * e is n1 by n1; it should be the identity
        // eprintln!("{degree}: reduction of elevation - must be the identity");
        let mut test = vec![0.0_f32.into(); n1 * n1];
        geo_nd::matrix::multiply_dyn(n1, n2, n1, &r, &e, &mut test);
        // utils::display::eprintln_matrix(&test, n1);
        utils::assert_near_identity(n1, &test);

        // Elevated reduction minus identity
        eprintln!("Elevated reduce of degree {n1} to {degree} minus identity");
        let mut test = vec![0.0_f32.into(); n2 * n2];
        geo_nd::matrix::multiply_dyn(n2, n1, n2, &e, &r, &mut test);
        for i in 0..n2 {
            test[i * (n2 + 1)] -= 1_f32.into();
        }
        utils::display::eprintln_rational_matrix(&test, degree + 2);
    }
}
