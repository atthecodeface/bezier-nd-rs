mod utils;
use bezier_nd::{bernstein_fns, constants};

/// Generate a Bernstein reduction matrix *to* degree (using ts provided)
///
/// ts must be of length degree+1, and the reduction will keep the values
/// of the reduced Bezier at these values of t *unchanged* compared to the
/// original Bezier

const MAX_DEGREE: usize = 10;

#[test]
fn elevate() {
    type F = f64;
    let display_matrix = utils::display::string_of_f64_matrix_f64;
    let accuracy = 9E-14_f64;
    let mut elevate_by_one = vec![];
    for degree in 1..MAX_DEGREE {
        let n2 = degree + 2;
        let n1 = degree + 1;
        let (scale, mut elevate) =
            bernstein_fns::elevate_reduce_matrix::generate_elevate_by_one_matrix::<F>(degree);
        for e in elevate.iter_mut() {
            *e /= scale;
        }
        elevate_by_one.push(elevate);
    }

    let mut elevate_by_one_string = String::new();
    elevate_by_one_string += "pub const ELEVAT_BY_ONE: &[&[f64]] = &[\n";
    for r in &elevate_by_one {
        elevate_by_one_string += "&";
        elevate_by_one_string += &display_matrix(r);
        elevate_by_one_string += ",\n";
    }
    elevate_by_one_string += "];\n";
    eprintln!("{elevate_by_one_string}");

    assert!(constants::ELEVATE_BY_ONE.len() == elevate_by_one.len());
    for (i, (r, c)) in elevate_by_one
        .iter()
        .zip(constants::ELEVATE_BY_ONE.iter())
        .enumerate()
    {
        assert_eq!(r.len(), c.len());
        for (j, (r, c)) in r.iter().zip(c.iter()).enumerate() {
            let r = f64::try_from(*r).unwrap();
            let diff = r - *c;
            assert!(
                diff.abs() < accuracy,
                "ELEVATE_BY_ONE {i}:{j} {r} {c} {diff}"
            );
        }
    }
}

#[test]
fn reduce_uniform() {
    // type F = RationalN<16>;
    // let display_matrix = utils::display::string_of_f64_matrix;
    // let accuracy = 1E-18_f64;
    type F = f64;
    let display_matrix = utils::display::string_of_f64_matrix_f64;
    let accuracy = 9E-14_f64;

    let mut reduce_by_one = vec![];
    let mut elevated_reduce_by_one = vec![];
    for degree in 1..MAX_DEGREE {
        let n2 = degree + 2;
        let n1 = degree + 1;

        // r is n1 by n2
        // e is n2 by n1
        eprintln!("Reduce degree {n1} to {degree}");
        let r =
            bernstein_fns::elevate_reduce_matrix::reduce_uniform_matrix::<F>(degree + 1, degree);
        // utils::display::eprintln_rational_matrix(&r, n1);
        reduce_by_one.push(r.clone());

        // eprintln!("{degree}: elevate");
        let (scale, mut e) =
            bernstein_fns::elevate_reduce_matrix::generate_elevate_by_one_matrix::<F>(degree);
        for e in e.iter_mut() {
            *e /= scale;
        }
        eprintln!("Elevation-by-one matrix {degree}");
        utils::display::eprintln_matrix(&e, n1);

        // Reduction * Elevation * Point == Point
        //
        // r * e is n1 by n1; it should be the identity
        // eprintln!("{degree}: reduction of elevation - must be the identity");
        let mut test = vec![0.0_f32.into(); n1 * n1];
        geo_nd::matrix::multiply_dyn(n1, n2, n1, &r, &e, &mut test);
        eprintln!("Test should be identity r * e");
        utils::display::eprintln_matrix(&test, n1);
        utils::assert_near_identity(n1, &test);

        // Elevated reduction minus identity
        eprintln!("Elevated reduce of degree {n1} to {degree} minus identity");
        let mut test = vec![0.0_f32.into(); n2 * n2];
        geo_nd::matrix::multiply_dyn(n2, n1, n2, &e, &r, &mut test);
        for i in 0..n2 {
            test[i * (n2 + 1)] -= <f32 as Into<F>>::into(1_f32);
        }
        // utils::display::eprintln_rational_matrix(&test, degree + 2);
        elevated_reduce_by_one.push(test);
    }

    let mut reduce_by_one_string = String::new();
    reduce_by_one_string += "pub const REDUCE_BY_ONE_UNIFORM: &[&[f64]] = &[\n";
    for r in &reduce_by_one {
        reduce_by_one_string += "&";
        reduce_by_one_string += &display_matrix(r);
        reduce_by_one_string += ",\n";
    }
    reduce_by_one_string += "];\n";
    eprintln!("{reduce_by_one_string}");

    let mut elevated_reduce_by_one_string = String::new();
    elevated_reduce_by_one_string += "pub const ER_UNIFORM_MINUS_I: &[&[f64]] = &[\n";
    for r in &elevated_reduce_by_one {
        elevated_reduce_by_one_string += "&";
        elevated_reduce_by_one_string += &display_matrix(r);
        elevated_reduce_by_one_string += ",\n";
    }
    elevated_reduce_by_one_string += "];\n";
    eprintln!("{elevated_reduce_by_one_string}");

    assert!(constants::REDUCE_BY_ONE_UNIFORM.len() == reduce_by_one.len());
    for (i, (r, c)) in reduce_by_one
        .iter()
        .zip(constants::REDUCE_BY_ONE_UNIFORM.iter())
        .enumerate()
    {
        assert_eq!(r.len(), c.len());
        for (j, (r, c)) in r.iter().zip(c.iter()).enumerate() {
            let r = f64::try_from(*r).unwrap();
            let diff = r - *c;
            assert!(
                diff.abs() < accuracy,
                "REDUCE_BY_ONE_UNIFORM {i}:{j} {r} {c} {diff}"
            );
        }
    }

    assert!(constants::ER_UNIFORM_MINUS_I.len() == reduce_by_one.len());
    for (i, (r, c)) in elevated_reduce_by_one
        .iter()
        .zip(constants::ER_UNIFORM_MINUS_I.iter())
        .enumerate()
    {
        assert_eq!(r.len(), c.len());
        for (j, (r, c)) in r.iter().zip(c.iter()).enumerate() {
            let r = f64::try_from(*r).unwrap();
            let diff = r - *c;
            assert!(
                diff.abs() < accuracy,
                "ER_UNIFORM_MINUS_I {i}:{j} {r} {c} {diff}"
            );
        }
    }
}

#[test]
fn reduce_min_l2() {
    // type F = RationalN<16>;
    // let display_matrix = utils::display::string_of_f64_matrix;
    // let accuracy = 1E-18_f64;
    type F = f64;
    let display_matrix = utils::display::string_of_f64_matrix_f64;
    let accuracy = 1E-14_f64;

    let mut reduce_by_one = vec![];
    let mut elevated_reduce_by_one = vec![];
    for degree in 1..MAX_DEGREE {
        let n2 = degree + 2;
        let n1 = degree + 1;

        // r is n1 by n2
        // e is n2 by n1
        eprintln!("Reduce degree {n1} to {degree}");
        let r = bernstein_fns::elevate_reduce_matrix::reduce_l2_min_by_one_matrix::<F>(degree);
        // utils::display::eprintln_rational_matrix(&r, degree + 2);
        reduce_by_one.push(r.clone());

        // eprintln!("{degree}: elevate");
        let (scale, mut e) =
            bernstein_fns::elevate_reduce_matrix::generate_elevate_by_one_matrix::<F>(degree);
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
            test[i * (n2 + 1)] -= <f32 as Into<F>>::into(1_f32);
        }
        // utils::display::eprintln_rational_matrix(&test, degree + 2);
        elevated_reduce_by_one.push(test);
    }

    let mut reduce_by_one_string = String::new();
    reduce_by_one_string += "pub const REDUCE_BY_ONE_LSQ: &[&[f64]] = &[\n";
    for r in &reduce_by_one {
        reduce_by_one_string += "&";
        reduce_by_one_string += &display_matrix(r);
        reduce_by_one_string += ",\n";
    }
    reduce_by_one_string += "];\n";
    eprintln!("{reduce_by_one_string}");

    let mut elevated_reduce_by_one_string = String::new();
    elevated_reduce_by_one_string += "pub const ER_LSQ_MINUS_I: &[&[f64]] = &[\n";
    for r in &elevated_reduce_by_one {
        elevated_reduce_by_one_string += "&";
        elevated_reduce_by_one_string += &display_matrix(r);
        elevated_reduce_by_one_string += ",\n";
    }
    elevated_reduce_by_one_string += "];\n";
    eprintln!("{elevated_reduce_by_one_string}");

    assert!(constants::REDUCE_BY_ONE_LSQ.len() == reduce_by_one.len());
    for (i, (r, c)) in reduce_by_one
        .iter()
        .zip(constants::REDUCE_BY_ONE_LSQ.iter())
        .enumerate()
    {
        assert_eq!(r.len(), c.len());
        for (j, (r, c)) in r.iter().zip(c.iter()).enumerate() {
            let r = f64::try_from(*r).unwrap();
            let diff = r - *c;
            assert!(
                diff.abs() < accuracy,
                "REDUCE_BY_ONE_LSQ {i}:{j} {r} {c} {diff}"
            );
        }
    }

    for (i, (r, c)) in elevated_reduce_by_one
        .iter()
        .zip(constants::ER_LSQ_MINUS_I.iter())
        .enumerate()
    {
        assert_eq!(r.len(), c.len());
        for (j, (r, c)) in r.iter().zip(c.iter()).enumerate() {
            let r = f64::try_from(*r).unwrap();
            let diff = r - *c;
            assert!(
                diff.abs() < accuracy,
                "ER_LSQ_MINUS_I {i}:{j} {r} {c} {diff}"
            );
        }
    }
}
