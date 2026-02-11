mod utils;
use bezier_nd::{bernstein_fns, Num};
use utils::test_beziers_approx_eq;

fn elevate_by_one_matrix<F: Num>(degree: usize) -> Vec<F> {
    // Elevate by one for degree N is and N+2 by N+1 matrix
    let n = degree + 1;
    let n1 = degree + 2;
    let mut elevate_by_one = vec![F::ZERO; n * n1];
    let f =
        bernstein_fns::elevate_matrix::generate_elevate_by_one_matrix(&mut elevate_by_one, degree);
    // eprintln!("{f} &{elevate_by_one:?},");
    elevate_by_one
}

#[test]
fn elevate_constants() {
    let degree = 1;
    let n2 = degree + 2;
    let n1 = degree + 1;
    let e = elevate_by_one_matrix::<f32>(degree);
    let mut e_t = e.clone();
    geo_nd::matrix::transpose_dyn(n2, n1, &e, &mut e_t);
    let mut e_e_t = vec![0.0; n1 * n1];
    geo_nd::matrix::multiply_dyn(n1, n2, n1, &e, &e_t, &mut e_e_t);
    let mut e_e_t_inv = e_e_t.clone();
    let mut lu = e_e_t.clone();
    let mut pivot = vec![0; n1];
    let mut tr1 = vec![0.0; n1];
    let mut tr2 = vec![0.0; n1];
    let det = geo_nd::matrix::lup_decompose(n1, &e_e_t, &mut lu, &mut pivot);
    eprintln!("{det} {e_e_t:?}");
    assert!(geo_nd::matrix::lup_invert(
        n1,
        &lu,
        &pivot,
        &mut e_e_t_inv,
        &mut tr1,
        &mut tr2
    ));
    eprintln!("{e_e_t:?} {e_e_t_inv:?}");
    let mut r = vec![0.0; n1 * n2];
    geo_nd::matrix::multiply_dyn(n1, n1, n2, &e_e_t_inv, &e_t, &mut r);
    eprintln!("{r:?}");
    assert!(false);
}
