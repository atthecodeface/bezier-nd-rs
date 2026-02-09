//a Imports
use std::collections::HashMap;

use bezier_nd::{Float, PolyFindRoots, PolyNewtonRaphson, Polynomial};

mod utils;

fn test_polynomial_roots<F: Float + num_traits::AsPrimitive<isize>, P>(poly: &P, roots: &[F])
where
    P: PolyFindRoots<F>
        + PolyNewtonRaphson<F>
        + Polynomial<F>
        + std::fmt::Debug
        + AsMut<[F]>
        + Clone,
{
    let mut roots_to_find = HashMap::new();
    for r in roots {
        let r: isize = r.as_();
        if let Some(n) = roots_to_find.get_mut(&r) {
            *n += 1;
        } else {
            roots_to_find.insert(r, 1_usize);
        }
    }
    let roots_to_find = &mut roots_to_find;

    let mut result = [F::ZERO; 5];
    result[0] = F::ONE;

    match poly.degree() {
        1 => {
            let r0 = poly.find_roots_linear();
            utils::assert_near_equal_sorted_scale(&[r0.unwrap()], &[roots[0]], 1.0);
            let (r0, r1) = poly.find_roots_quad();
            assert!(r1.is_none());
            utils::assert_near_equal_sorted_scale(&[r0.unwrap()], &[roots[0]], 1.0);
            let (r0, r1, r2) = poly.find_roots_cubic();
            assert!(r1.is_none());
            assert!(r2.is_none());
            utils::assert_near_equal_sorted_scale(&[r0.unwrap()], &[roots[0]], 1.0);
        }
        2 => {
            let (r0, r1) = poly.find_roots_quad();
            if roots_to_find.is_empty() {
                assert!(r0.is_none());
                assert!(r1.is_none());
            } else {
                utils::assert_near_equal_sorted_scale(
                    &[r0.unwrap(), r1.unwrap()],
                    &[roots[0], roots[1]],
                    1.0,
                );
                let (r0, r1, r2) = poly.find_roots_cubic();
                assert!(r2.is_none());
                utils::assert_near_equal_sorted_scale(
                    &[r0.unwrap(), r1.unwrap()],
                    &[roots[0], roots[1]],
                    1.0,
                );
            }
        }
        3 => {
            let (r0, r1, r2) = poly.find_roots_cubic();
            if roots_to_find.len() == 1 {
            } else {
                utils::assert_near_equal_sorted_scale(
                    &[r0.unwrap(), r1.unwrap(), r2.unwrap()],
                    &[roots[0], roots[1], roots[2]],
                    1.0,
                );
            }
        }
        _ => {}
    }

    fn found_root<F: Float + num_traits::AsPrimitive<isize>>(
        root_set: &mut HashMap<isize, usize>,
        root: F,
    ) {
        eprintln!("Found root {root} from {root_set:?}");
        let root_i = root.round();
        let root_f = (root_i - root).abs();
        let root_i: isize = root_i.as_();
        let Some(root_count) = root_set.get(&root_i) else {
            panic!(
                "Root set {root_set:?} does not contain the root {root} {root_i} that was found"
            );
        };
        if *root_count > 1 {
            assert!(
                root_f < 0.1_f32.into(),
                "Repeated root {root} not close enough to an integer"
            );
        } else {
            assert!(
                root_f < 1E-5_f32.into(),
                "Root {root} not close enough to an integer"
            );
        }
        if *root_count == 1 {
            root_set.remove(&root_i);
        } else {
            *(root_set.get_mut(&root_i).unwrap()) = root_count - 1;
        }
    }
    loop {
        let mut poly_remaining = poly.clone();
        let divide_without_remainder =
            poly_remaining.set_divide(poly.clone().as_mut(), &result, 0.0001_f32.into());
        assert!(divide_without_remainder, "Should divide without remainder");
        poly_remaining.normalize(1E-6_f32.into());
        let d = poly_remaining.degree();
        eprintln!("{poly_remaining:?}: {result:?}, {poly:?}, {d}");
        match d {
            0 => {
                assert!(
                    roots_to_find.values().all(|a| *a == 0),
                    "Polynomial is a constant but there are more roots to find {roots_to_find:?}"
                );
                assert!((poly_remaining.as_mut()[0] - F::ONE).abs()<1E-3_f32.into(),
                    "Remaining polynomial when no roots are left should be constant one but is {poly_remaining:?}");
                break;
            }
            1 => {
                let Some(root) = poly_remaining.find_roots_linear() else {
                    panic!(
                        "Failed to find a linear root when was expecting one of {roots_to_find:?}"
                    );
                };
                found_root(roots_to_find, root);
                result.set_multiply(&result.clone(), &[-root.round(), F::ONE]);
            }
            2 => {
                let (r0, r1) = poly_remaining.find_roots_quad();
                if roots_to_find.is_empty() {
                    assert!(r0.is_none());
                    assert!(r1.is_none());
                    break;
                } else {
                    assert!(
                        r0.is_some() && r1.is_some(),
                        "Expected two roots to the remaining quadratic"
                    );
                    let r0 = r0.unwrap();
                    let r1 = r1.unwrap();
                    found_root(roots_to_find, r0);
                    found_root(roots_to_find, r1);
                    result.set_multiply(&result.clone(), &[-r0, F::ONE]);
                    result.set_multiply(&result.clone(), &[-r1, F::ONE]);
                }
            }
            _ => {
                let Some(root) = poly_remaining.find_root_nr(10.0_f32.into(), 1E-5_f32.into())
                else {
                    panic!("Failed to find root for polynomial {poly_remaining:?}");
                };
                assert!(
                    poly_remaining.calc(root).abs() < 0.1_f32.into(),
                    "Expected polynomial to be roughly 0.0 at root but was {}",
                    poly_remaining.calc(root)
                );
                result.set_multiply(&result.clone(), &[-root.round(), F::ONE]);
                found_root(roots_to_find, root);
            }
        }
    }
}
fn test_polynomial(roots: &[isize]) {
    assert!(roots.len() <= 4, "Can only test up to four roots");
    let mut poly = [0.0_f32; 5];
    assert_eq!(poly.degree(), 0);
    poly.differentiate();
    assert_eq!(&poly, &[0.0_f32; 5]);

    poly[0] = 1.0;
    assert_eq!(poly.degree(), 0);

    let mut roots_vec = Vec::new();
    let mut roots_to_find = HashMap::new();
    for r in roots {
        poly.set_multiply(&poly.clone(), &[(-*r) as f32, 1.]);
        roots_vec.push(*r as f32);
        if let Some(n) = roots_to_find.get_mut(r) {
            *n += 1;
        } else {
            roots_to_find.insert(*r, 1_usize);
        }
    }
    test_polynomial_roots(&poly, &roots_vec);
}

#[test]
fn test_differentiate() {
    let mut poly: [f32; 0] = [];
    poly.differentiate();

    let mut poly = [0.0_f32; 5];
    poly[0] = 1.0;

    poly.set_multiply(&poly.clone(), &[2.0, 3.0]);
    poly.set_multiply(&poly.clone(), &[1.0, 4.0]);
    // poly = (3x+2)*(4x+1) = 12x^2 + 11x + 2 -> 24x + 11
    let g = poly.gradient(0.3);
    poly.differentiate();
    (&mut poly[0..0]).differentiate();
    assert_eq!(poly[0], 11.0);
    assert_eq!(poly[1], 24.0);
    assert_eq!(poly[2], 0.0);
    poly.normalize(0.001_f32.into());
    assert_eq!(poly.degree(), 1);
    let g2 = poly.calc(0.3);
    assert_eq!(g, g2);

    let mut poly = [0.0_f32; 5];
    poly[0] = 1.0;
    poly.set_multiply(&poly.clone(), &[2.0, 3.0]);
    poly.set_multiply(&poly.clone(), &[1.0, 4.0]);
    // poly = (3x+2)*(4x+1) = 12x^2 + 11x + 2 -> 24x + 11
    let poly2 = poly.as_mut();
    let g = poly2.gradient(0.3);
    poly2.differentiate();
    assert_eq!(poly2[0], 11.0);
    assert_eq!(poly2[1], 24.0);
    assert_eq!(poly2[2], 0.0);
    poly2.normalize(0.001_f32.into());
    assert_eq!(poly2.degree(), 1);
    let g2 = poly2.calc(0.3);
    assert_eq!(g, g2);
}

#[test]
fn test_multiply() {
    let mut poly = [0.0_f32; 5];
    poly[0] = 1.0;
    poly.set_multiply(&poly.clone(), &[2.0, 3.0]);
    assert_eq!(poly[0], 2.0);
    assert_eq!(poly[1], 3.0);
    assert_eq!(poly[2], 0.0);

    poly.set_multiply(&poly.clone(), &[1.0, 4.0]);
    assert_eq!(poly[0], 2.0);
    assert_eq!(poly[1], 11.0);
    assert_eq!(poly[2], 12.0);
    assert_eq!(poly[3], 0.0);

    let mut poly2 = [0.0_f32; 5];
    poly2.as_mut().set_multiply(&poly, &[0.0, 2.0]);
    utils::assert_near_equal_scale(&poly[0..4], &poly2[1..5], 2.0);
}

#[test]
fn test_divide() {
    let mut poly = [0.0_f32; 5];
    poly[0] = 1.0;
    poly.set_multiply(&poly.clone(), &[2.0, 3.0]);
    poly.set_multiply(&poly.clone(), &[1.0, 4.0]);

    let mut poly_b = poly.clone();
    assert!(poly_b.set_divide(&mut poly.clone(), &[1.0, 4.0], 0.0001_f32.into()));
    assert_eq!(poly_b[0], 2.0);
    assert_eq!(poly_b[1], 3.0);
    assert_eq!(poly_b[2], 0.0);

    let mut poly_b = poly.clone();
    assert!(poly_b.set_divide(&mut poly.clone(), &[2.0, 3.0], 0.0001_f32.into()));
    assert_eq!(poly_b[0], 1.0);
    assert_eq!(poly_b[1], 4.0);
    assert_eq!(poly_b[2], 0.0);

    let mut poly_b = poly.clone();
    assert!(poly_b.set_divide(&mut poly.clone(), &poly.clone(), 0.0001_f32.into()));
    assert_eq!(poly_b[0], 1.0);
    assert_eq!(poly_b[1], 0.0);
    assert_eq!(poly_b[2], 0.0);

    let poly: [f32; 5] = [0., 0., 0., 1., 0.];
    let mut poly_b = poly.clone();
    assert!(poly_b.set_divide(&mut poly.clone(), &[0., 1., 0.], 0.0001_f32.into()));
    assert_eq!(poly_b[0], 0.0);
    assert_eq!(poly_b[1], 0.0);
    assert_eq!(poly_b[2], 1.0);
    assert_eq!(poly_b[3], 0.0);
    assert_eq!(poly_b[4], 0.0);

    let mut poly_c = poly_b.clone();
    assert!(poly_c.set_divide(&mut poly_b.clone(), &[0., 1., 0.], 0.0001_f32.into()));
    assert_eq!(poly_c[0], 0.0);
    assert_eq!(poly_c[1], 1.0);
    assert_eq!(poly_c[2], 0.0);
    assert_eq!(poly_c[3], 0.0);
    assert_eq!(poly_c[4], 0.0);

    let mut poly2 = [0.0_f32; 5];
    poly2
        .as_mut()
        .set_divide(poly.clone().as_mut(), &[0., 1., 0.], 0.0001_f32.into());
    utils::assert_near_equal_scale(&poly[1..5], &poly2[0..4], 1.0);
}

#[test]
fn test_roots() {
    test_polynomial(&[1]);
    test_polynomial(&[1, 2]);
    test_polynomial(&[6, 6]);
    test_polynomial(&[1, 2, 3]);
    test_polynomial(&[4, 4, 2]);
    test_polynomial(&[4, 4, 4]);
    test_polynomial(&[1, 2, 3, 4]);
    test_polynomial(&[4, 3, 2, 1]);

    test_polynomial_roots(&[1.], &[]);
    test_polynomial_roots(&[1., 0., 1.], &[]);
    test_polynomial_roots(&[1., 0., 0., 1.], &[-1.]);
}

fn test_specific_quad(poly: &[f32; 3]) {
    eprintln!("Test specific quad {poly:?}");
    let disc: f32 = poly[1] * poly[1] - 4.0 * poly[0] * poly[2];
    let x0 = (-poly[1] + disc.sqrt()) / 2.0 / poly[2];
    let x1 = (-poly[1] - disc.sqrt()) / 2.0 / poly[2];
    assert!(
        poly.calc(x0).abs() < 1E-4,
        "Expected roots to be x0 {x0} & x1 {x1} but poly[x0] = {}",
        poly.calc(x0)
    );
    assert!(
        poly.calc(x1).abs() < 1E-4,
        "Expected roots to be x0 {x0} & x1 {x1} but poly[x1] = {}",
        poly.calc(x1)
    );

    let root0 = poly.find_root_nr(0.0, 1E-6).unwrap();
    assert!(
        (root0 - x1).abs() < 1E-6 || (root0 - x0).abs() < 1E-6,
        "Newton-Raphson should have found root x0 {x0} or x1 {x1} but found {root0}"
    );

    let (r0, r1) = poly.find_roots_quad();
    eprintln!("Poly {poly:?} found quad roots {r0:?} {r1:?}, known to be {x0}, {x1}");
    let r0 = r0.expect("Two roots should have been found, first was None");
    let r1 = r1.expect("Two roots should have been found, first was None");
    utils::assert_near_equal_sorted_scale(&[r0, r1], &[x1, x0], 1.0);
}

fn test_specific_cubic(poly: &[f32; 4]) {
    eprintln!("Test specific cubic {poly:?}");
    let root0 = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
        .into_iter()
        .find_map(|x| poly.find_root_nr(x, 1E-6))
        .unwrap();
    let (r0, r1, r2) = poly.find_roots_cubic();
    eprintln!(
        "Poly {poly:?} found cubic roots {r0:?} {r1:?} {r2:?}, known to be (from NR) {root0},..."
    );
    let r0 = r0.expect("Cubics have at least one real root");
    let v = poly.calc(r0);
    assert!(
        v.abs() < 1E-5,
        "Root of specific cubic {poly:?} should have poly value 0 but got {v}"
    );
    //     let r1 = r1.expect("Two roots should have been found, first was None");
    // utils::assert_near_equal_sorted_scale(&[r0, r1], &[x1, x0], 1.0);
}

#[test]
fn specific_tests() {
    test_specific_quad(&[12.0_f32, -25.0, 4.0]);
    test_specific_quad(&[11.903965, -25.258091, 3.5866723]);
    test_specific_cubic(&[9.146291, -4.4568596, -97.76549, 89.3183]);
    test_specific_cubic(&[-0.024026299, 0.32613635, -1.094574, 1.0]);
    test_specific_cubic(&[-0.038705796, 0.38514474, -1.1151307, 1.0]);
    test_specific_cubic(&[-0.26665005, 1.2261595, -1.9053915, 1.0]);
    //    assert!(false);
}
