//a Imports
use std::collections::HashSet;

use bezier_nd::{PolyFindRoots, Polynomial};

fn test_polynomial(roots: &[isize]) {
    assert!(roots.len() <= 4, "Can only test up to four roots");
    let mut poly = [0.0_f32; 5];
    let mut result = [0.0_f32; 5];
    poly[0] = 1.0;
    result[0] = 1.0;

    let mut roots_to_find = HashSet::new();
    for r in roots {
        poly.set_multiply(&poly.clone(), &[(-*r) as f32, 1.]);
        roots_to_find.insert(*r);
    }

    fn found_root(root_set: &mut HashSet<isize>, root: f32) {
        eprintln!("Found root {root} from {root_set:?}");
        let root_i = root.round();
        let root_f = (root_i - root).abs();
        let root_i = root_i as isize;
        assert!(root_f < 1E-3, "Root {root} not close enough to an integer");
        assert!(
            root_set.contains(&root_i),
            "Root set {root_set:?} does not contain the root {root_i} that was found"
        );
        root_set.remove(&root_i);
    }
    loop {
        let mut poly_remaining = poly.clone();
        dbg!(&poly_remaining, &result, &poly);
        assert!(
            poly_remaining.set_divide(&mut poly.clone(), &result, (0.0001_f32.into())),
            "Should divide without remainder"
        );
        poly_remaining.normalize(1E-6);
        let d = poly_remaining.degree();
        dbg!(&poly_remaining, &result, &poly, &d);
        match d {
            0 => {
                assert!(
                    roots_to_find.is_empty(),
                    "Polynomial is a constant but there are more roots to find {roots_to_find:?}"
                );
                assert!((poly_remaining[0] - 1.0).abs()<1E-3,
                    "Remaining polynomial when no roots are left should be constant one but is {poly_remaining:?}");
                break;
            }
            1 => {
                let Some(root) = poly_remaining.find_roots_linear() else {
                    panic!(
                        "Failed to find a linear root when was expecting one of {roots_to_find:?}"
                    );
                };
                found_root(&mut roots_to_find, root);
                result.set_multiply(&result.clone(), &[-root, 1.0]);
            }
            2 => {
                let (Some(r0), Some(r1)) = poly_remaining.find_roots_quad() else {
                    panic!("Failed to find real roots from a quadratic when expecting {roots_to_find:?}");
                };
                dbg!(r0, r1);
                found_root(&mut roots_to_find, r0);
                found_root(&mut roots_to_find, r1);
                result.set_multiply(&result.clone(), &[-r0, 1.0]);
                result.set_multiply(&result.clone(), &[-r1, 1.0]);
            }
            _ => {
                let Some(root) = poly_remaining.find_root_nr(100.0, 1E-7) else {
                    panic!("Failed to find root for polynomial {poly_remaining:?}");
                };
                result.set_multiply(&result.clone(), &[-root, 1.0]);
                found_root(&mut roots_to_find, root);
            }
        }
    }
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
}

#[test]
fn test_divide() {
    let mut poly = [0.0_f32; 5];
    poly[0] = 1.0;
    poly.set_multiply(&poly.clone(), &[2.0, 3.0]);
    poly.set_multiply(&poly.clone(), &[1.0, 4.0]);

    let mut poly_b = poly.clone();
    assert!(poly_b.set_divide(&mut poly.clone(), &[1.0, 4.0], (0.0001_f32.into())));
    assert_eq!(poly_b[0], 2.0);
    assert_eq!(poly_b[1], 3.0);
    assert_eq!(poly_b[2], 0.0);

    let mut poly_b = poly.clone();
    assert!(poly_b.set_divide(&mut poly.clone(), &[2.0, 3.0], (0.0001_f32.into())));
    assert_eq!(poly_b[0], 1.0);
    assert_eq!(poly_b[1], 4.0);
    assert_eq!(poly_b[2], 0.0);

    let mut poly_b = poly.clone();
    assert!(poly_b.set_divide(&mut poly.clone(), &poly.clone(), (0.0001_f32.into())));
    assert_eq!(poly_b[0], 1.0);
    assert_eq!(poly_b[1], 0.0);
    assert_eq!(poly_b[2], 0.0);
}

#[test]
fn test_roots() {
    test_polynomial(&[1]);
    test_polynomial(&[1, 2]);
    test_polynomial(&[1, 2, 3]);
    test_polynomial(&[1, 2, 3, 4]);
    test_polynomial(&[4, 3, 2, 1]);
}
