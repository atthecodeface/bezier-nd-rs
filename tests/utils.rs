use geo_nd::Float;

pub fn float_iter<F: Float>(t0: F, t1: F, n: usize) -> impl Iterator<Item = F> {
    assert!(
        n >= 2,
        "Float iterator must have at least two steps for begin and end"
    );
    let r = t1 - t0;
    assert!(r >= F::zero(), "Float range must be positive");
    (0..n).map(move |i: usize| (<f32 as Into<F>>::into((i as f32))) / (r * ((n - 1) as f32).into()))
}

#[track_caller]
pub fn assert_near_identity<F: Float>(n: usize, m: &[F]) {
    let eps = (1E-4_f32).into();
    for i in 0..n * n {
        if i % (n + 1) == 0 {
            assert!((m[i] - F::one()).abs() < eps);
        } else {
            assert!(m[i].abs() < eps);
        }
    }
}

#[track_caller]
pub fn assert_near_equal<F: Float>(m0: &[F], m1: &[F]) {
    let eps = (1E-4_f32).into();

    for (i, (v0, v1)) in m0.iter().zip(m1.iter()).enumerate() {
        assert!(
            (*v0 - *v1).abs() < eps,
            "Data {i} is mismatch in {m0:?} <> {m1:?}"
        );
    }
}
