mod utils;
use bezier_nd::bernstein_fns::elevate_reduce_matrix::generate_elevate_by_one_matrix;
use bezier_nd::bignum::RationalN;
use bezier_nd::{
    BezierElevate, BezierEval, BezierMetric, BezierND, BezierReduce, BezierReduction, Num,
};
use geo_nd::{matrix, vector};
use utils::*;

#[test]
fn bernstein_matrix_br() {
    type N = RationalN<8>;
    let mut bern_n = [N::default(); 100];
    generate_bernstein_matrix_br(&mut bern_n[0..3], 2, &[0_i64.into()]);
    assert_eq!(bern_n[0], 1_i64.into());
    assert_eq!(bern_n[1], 0_i64.into());
    assert_eq!(bern_n[2], 0_i64.into());
    generate_bernstein_matrix_br(&mut bern_n[0..3], 2, &[(1_i64, 2_u64).into()]);
    assert_eq!(bern_n[0], (1_i64, 4_u64).into());
    assert_eq!(bern_n[1], (1_i64, 2_u64).into());
    assert_eq!(bern_n[2], (1_i64, 4_u64).into());
    generate_bernstein_matrix_br(&mut bern_n[0..3], 2, &[1_i64.into()]);
    assert_eq!(bern_n[0], 0_i64.into());
    assert_eq!(bern_n[1], 0_i64.into());
    assert_eq!(bern_n[2], 1_i64.into());
}

/// Generate a Bernstein reduction matrix *to* degree (using ts provided)
///
/// ts must be of length degree+1, and the reduction will keep the values
/// of the reduced Bezier at these values of t *unchanged* compared to the
/// original Bezier
fn generate_bs_reduce_matrix<N: Num + Into<f64>>(
    reduce_from_degree: usize,
    degree: usize,
    ts: &[N],
    reduce: &mut [N],
    elevated_reduce: &mut [N],
) {
    assert_eq!(
        ts.len(),
        degree + 1,
        "Generated a reduction matrix to {degree} requires {} points",
        degree + 1
    );

    let n2 = (degree + 1) * (degree + 1);

    // bern_n is a (degree+1) * (degree+1) matrix
    let mut bern_n = vec![N::ZERO; n2];
    let mut bern_np1 = vec![N::ZERO; (degree + 1) * (reduce_from_degree + 1)];
    let mut tr0 = vec![N::ZERO; degree + 1];
    let mut tr1 = vec![N::ZERO; degree + 1];
    let mut pivot = vec![0; degree + 1];
    let mut lu = bern_n.clone();

    // bern_n is the matrix that can be applied to a Bezier of 'degree' to
    // generate the points at the values of t given in 'ts'
    //
    // bern_n is (degree+1) * (degree+1)
    generate_bernstein_matrix_br(&mut bern_n[0..n2], degree, &ts);

    // bern_np1 is the matrix that can be applied to a Bezier of 'reduce_from_degree' to
    // generate the points at the values of t given in 'ts'
    //
    // bern_np1 is (degree+1) * (reduce_from_degree+1)
    generate_bernstein_matrix_br(
        &mut bern_np1[0..(reduce_from_degree + 1) * (degree + 1)],
        reduce_from_degree,
        &ts,
    );

    // bern_n_inverse is such that bern_n_inverse * bern_n == Identity
    //
    // bern_n_inverse is (degree+1) * (degree+1)
    let mut bern_n_inverse = bern_n.clone();
    assert_ne!(
        matrix::lup_decompose(degree + 1, &bern_n, &mut lu, &mut pivot),
        N::ZERO,
        "Matrix is invertible"
    );
    assert!(
        matrix::lup_invert(
            degree + 1,
            &lu,
            &pivot,
            &mut bern_n_inverse,
            &mut tr0,
            &mut tr1
        ),
        "Matrix is invertible"
    );

    // reduce = bern_n_inverse * bern_np1
    //
    // (degree+1) * (degree+1) times (degree+1) * (reduce_from_degree+1) yields...
    //
    // reduce is (degree+1) *(reduce_from_degree+1)
    matrix::multiply_dyn(
        degree + 1,
        degree + 1,
        reduce_from_degree + 1,
        &bern_n_inverse,
        &bern_np1,
        reduce,
    );

    // Find reduction of elevation
    let mut test = [N::zero(); 100];
    let elevate = generate_elevate_by_n_matrix(degree, reduce_from_degree - degree);

    // Validate that reduce of elevation is the identity (approx)
    matrix::multiply_dyn(
        degree + 1,
        reduce_from_degree + 1,
        degree + 1,
        &reduce,
        &elevate,
        &mut test,
    );
    eprintln!(
        "test {:?}",
        test.iter()
            .take((degree + 1) * (degree + 1))
            .map(|s| s.to_string())
            .collect::<Vec<_>>()
    );
    assert_near_identity(degree + 1, &test);

    // Generate the elevated_reduce matrix
    //
    matrix::multiply_dyn(
        reduce_from_degree + 1,
        degree + 1,
        reduce_from_degree + 1,
        &elevate,
        &reduce,
        elevated_reduce,
    );

    eprintln!("Bernstein reduce for ts from {reduce_from_degree} to degree {degree}");
    eprintln!(
        "  &[{}]",
        reduce
            .iter()
            .take((degree + 1) * (reduce_from_degree + 1))
            .fold(String::new(), |string, s| string
                + &(<N as Into<f64>>::into(*s)).to_string()
                + ", ")
    );
    eprintln!(
        "  &[{}]",
        reduce
            .iter()
            .take((degree + 1) * (reduce_from_degree + 1))
            .fold(String::new(), |string, s| string + &s.to_string() + ", ")
    );
    eprintln!("Bernstein elevated reduce for ts {reduce_from_degree} to degree {degree}");
    eprintln!(
        "  &[{}]",
        elevated_reduce
            .iter()
            .take((reduce_from_degree + 1) * (reduce_from_degree + 1))
            .fold(String::new(), |string, s| string
                + &(<N as Into<f64>>::into(*s)).to_string()
                + ", ")
    );
    eprintln!(
        "  &[{}]",
        elevated_reduce
            .iter()
            .take((reduce_from_degree + 1) * (reduce_from_degree + 1))
            .fold(String::new(), |string, s| string + &s.to_string() + ", ")
    );
}

#[test]
fn validate_bs_reductions_to_quad() {
    let mut reduce = [RationalN::<8>::default(); 200];
    let mut elevated_reduce = [RationalN::<8>::default(); 200];

    for degree in [1, 2, 3] {
        let ts: Vec<RationalN<8>> = (0..(degree + 1))
            .map(|i| (i as i64, degree).into())
            .collect();

        for reduce_from_degree in (degree + 1)..11 {
            generate_bs_reduce_matrix(
                reduce_from_degree as usize,
                degree as usize,
                &ts,
                &mut reduce,
                &mut elevated_reduce,
            );
            let (i, _lcm) = RationalN::<8>::with_common_denom(reduce.iter());
            eprintln!(
                "\n\n****\nreduce = [{}]",
                i.fold(String::new(), |s, a| s + &a.to_string() + ", ")
            );
            let (i, lcm) = RationalN::<8>::with_common_denom(elevated_reduce.iter());
            eprintln!(
                "elevated_reduce = {lcm} : [{}]",
                i.fold(String::new(), |s, a| s + &a.to_string() + ", ")
            );
        }
    }
    // assert!(false);
}

#[test]
fn bernstein_reduce_matrix_c_to_q() {
    let mut bern_n = [0.0f64; 100];
    let mut bern_np1 = [0.0f64; 100];

    // Find reduction from cubic (degree 3, 4 pts) to quadratic using ts of 0, 1/2, 1
    let ts = [0.0, 0.5, 1.0];
    generate_bernstein_matrix(&mut bern_n[0..9], 2, &ts);
    assert_near_equal(&bern_n, &[1.0, 0., 0., 0.25, 0.5, 0.25, 0., 0., 1.]);
    generate_bernstein_matrix(&mut bern_np1[0..12], 3, &ts);
    assert_near_equal(
        &bern_np1,
        &[1.0, 0., 0., 0., 0.125, 0.375, 0.375, 0.125, 0., 0., 0., 1.],
    );
    let mut bern_n_inverse = bern_n.clone();
    assert!(
        matrix::invert::<_, 9, 3>(&mut bern_n_inverse.first_chunk_mut::<9>().unwrap()),
        "Bernstein basis must be invertible"
    );
    assert_near_equal(&bern_n_inverse, &[1.0, 0., 0., -0.5, 2.0, -0.5, 0., 0., 1.]);

    // Reduction matrix is 3x4
    let reduce = matrix::multiply::<_, 9, 12, 12, 3, 3, 4>(
        bern_n_inverse.first_chunk::<9>().unwrap(),
        bern_np1.first_chunk::<12>().unwrap(),
    );
    eprintln!("bernstein_reduce[3] = {reduce:?}");
    assert_near_equal(
        &reduce,
        &[
            1.0, 0.0, 0., 0., -0.25, 0.75, 0.75, -0.25, 0.0, 0.0, 0., 1.0,
        ],
    );

    // Find elevation from quadratic to cubic
    let (scale, mut elevate) = generate_elevate_by_one_matrix(2);
    for e in elevate.iter_mut() {
        *e /= scale;
    }

    // Find reduction of elevation
    let reduced_elevate = matrix::multiply::<_, 12, 12, 9, 3, 4, 3>(
        reduce.first_chunk::<12>().unwrap(),
        elevate.first_chunk::<12>().unwrap(),
    );
    dbg!(&reduced_elevate);
    assert_near_equal(&reduced_elevate, &matrix::identity3());
}

#[test]
fn reduce_and_elevate_cubic() {
    let p0 = [0., 0.];
    let p1 = [10., -1.];
    let p2 = [6., 1.];
    let p3 = [20., 5.];
    let max_dc = 0.03;

    let b = BezierND::<f64, 4, 2>::new(&[p0, p1, p2, p3]);
    assert_eq!(b.degree(), 3);
    let b_split = b.reduce_and_split_iter(
        &bezier_nd::REDUCE_BY_ONE_UNIFORM[1],
        &bezier_nd::ER_UNIFORM_MINUS_I[1],
        2,
        max_dc,
    );
    let mut t = 0.0;
    for (ns, bs) in b_split {
        let t0 = t;
        let dt = 1.0 / (2.0_f64).powi(ns as i32);
        let t1 = t + dt;
        let p0 = b.point_at(t0);
        let p1 = b.point_at(t1);
        let ps0 = bs.point_at(0.0);
        let ps1 = bs.point_at(1.0);
        assert!(
            vector::distance_sq(&p0, &ps0) < 1E-4,
            "Split should match at pt0 {t0} {p0:?} {ps0:?}"
        );
        assert!(
            vector::distance_sq(&p1, &ps1) < 1E-4,
            "Split should match at pt1 {t1} {p1:?} {ps1:?}"
        );
        for tp in float_iter(100) {
            let p = b.point_at(t0 + tp * dt);
            let ps = bs.point_at(tp);
            assert!(
                vector::distance_sq(&p, &ps) < max_dc,
                "Split should be no more than {max_dc} apartat pt {tp} of {t0}->{t1} {p:?} {ps:?}"
            );
        }
        t = t1;
    }
}

#[track_caller]
fn test_reduce_least_squares<F, const D: usize, B>(bezier: &B)
where
    B: BezierReduce<F, [F; D]>,
    F: Num,
{
    assert!(
        bezier.can_reduce(BezierReduction::LeastSquares),
        "Must be able to reduce by least squares"
    );
    let reduced = bezier
        .reduce(BezierReduction::LeastSquares)
        .expect("Must be able to reduce as it committed to it");
    let reduced_dc_sq = bezier.dc_sq_from_reduction(BezierReduction::LeastSquares);
    let r_vec: Vec<_> = reduced.control_points().into();
    let e = r_vec.elevate_by_one().unwrap();
    let dc_sq = bezier
        .metric_from(Some(e.control_points()), BezierMetric::MaxControlSquared)
        .unwrap();
    eprintln!("{reduced_dc_sq} {dc_sq}");
    let d_sq = bezier
        .metric_from(
            Some(e.control_points()),
            BezierMetric::SumDistanceSquared(1000),
        )
        .unwrap();
    let m_sq = bezier
        .metric_from(
            Some(e.control_points()),
            BezierMetric::MaxDistanceSquared(1000),
        )
        .unwrap();
    for i in 0..r_vec.len() {
        let mut mr = r_vec.clone();
        if vector::length_sq(&mr[i]).is_unreliable_divisor() {
            continue;
        }
        mr[i] = vector::scale(mr[i], 1.01_f32.into());
        let me = mr.elevate_by_one().unwrap();
        let md_sq = bezier
            .metric_from(
                Some(me.control_points()),
                BezierMetric::SumDistanceSquared(1000),
            )
            .unwrap();
        let mm_sq = bezier
            .metric_from(
                Some(me.control_points()),
                BezierMetric::MaxDistanceSquared(1000),
            )
            .unwrap();
        assert!(md_sq*1.01_f32.into() > d_sq, "Expected tweaking to least-squared reduction to have a larger total distance {md_sq} than the actual least squared reduction {d_sq}");
        assert!(mm_sq*1.01_f32.into() > m_sq, "Expected tweaking to least-squared reduction to have a larger max distance {mm_sq} than the actual least squared reduction {m_sq}");
    }
}

#[track_caller]
fn test_reduce_preserve_uniform<F, const D: usize, B>(bezier: &B)
where
    B: BezierReduce<F, [F; D]> + BezierEval<F, [F; D]> + std::fmt::Debug,
    <B as BezierReduce<F, [F; D]>>::Reduced: std::fmt::Debug,
    F: Num,
{
    assert!(
        bezier.can_reduce(BezierReduction::UniformT),
        "Must be able to reduce by preserving uniform ts"
    );
    let reduced = bezier
        .reduce(BezierReduction::UniformT)
        .expect("Must be able to reduce as it committed to it");
    let reduced_dc_sq = bezier.dc_sq_from_reduction(BezierReduction::UniformT);

    let r_vec: Vec<_> = reduced.control_points().into();
    let e = r_vec.elevate_by_one().unwrap();
    let dc_sq = bezier
        .metric_from(Some(e.control_points()), BezierMetric::MaxControlSquared)
        .unwrap();
    eprintln!("Bezier is {bezier:?}");
    eprintln!("Reduction is {reduced:?}");
    eprintln!("Elevated reduction is {e:?}");
    eprintln!("{reduced_dc_sq} {dc_sq}");
    for t in utils::float_iter(r_vec.len()) {
        assert_near_equal(&bezier.point_at(t), &e.point_at(t));
    }
}

#[test]
fn least_squares() {
    let mut rng = utils::make_random("test_cubics_seed");
    let distribution = rand::distr::Uniform::new(-10.0_f32, 10.0).unwrap();

    // Cannot reduce linear beziers - nor points
    for _ in 0..10 {
        let mut pts = [[0.0_f32; 2]; 3];
        utils::set_random_bezier(&mut rng, &distribution, &mut pts);
        test_reduce_least_squares(&pts);

        let mut pts = [[0.0_f32; 2]; 4];
        utils::set_random_bezier(&mut rng, &distribution, &mut pts);
        test_reduce_least_squares(&pts);

        let mut pts = vec![[0.0_f32; 2]; 3];
        utils::set_random_bezier(&mut rng, &distribution, &mut pts);
        test_reduce_least_squares(&pts);

        let mut pts = vec![[0.0_f32; 2]; 4];
        utils::set_random_bezier(&mut rng, &distribution, &mut pts);
        test_reduce_least_squares(&pts);

        let mut pts = vec![[0.0_f32; 2]; 5];
        utils::set_random_bezier(&mut rng, &distribution, &mut pts);
        test_reduce_least_squares(&pts);

        let mut pts = vec![[0.0_f32; 2]; 6];
        utils::set_random_bezier(&mut rng, &distribution, &mut pts);
        test_reduce_least_squares(&pts);

        let mut pts = vec![[0.0_f32; 2]; 10];
        utils::set_random_bezier(&mut rng, &distribution, &mut pts);
        test_reduce_least_squares(&pts);
    }
}

#[test]
fn preserve_uniform() {
    let mut rng = utils::make_random("test_preserve_uniform");
    let distribution = rand::distr::Uniform::new(-10.0_f32, 10.0).unwrap();

    // Cannot reduce linear beziers - nor points
    for _ in 0..10 {
        let mut pts = [[0.0_f32; 2]; 3];
        utils::set_random_bezier(&mut rng, &distribution, &mut pts);
        test_reduce_preserve_uniform(&pts);

        let mut pts = [[0.0_f32; 2]; 4];
        utils::set_random_bezier(&mut rng, &distribution, &mut pts);
        test_reduce_preserve_uniform(&pts);

        let mut pts = vec![[0.0_f32; 2]; 3];
        utils::set_random_bezier(&mut rng, &distribution, &mut pts);
        test_reduce_preserve_uniform(&pts);

        let mut pts = vec![[0.0_f32; 2]; 4];
        utils::set_random_bezier(&mut rng, &distribution, &mut pts);
        test_reduce_preserve_uniform(&pts);

        let mut pts = vec![[0.0_f32; 2]; 5];
        utils::set_random_bezier(&mut rng, &distribution, &mut pts);
        test_reduce_preserve_uniform(&pts);

        let mut pts = vec![[0.0_f32; 2]; 6];
        utils::set_random_bezier(&mut rng, &distribution, &mut pts);
        test_reduce_preserve_uniform(&pts);

        let mut pts = vec![[0.0_f32; 2]; 10];
        utils::set_random_bezier(&mut rng, &distribution, &mut pts);
        test_reduce_preserve_uniform(&pts);
    }
}
