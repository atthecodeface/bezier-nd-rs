//a Documentation
/*!

Polynomial of best fit

given data xi, yi, we want E = Sum((yi-Sum(aj.xi^j))^2) to be a minimum

E = Sum((yi-Sum(aj.xi^j))^2)

Now d/dx f(g(x)) = g'(x) . f'(g(x))

  g(x) = yi-Sum(aj.xi^j)
  f(z) = z^2
  f'(z) = 2z

dE/daj = Sum( d/daj((yi-Sum(aj.xi^j))) . 2(yi-Sum(ak.xi^k)) )
       = Sum( 2(yi-Sum(ak.xi^k)) . (-xi^j) )

e.g.
Sum((yi - a - b.xi)^2) = Sum(yi*2 +a^2 +b^2.xi^2 - 2a.yi - 2.b.xi.yi + 2a.b.xi)
d/da(E) = Sum(2a - 2yi +2b.xi)
        = Sum(-2(yi -a - b.xi))
d/db(E) = Sum(2b.xi^2 - 2xi.yi +2a.xi)
        = Sum(-2xi.(yi - a - b.xi ))

dE/daj = Sum( 2(yi-Sum(aj.xi^j)) . (-xi^j) )
       = 0 for all aj at the minimum square error
Hence
Sum( 2(yi-Sum(ak.xi^k)) . (-xi^j) ) = 0 for all j
Sum( xi^j.yi ) = Sum(xi^j.Sum(ak.xi^k)) for all j

i.e. Xt.y = Xt.(X.a) (where Xt is X transpose)
or  (Xt.X).a = Xt.y
or         a = (Xt.X)' . Xt.y (where M' = inverse of M)

!*/

//a Imports
use crate::Num;

fn poly_calc<F: Num>(poly: &[F], x: F) -> F {
    let mut r = F::ZERO;
    let mut xn = F::ONE;
    for p in poly.iter() {
        r += (*p) * xn;
        xn *= x;
    }
    r
}

fn poly_gradient<F: Num>(poly: &[F], x: F) -> F {
    let mut r = F::ZERO;
    let mut xn = F::ONE;
    for (i, p) in poly.iter().enumerate().skip(1) {
        r += (*p) * xn * F::of_usize(i);
        xn *= x;
    }
    r
}

fn poly_d2f<F: Num>(poly: &[F], x: F) -> F {
    let mut r = F::ZERO;
    let mut xn = F::ONE;
    for (i, p) in poly.iter().enumerate().skip(2) {
        r += (*p) * xn * F::of_usize(((i - 1) * i));
        xn *= x;
    }
    r
}

fn poly_differentiate<F: Num>(poly: &mut [F]) {
    let n = poly.len();
    if n < 1 {
        return;
    }
    for i in 0..n - 1 {
        poly[i] = poly[i + 1] * F::of_usize(i + 1);
    }
    poly[n - 1] = F::ZERO;
}

fn poly_degree<F: Num>(poly: &[F]) -> usize {
    for (i, v) in poly.iter().enumerate().rev() {
        if !v.is_zero() {
            return i;
        }
    }
    0
}

fn poly_normalize<F: Num>(poly: &mut [F], eps: F) {
    for v in poly.iter_mut() {
        if (*v > -eps) && (*v < eps) {
            *v = F::ZERO;
        }
    }
}

/// Multiply poly by multiplicand and store the result
///
/// If the multiplier has degree 0 then it is a constant; two polys of degree 0
/// produce a polynomial which has a degree of 0, and hence an array length of 1
///
/// If the multiplier has degree 1 then it is linear; two polys of degree 1
/// produce a polynomial which has a degree of 2, and hence an array length of 3
#[track_caller]
fn poly_multiply<F: Num>(poly: &[F], multiplicand: &[F], result: &mut [F]) {
    let mul_deg = poly_degree(multiplicand);
    let poly_deg = poly_degree(poly);
    for p in result.iter_mut() {
        *p = F::ZERO;
    }
    assert!(
        result.len() > mul_deg + poly_deg,
        "Result of multiplication must be large enough to store result"
    );
    for (m_i, m) in multiplicand.iter().take(mul_deg + 1).enumerate() {
        for (p_i, p) in poly.iter().take(poly_deg + 1).enumerate() {
            result[m_i + p_i] += (*m) * (*p);
        }
    }
}

// Return true if the remeainder is all less than or equal to eps
#[track_caller]
fn poly_divide<F: Num>(poly: &mut [F], divisor: &[F], result: &mut [F], eps: F) -> bool {
    let div_deg = poly_degree(divisor);
    let poly_deg = poly_degree(poly);
    for p in result.iter_mut() {
        *p = F::ZERO;
    }
    if poly_deg >= div_deg {
        let n_iter = poly_deg + 1 - div_deg;
        assert!(
            result.len() >= n_iter,
            "Result of division must be large enough to store dividend"
        );
        let div_max = divisor[div_deg];
        for i in 0..n_iter {
            let p_max = poly[poly_deg - i];
            // Subtract p_max/div_max * divisor from poly
            let amount = p_max / div_max;
            result[n_iter - i - 1] = amount;
            for j in 0..(div_deg + 1) {
                poly[poly_deg - i - j] -= amount * divisor[div_deg - j];
            }
        }
    }
    poly.iter().all(|v| *v > -eps && *v < eps)
}

//tt Polynomial
/// A collections of methods for polynomials
pub trait Polynomial<F: Num> {
    /// Differentiate the polynomial
    fn differentiate(&mut self);
    /// Get the degree of the polynomial (largest non-zero coefficient)
    fn degree(&self) -> usize;
    /// Normalize - zero any coefficient whose absolute value is less than epsilon
    fn normalize(&mut self, epsilon: F);
    /// Calculate the value at x
    fn calc(&self, x: F) -> F;
    /// Calculate the gradient at x
    fn gradient(&self, x: F) -> F;
    /// Calculate the second derivative at x
    fn d2f(&self, x: F) -> F;
    /// Multiply by a polynomial
    fn set_multiply(&mut self, poly: &[F], multiplicand: &[F]);
    /// Set to result of polynomial divided by divisor, leaving the remainder in poly
    ///
    /// Return false if the polynomial division left a remainder polynomial with any
    /// coefficient greater than 'eps'
    fn set_divide(&mut self, poly: &mut [F], divisor: &[F], eps: F) -> bool;
}

//ip Polynomial for [F; N]
impl<F: Num, const N: usize> Polynomial<F> for [F; N] {
    fn degree(&self) -> usize {
        poly_degree(self)
    }
    fn normalize(&mut self, epsilon: F) {
        poly_normalize(self, epsilon)
    }
    fn calc(&self, x: F) -> F {
        poly_calc(self, x)
    }
    fn gradient(&self, x: F) -> F {
        poly_gradient(self, x)
    }
    fn d2f(&self, x: F) -> F {
        poly_d2f(self, x)
    }
    fn differentiate(&mut self) {
        poly_differentiate(self);
    }
    /// Multiply by a polynomial
    #[track_caller]
    fn set_multiply(&mut self, poly: &[F], multiplicand: &[F]) {
        poly_multiply(poly, multiplicand, self);
    }
    /// Set to result of polynomial divided by divisor, leaving the remainder in poly
    #[track_caller]
    fn set_divide(&mut self, poly: &mut [F], divisor: &[F], eps: F) -> bool {
        poly_divide(poly, divisor, self, eps)
    }
}

//ip Polynomial for [F; N]
impl<F: Num> Polynomial<F> for [F] {
    fn degree(&self) -> usize {
        poly_degree(self)
    }
    fn normalize(&mut self, epsilon: F) {
        poly_normalize(self, epsilon)
    }
    fn calc(&self, x: F) -> F {
        poly_calc(self, x)
    }
    fn gradient(&self, x: F) -> F {
        poly_gradient(self, x)
    }
    fn d2f(&self, x: F) -> F {
        poly_d2f(self, x)
    }
    fn differentiate(&mut self) {
        poly_differentiate(self);
    }
    /// Multiply by a polynomial
    fn set_multiply(&mut self, poly: &[F], multiplicand: &[F]) {
        poly_multiply(poly, multiplicand, self);
    }
    /// Set to result of polynomial divided by divisor, leaving the remainder in poly
    fn set_divide(&mut self, poly: &mut [F], divisor: &[F], eps: F) -> bool {
        poly_divide(poly, divisor, self, eps)
    }
}

/// Improve an estimate for a root; as we get closer the dx should get smaller
///
/// If the dx is *larger* than the last dx then stop
fn improve_root<F: Num>(poly: &[F], x: F) -> Option<(F, F)> {
    let f = poly.calc(x);
    let df = poly.gradient(x);
    let d2f = poly.d2f(x);
    let denom = df * df - f * d2f * F::frac(1, 2);
    if denom.is_unreliable_divisor() {
        if df.is_unreliable_divisor() {
            None
        } else {
            let new_x = x - f / df;
            Some((new_x, (new_x - x).nabs()))
        }
    } else {
        let new_x = x - f * df / denom;
        Some((new_x, (new_x - x).nabs()))
    }
}

pub fn find_real_roots_linear<F: Num>(poly: &[F]) -> Option<F> {
    assert!(
        poly.len() >= 2,
        "Root of a linear polynomial requires at least two coefficients (and [2..] should be zero)"
    );
    if poly[1].is_unreliable_divisor() {
        None
    } else {
        Some(-poly[0] / poly[1])
    }
}

pub fn find_real_roots_quad<F: Num>(poly: &[F]) -> (Option<F>, Option<F>) {
    assert!(
        poly.len() >= 3,
        "Root of a quadratic polynomial requires at least three coefficients (and [3..] should be zero)"
    );
    if poly[2].is_unreliable_divisor() {
        (find_real_roots_linear(poly), None)
    } else {
        let two_a = poly[2] + poly[2];
        let b = poly[1];
        let c = poly[0];
        let disc = b * b - F::of_i32(2) * two_a * c;
        if disc < F::zero() {
            (None, None)
        } else {
            let disc_sq = disc.sqrt_est();
            (Some((-b + disc_sq) / two_a), Some((-b - disc_sq) / two_a))
        }
    }
}

pub fn find_real_roots_cubic<F: Num + num_traits::Float>(
    poly: &[F],
) -> (Option<F>, Option<F>, Option<F>) {
    assert!(
        poly.len() >= 4,
        "Root of a cubic polynomial requires at least four coefficients (and [4..] should be zero)"
    );
    if poly[3].is_unreliable_divisor() {
        let (root_a, root_b) = find_real_roots_quad(poly);
        (root_a, root_b, None)
    } else {
        let a = poly[3];
        let b = poly[2];
        let c = poly[1];
        let d = poly[0];

        let two = F::of_usize(2);
        let three = F::of_usize(3);
        let sqrt3 = three.sqrt();

        let delta_0 = b * b - a * c * three;
        let delta_1 = b * b * b * two - a * b * c * F::of_usize(9) + a * a * d * F::of_usize(27);
        let disc = delta_1 * delta_1 - delta_0 * delta_0 * delta_0 * F::of_usize(4);
        // dbg!(delta_0, delta_1, disc);
        if delta_0.is_unreliable_divisor() {
            // three identical roots if delta_1 is zero
            //
            // If delta_1 is nonzero then one real root (?)
            let big_c = (delta_1 / two).cbrt();
            let x = (b + big_c) / (-a * three);
            if delta_1.is_unreliable_divisor() {
                (Some(x), Some(x), Some(x))
            } else {
                (Some(x), None, None)
            }
        } else if disc > F::zero() {
            // discriminant is +ve, so square root of discriminant is real
            let big_c = (({
                if delta_1 < F::zero() {
                    delta_1 - disc.sqrt()
                } else {
                    delta_1 + disc.sqrt()
                }
            }) / two)
                .cbrt();
            let x = (b + big_c + delta_0 / big_c) / (-a * three);
            (Some(x), None, None)
        } else {
            let mut r0 = None;
            let mut r1 = None;
            let mut r2 = None;
            let small = F::frac(1, 1000);

            // Note there is *always* one real root

            // Three real roots if 1 - delta_0/cbrt_mag^2 == 0 i.e. delta_0 == cbrt_mag^2
            //
            // i.e. delta_0 == |big_c_cubed|.cbrt()
            //
            // discriminant is -ve, so square root of discriminant is imaginary
            let big_c_cubed_i = (-disc).sqrt() / two;
            let big_c_cubed_r = delta_1 / two;
            // dbg!(big_c_cubed_r, big_c_cubed_i);
            let cbrt_theta = big_c_cubed_i.atan2(big_c_cubed_r) / three;
            let cbrt_mag =
                (big_c_cubed_i * big_c_cubed_i + big_c_cubed_r * big_c_cubed_r).powf(F::frac(1, 6));
            let big_c_r = cbrt_theta.cos() * cbrt_mag;
            let big_c_i = cbrt_theta.sin() * cbrt_mag;
            // Note that / big_C is the same as * big_C comp / |C|^2
            let thing = big_c_i - delta_0 * big_c_i / (cbrt_mag * cbrt_mag);
            // dbg!(thing);
            if thing.nabs() < small {
                let x = (b + big_c_r + delta_0 * big_c_r / (cbrt_mag * cbrt_mag)) / (-a * three);
                // eprintln!("Value at {x} {}", poly.calc(x));
                r0 = Some(x);
            }
            let new_big_c_r = (-big_c_r - big_c_i * sqrt3) / two;
            let new_big_c_i = (-big_c_i + big_c_r * sqrt3) / two;
            let big_c_r = new_big_c_r;
            let big_c_i = new_big_c_i;
            let thing = big_c_i - delta_0 * big_c_i / (cbrt_mag * cbrt_mag);
            // dbg!(thing);
            if thing.nabs() < small {
                let x = (b + big_c_r + delta_0 * big_c_r / (cbrt_mag * cbrt_mag)) / (-a * three);
                // eprintln!("Value at {x} {}", poly.calc(x));
                r1 = Some(x);
            }
            let new_big_c_r = (-big_c_r - big_c_i * sqrt3) / two;
            let new_big_c_i = (-big_c_i + big_c_r * sqrt3) / two;
            let big_c_r = new_big_c_r;
            let big_c_i = new_big_c_i;
            let thing = big_c_i - delta_0 * big_c_i / (cbrt_mag * cbrt_mag);
            // dbg!(thing);
            if thing.nabs() < small {
                let x = (b + big_c_r + delta_0 * big_c_r / (cbrt_mag * cbrt_mag)) / (-a * three);
                // eprintln!("Value at {x} {}", poly.calc(x));
                r2 = Some(x);
            }
            (r0, r1, r2)
        }
    }
}

//a PolyFindRoots
//tt PolyFindRoots
/// Types that provide PolyFindRoots will provide accurate calculation of
/// polynomial roots if the polynomial they represent is linear, quadratic
/// or cubic
pub trait PolyFindRoots<F: Num + num_traits::Float> {
    /// Assume the polynomial is linear, and solve
    fn find_roots_linear(&self) -> Option<F>;
    /// Assume the polynomial is quadratic, and solve
    fn find_roots_quad(&self) -> (Option<F>, Option<F>);
    /// Assume the polynomial is cubic, and solve
    fn find_roots_cubic(&self) -> (Option<F>, Option<F>, Option<F>);
}

/// A type that provide PolyNewtonRaphson allow for finding
/// roots using Newton-Raphson; this does not require sqrt or cbrt,
/// and so only requires 'Num' not 'Float'
pub trait PolyNewtonRaphson<F: Num> {
    /// Improve a root using Newton-Raphson
    fn improve_root(&self, x: F) -> Option<(F, F)>;
    /// Find a root using Newton-Raphson given a starting guess, and minimum gradient (in case of root multiplicity)
    fn find_root_nr_with_err(&self, mut x: F, min_dx: F, mut iters: usize) -> (F, F)
    where
        Self: Polynomial<F>,
    {
        let mut last_dx = F::of_i32(i32::MAX);
        while let Some((improved_x, improved_dx)) = self.improve_root(x) {
            // eprintln!("x:{x}, {improved_x} {improved_dx}");
            x = improved_x;
            last_dx = improved_dx;
            if iters == 0 || improved_dx < min_dx {
                break;
            }
            iters -= 1;
        }
        (x, last_dx)
    }

    /// Find a root using Newton-Raphson given a starting guess, and minimum gradient (in case of root multiplicity)
    fn find_root_nr(&self, x: F, min_dx: F) -> Option<F>
    where
        Self: Polynomial<F>,
    {
        let (x, last_dx) = self.find_root_nr_with_err(x, min_dx, 100);
        if last_dx > min_dx {
            None
        } else {
            Some(x)
        }
    }
}

//ip PolyFindRoots for [F; N]
impl<F: Num + num_traits::Float, const N: usize> PolyFindRoots<F> for [F; N] {
    fn find_roots_linear(&self) -> Option<F> {
        find_real_roots_linear(self.as_slice())
    }
    fn find_roots_quad(&self) -> (Option<F>, Option<F>) {
        find_real_roots_quad(self.as_slice())
    }
    fn find_roots_cubic(&self) -> (Option<F>, Option<F>, Option<F>) {
        find_real_roots_cubic(self.as_slice())
    }
}
impl<F: Num, const N: usize> PolyNewtonRaphson<F> for [F; N] {
    fn improve_root(&self, x: F) -> Option<(F, F)> {
        improve_root(self, x)
    }
}
