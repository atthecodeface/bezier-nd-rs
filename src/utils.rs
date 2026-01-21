use geo_nd::vector;
use geo_nd::Float;

// Calculate the length squared of a line, the projection of the point onto the line
/// and if the line is approximately zero length
pub fn relative_to_line<F: Float, const D: usize>(
    pt: &[F; D],
    l0: &[F; D],
    l1: &[F; D],
) -> (F, F, bool) {
    // Make everything relative to line.0
    //
    // so point is P and line is L (i.e. on the line = k*L for some 0 <= k <= 1)
    //
    // Note P = kL + N for some k and normal to the line N (N.L = 0)
    //
    // P.L = k(L.L) + N.L = k|L|^2
    let l1_m_l0 = vector::sub(*l1, l0, F::one());
    // len_l1_m_l0_sq = |L|^2
    let len_l1_m_l0_sq = vector::length_sq(&l1_m_l0);

    let pt_m_l0 = vector::sub(*pt, l0, F::one());

    if len_l1_m_l0_sq < (f32::EPSILON).into() {
        return (F::ZERO, len_l1_m_l0_sq, false);
    }

    // pt_along_line = P . L = k|L|^2
    let pt_along_line = vector::dot(&pt_m_l0, &l1_m_l0);
    (pt_along_line / len_l1_m_l0_sq, len_l1_m_l0_sq, true)
}
