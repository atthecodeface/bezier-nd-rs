use super::{Fixed, IsFixed};

trait TestKind<I>:
    num_traits::Num
    + num_traits::ConstOne
    + num_traits::ConstZero
    + std::fmt::Debug
    + num_traits::FromPrimitive
    + Ord
    + Copy
    + std::borrow::Borrow<I>
    + std::borrow::BorrowMut<I>
{
}
impl<T, I> TestKind<I> for T where
    T: num_traits::Num
        + num_traits::ConstOne
        + num_traits::ConstZero
        + std::fmt::Debug
        + num_traits::FromPrimitive
        + Ord
        + Copy
        + std::borrow::Borrow<I>
        + std::borrow::BorrowMut<I>
{
}

trait BackingKind:
    num_traits::Num
    + num_traits::ConstOne
    + num_traits::ConstZero
    + std::fmt::Debug
    + num_traits::FromPrimitive
    + num_traits::ToPrimitive
    + Ord
    + Copy
{
}

impl<T> BackingKind for T where
    T: num_traits::Num
        + num_traits::ConstOne
        + num_traits::ConstZero
        + std::fmt::Debug
        + num_traits::FromPrimitive
        + num_traits::ToPrimitive
        + Ord
        + Copy
{
}

fn simple_constants<T, I, const N: usize>()
where
    T: TestKind<I>,
    I: BackingKind,
{
    let raw = |t: &T| (t.borrow().to_i64().unwrap() as u64);
    let set_raw = |t: &mut T, value: u64| {
        *(t.borrow_mut()) = I::from_u64(value).unwrap();
    };
    let int_bits = std::mem::size_of::<I>() * 8 - N - 1;

    let z = T::ZERO;
    let o = T::ONE;

    eprintln!(
        "Simple_constants test on type {} with {N} fractional bits {int_bits} integer bits",
        std::any::type_name::<I>()
    );

    assert!(z.is_zero(), "Zero must be zero");
    assert!(o.is_one(), "One must be one");

    eprintln!("Zero {z:?} One {o:?}");

    assert_ne!(z, o, "One is not zero");

    assert_eq!(z, T::from_u8(0).unwrap(), "Zero must be zero");
    assert_eq!(o, T::from_u8(1).unwrap(), "One must be one");

    assert_eq!(raw(&z), 0, "Zero encoded as 0");
    assert_eq!(raw(&o), 1 << N, "One encoded as 1<<N");

    if int_bits > 1 {
        let t = o + o;
        assert_eq!(t, T::from_u8(2).unwrap(), "Two must be two");

        eprintln!("Zero {z:?} One {o:?} Two {t:?}");

        assert!(t > o, "2>1");
        assert!(o < t, "1<2");

        if int_bits > 2 {
            let f = t + t;
            assert_eq!(f, T::from_u8(4).unwrap(), "Four must be four");

            assert!(f > t, "4>2");
            assert!(t < f, "2<4");

            let s = t + f;
            assert_eq!(s, T::from_u8(6).unwrap(), "Six must be six");

            assert!(s > f, "6>4");
            assert!(f < s, "4<6");
        }
    }

    let mut half = T::ZERO;
    set_raw(&mut half, 1 << (N - 1));
    eprintln!("Half {half:?}");
    assert_eq!(half + half, o, "Two halves added is one");
    assert_eq!(raw(&half), 1 << (N - 1), "Half encoded as 1<<(N-1)");

    let mut smallest = T::ZERO;
    set_raw(&mut smallest, 1);
    assert_eq!(raw(&smallest), 1, "Smallest value must be 1");
}

#[test]
fn basic() {
    simple_constants::<Fixed<i8, 4>, i8, 4>();
    simple_constants::<Fixed<i8, 2>, i8, 2>();
    simple_constants::<Fixed<i8, 6>, i8, 6>();

    simple_constants::<Fixed<i32, 16>, i32, 16>();

    assert!(false, "Forced failure");
}
