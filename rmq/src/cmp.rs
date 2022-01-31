/// Takes (a,b,op(a,b)) and lift it
macro_rules! lift_binop {
    ($a:ident, $b:ident, $expr:expr) => {
        match (&$a, &$b) {
            (&None, &None) => None,
            (&Some(_), &None) => $a,
            (&None, &Some(_)) => $b,
            (&Some($a), &Some($b)) => Some($expr),
        }
    };
}

pub trait Min {
    fn min(a: Self, b: Self) -> Self;
}

#[inline]
pub fn min<T: Min>(a: T, b: T) -> T {
    T::min(a, b)
}

#[inline]
pub fn min3<T: Min>(a: T, b: T, c: T) -> T {
    min(min(a, b), c)
}

// Lift min to Option<T>
impl<T> Min for Option<T>
where
    T: Min + Copy,
{
    #[inline]
    fn min(a: Self, b: Self) -> Self {
        lift_binop!(a, b, T::min(a, b))
    }
}
