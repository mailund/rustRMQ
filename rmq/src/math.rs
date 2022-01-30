//! The log functions in this module requires that the index is at least 1
//! since log(0) is not well-defined. This works out well, though, since we
//! use them for j in intervals [i,j) where j > i, so j >= 1. For allocating
//! memory for tables that address in log-space, however, it means that
//! if you want to index k (j == 2^k+jj) you need at least k entries, so
//! you cannot use log2_up when j is a power of two (log2_up(2^k) == k)
//! as this will be one too little. To allocate memory, use log_table_size(j)
//! which will be log2_down(j) + 1. When j is not a power of two, this is
//! the same as log2_up(j), but for powers of two, it adds the extra entry.

use super::Idx;

/// Tests if x is a power of two, x=2^k.
pub fn power_of_two(x: Idx) -> bool {
    (x == 0) || ((x & (x - 1)) == 0)
}

/// Type for powers of two, 2^k. Contains k, but wrapped in
/// a type so we don't confuse log-space with linear space.
#[derive(Debug, Clone, Copy)]
pub struct Pow(pub Idx);

impl std::cmp::PartialEq for Pow {
    #[inline]
    fn eq(&self, other: &Pow) -> bool {
        self.0 == other.0
    }
}

impl std::ops::Add for Pow {
    type Output = Pow;
    fn add(self, rhs: Self) -> Self::Output {
        let (Pow(k), Pow(kk)) = (self, rhs);
        Pow(k + kk)
    }
}

impl std::fmt::Display for Pow {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "2^{}", self.0)
    }
}

impl Pow {
    /// for a power Pow(k) get 2^k.
    #[inline]
    pub fn value(&self) -> Idx {
        1 << self.0
    }
    /// for a power Pow(k) get k.
    #[inline]
    pub fn exponent(&self) -> Idx {
        self.0
    }
}

/// Get k such that 2**k is j rounded down to the
/// nearest power of 2.
/// j=1=2^0 => 0
/// j=2=2^1 => 1
/// j=3=2^1+1 => 1
/// j=4=2^2 => 2
/// and so on.
pub fn log2_down(j: Idx) -> Pow {
    assert!(j != 0); // not defined for zero

    // Rounded down means finding the index of the first
    // 1 in the bit-pattern. If j = 00010101110
    // then 00010000000 (only first bit) is the closest
    // power of two, and we want the position of that bit.
    // j.leading_zeros() counts the number of leading zeros
    // and we get the index by subtracting this
    // from the total number of bits minus one.
    Pow((Idx::BITS - j.leading_zeros() - 1) as Idx)
    // Idx::BITS and j.leading_zeros() will be u32, so
    // we cast the result back to Idx.
}

/// We always have to add one to the exponent, because in log-space
/// we are working with 1-indexed (0-indexed in log-space) values,
/// so to have a table that can handle maximum value k, we need k+1
/// entires. That is what this function gives us.
pub fn log_table_size(n: Idx) -> Pow {
    let Pow(k) = log2_down(n);
    Pow(k + 1)
}

/// For n, get (rounded up) log2(n).
pub fn log2_up(n: Idx) -> Pow {
    // log_table_size(n) with n=2^k+m will always give us 2^{k+1},
    // whether m is zero or not. We want 2^{k+1} when m > 0 and 2^k
    // when m is zero, i.e. when n is a power of two.
    // So we should subtract one from the exponent if n is a power of two.
    let Pow(k) = log_table_size(n);
    Pow(k - power_of_two(n) as Idx)
}

/// For n and block size bs, compute (r,r*n) where r
/// is n/bs rounded down. That is, r is n divided by bs
/// rounded down, and r*bs is n adjusted downwards to the
/// closest multiple of bs.
pub fn round_down(n: Idx, bs: usize) -> (Idx, Idx) {
    let r = n / bs;
    (r, r * bs)
}

/// For n and block size bs, compute (r,r*n) where r
/// is n/bs rounded up. That is, r is n divided by bs
/// rounded down, and r*bs is n adjusted upwards to the
/// closest multiple of bs.
pub fn round_up(n: Idx, bs: usize) -> (Idx, Idx) {
    let r = (n + bs - 1) / bs;
    (r, r * bs)
}

#[cfg(test)]
mod tests {
    use super::*;
    use more_asserts::*;

    #[test]
    fn test_power_of_two() {
        assert!(power_of_two(0));
        assert!(power_of_two(1));
        assert!(power_of_two(2));
        assert!(!power_of_two(3));
        assert!(power_of_two(4));
        assert!(!power_of_two(5));
        assert!(!power_of_two(6));
        assert!(!power_of_two(7));
        assert!(power_of_two(8));
        assert!(!power_of_two(9));
    }

    #[test]
    fn test_log2_down() {
        assert_eq!(Pow(0), log2_down(1));
        assert_eq!(Pow(1), log2_down(2));
        assert_eq!(Pow(1), log2_down(3));
        assert_eq!(Pow(2), log2_down(4));
        assert_eq!(Pow(2), log2_down(5));
        assert_eq!(Pow(2), log2_down(6));
        assert_eq!(Pow(2), log2_down(7));
        assert_eq!(Pow(3), log2_down(8));
        assert_eq!(Pow(3), log2_down(9));
        for i in 1..100 {
            let k = log2_down(i);
            assert_le!(k.value(), i);
            if power_of_two(i) {
                assert_eq!(i, k.value());
            }
        }
    }

    #[test]
    fn test_log2_up() {
        assert_eq!(0, log2_up(1).exponent());
        assert_eq!(1, log2_up(2).exponent());
        assert_eq!(2, log2_up(3).exponent());
        assert_eq!(2, log2_up(4).exponent());
        assert_eq!(3, log2_up(5).exponent());
        assert_eq!(3, log2_up(6).exponent());
        assert_eq!(3, log2_up(7).exponent());
        assert_eq!(3, log2_up(8).exponent());
        assert_eq!(4, log2_up(9).exponent());
        for i in 1..100 {
            let k = log2_up(i);
            assert_le!(i, k.value());
            if i > 1 && power_of_two(i) {
                assert_eq!(i, k.value());
            }
        }
        for k in 2..10 {
            let i = 1 << k;
            assert_eq!(log2_up(i).exponent(), k);
            assert_eq!(log2_down(i).exponent(), k);
        }
    }

    #[test]
    fn test_log_table_size() {
        assert_eq!(Pow(1), log_table_size(1));
        assert_eq!(Pow(2), log_table_size(2));
        assert_eq!(Pow(2), log_table_size(3));
        assert_eq!(Pow(3), log_table_size(4));
        assert_eq!(Pow(3), log_table_size(5));
        assert_eq!(Pow(3), log_table_size(6));
        assert_eq!(Pow(3), log_table_size(7));
        assert_eq!(Pow(4), log_table_size(8));
        assert_eq!(Pow(4), log_table_size(9));
    }

    #[test]
    fn test_round() {
        let bs = 4;
        assert_eq!((0, 0), round_down(0, bs));
        assert_eq!((0, 0), round_up(0, bs));
        assert_eq!((0, 0), round_down(1, bs));
        assert_eq!((1, 4), round_up(1, bs));

        assert_eq!((0, 0), round_down(2, bs));
        assert_eq!((1, 4), round_up(2, bs));

        assert_eq!((0, 0), round_down(3, bs));
        assert_eq!((1, 4), round_up(3, bs));

        assert_eq!((1, 4), round_down(4, bs));
        assert_eq!((1, 4), round_up(4, bs));

        assert_eq!((1, 4), round_down(5, bs));
        assert_eq!((2, 8), round_up(5, bs));

        assert_eq!((1, 4), round_down(6, bs));
        assert_eq!((2, 8), round_up(6, bs));
    }
}
