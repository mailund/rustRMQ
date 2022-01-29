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
use std::cmp::max;

/// Tests if x is a power of two, x=2^k.
pub fn power_of_two(x: Idx) -> bool {
    (x == 0) || ((x & (x - 1)) == 0)
}

/// Get k such that 2**k is j rounded down to the
/// nearest power of 2.
/// j=1=2^0 => 0
/// j=2=2^1 => 1
/// j=3=2^1+1 => 1
/// j=4=2^2 => 2
/// and so on.
pub fn log2_down(j: Idx) -> Idx {
    assert!(j != 0); // not defined for zero

    // Rounded down means finding the index of the first
    // 1 in the bit-pattern. If j = 00010101110
    // then 00010000000 (only first bit) is the closest
    // power of two, and we want the position of that bit.
    // j.leading_zeros() counts the number of leading zeros
    // and we get the index by subtracting this
    // from the total number of bits minus one.
    (Idx::BITS - j.leading_zeros() - 1) as Idx
    // Idx::BITS and j.leading_zeros() will be u32, so
    // we cast the result back to Idx.
}

/// For n, get (rounded up) log2(n).
/// We need this function for computing the size of tables with
/// log(n) entries. Although n=1 = 2^0, we can't use zero
/// for a table that must contain at least one element, so we
/// always return at least 1.
pub fn log2_up(n: Idx) -> Idx {
    // Round down, but add one if n is not a power of two.
    // We have to always return at least 1, to handle arrays of length
    // one, even though 1 = 2^0 is a power of two.
    let k = log2_down(n);
    let add = !power_of_two(n) as Idx; // 1 if n is not a power of two.
    max(1, k + add)
}

pub fn log_table_size(n: Idx) -> Idx {
    log2_down(n) + 1
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
        assert_eq!(0, log2_down(1));
        assert_eq!(1, log2_down(2));
        assert_eq!(1, log2_down(3));
        assert_eq!(2, log2_down(4));
        assert_eq!(2, log2_down(5));
        assert_eq!(2, log2_down(6));
        assert_eq!(2, log2_down(7));
        assert_eq!(3, log2_down(8));
        assert_eq!(3, log2_down(9));
        for i in 1..100 {
            let k = log2_down(i);
            assert_le!(1 << k, i);
            if power_of_two(i) {
                assert_eq!(i, 1 << k);
            }
        }
    }

    #[test]
    fn test_log2_up() {
        assert_eq!(1, log2_up(1));
        assert_eq!(1, log2_up(2));
        assert_eq!(2, log2_up(3));
        assert_eq!(2, log2_up(4));
        assert_eq!(3, log2_up(5));
        assert_eq!(3, log2_up(6));
        assert_eq!(3, log2_up(7));
        assert_eq!(3, log2_up(8));
        assert_eq!(4, log2_up(9));
        for i in 1..100 {
            let k = log2_up(i);
            assert_le!(i, 1 << k);
            if i > 1 && power_of_two(i) {
                assert_eq!(i, 1 << k);
            }
        }
        for k in 2..10 {
            let i = 1 << k;
            assert_eq!(log2_up(i), k);
            assert_eq!(log2_down(i), k);
        }
    }

    #[test]
    fn test_log_table_size() {
        assert_eq!(1, log_table_size(1));
        assert_eq!(2, log_table_size(2));
        assert_eq!(2, log_table_size(3));
        assert_eq!(3, log_table_size(4));
        assert_eq!(3, log_table_size(5));
        assert_eq!(3, log_table_size(6));
        assert_eq!(3, log_table_size(7));
        assert_eq!(4, log_table_size(8));
        assert_eq!(4, log_table_size(9));
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
