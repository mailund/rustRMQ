mod brute_force;
mod cmp;
mod full_table;
mod inter_table;
mod interval;
mod math;
mod power_table;
mod powers;
mod reduce;
mod rmq_array;

#[allow(unused_imports)] // Importing RMQArray to expose it; I don't use it here
use rmq_array::RMQArray;
use rmq_array::RMQArray_;

/// Implements RMQ by running through the [i,j) interval
/// and finding the index with the smallest value.
/// O(1) preprocessing and O(n) access.
pub type BruteForceRMQ = RMQArray_<brute_force::BruteForceRMQImpl>;

/// Implements RMQ by table lookup. Has a complete table of all [i,j),
/// so uses O(n²) memory, takes O(n²) time preprocessing, but then
/// does RMQ in O(1).
pub type FullTabulateRMQ = RMQArray_<full_table::FullTabulateRMQImpl>;

/// RMQ table that tabulates all [i,i+2^k] ranges (there are O(n log n)),
/// form which we can get the RMQ from the table by splitting [i,j) into
/// two, [i,2^k) and [j-2^k,j) (where k is the largest such k). We can get
/// the RMQ from those two intervals with a table lookup in O(1) and then
/// pick the one of those with the smallest value.
/// The result is O(n log n) preprocessing and O(1) lookup.
pub type PowerRMQ = RMQArray_<powers::PowerRMQImpl>;

/// Reduces the input vector to one of length m = n/log(n); then preprocess
/// it with the PowerRMQ method (in time m log m = n/log(n) log(n/log(n)) = O(n)).
/// Lookup is now one constant time lookup in the PowerRMQ and two linear searches
/// in blocks of length log n. So, preprocessing O(n) and lookup O(log n).
pub type ReducedPowerRMQ = RMQArray_<reduce::ReducedPowerRMQImpl>;

#[cfg(test)]
mod tests {
    use super::interval::{Idx, Point};
    use super::*;

    fn check_same_index<R: RMQArray>(rmqa: &R, vals: &Vec<u32>) {
        assert_eq!(rmqa.len(), vals.len());
        for (i, &v) in vals.iter().enumerate() {
            assert_eq!(v, rmqa[i]);
        }
    }

    fn check_min_in_interval<R: RMQArray>(rmqa: &R, i: Idx, j: Idx) {
        let Point(k, _) = rmqa.rmq(i, j);
        assert!(i <= k);
        assert!(k < j);

        let v = rmqa[k];
        for l in i..k {
            assert!(rmqa[l] > v);
        }
        for l in k + 1..j {
            assert!(rmqa[l] >= v)
        }
    }

    fn check_min<R: RMQArray>(rmqa: &R) {
        for i in 0..rmqa.len() {
            for j in i + 1..rmqa.len() + 1 {
                check_min_in_interval(rmqa, i, j)
            }
        }
    }

    fn check_index<R: RMQArray>() {
        // Not power of two
        let v = vec![2, 1, 2, 5, 3, 6, 1, 3, 7, 4];
        let rmqa = R::new(v.clone());
        check_same_index(&rmqa, &v);
        // Power of two
        let v = vec![2, 1, 2, 5, 3, 6, 1, 3, 7, 4, 2, 6, 3, 4, 7, 9];
        let rmqa = R::new(v.clone());
        check_same_index(&rmqa, &v);
    }

    #[test]
    fn test_index_brute() {
        check_index::<BruteForceRMQ>()
    }

    #[test]
    fn test_index_full() {
        check_index::<FullTabulateRMQ>()
    }

    #[test]
    fn test_index_power() {
        check_index::<PowerRMQ>()
    }

    #[test]
    fn test_index_reduced_power() {
        check_index::<ReducedPowerRMQ>()
    }

    fn check_rmq<R: RMQArray>() {
        // Not power of two
        let v = vec![2, 1, 2, 5, 3, 6, 1, 3, 7, 4];
        let rmqa = R::new(v.clone());
        check_min(&rmqa);
        // Power of two
        let v = vec![2, 1, 2, 5, 3, 6, 1, 3, 7, 4, 2, 6, 3, 4, 7, 9];
        let rmqa = R::new(v.clone());
        check_min(&rmqa);
        // Not power of two
        let v = vec![2, 1, 2, 0, 2, 1, 3, 7, 4];
        let rmqa = R::new(v.clone());
        check_min(&rmqa);
        // Power of two
        let v = vec![2, 1, 2, 5, 3, 6, 1, 3];
        let rmqa = R::new(v.clone());
        check_min(&rmqa);
    }

    #[test]
    fn test_rmq_brute() {
        check_rmq::<BruteForceRMQ>()
    }

    #[test]
    fn test_rmq_full() {
        check_rmq::<FullTabulateRMQ>()
    }

    #[test]
    fn test_rmq_power() {
        check_rmq::<PowerRMQ>()
    }

    #[test]
    fn test_rmq_reduced_power() {
        check_rmq::<ReducedPowerRMQ>()
    }
}
