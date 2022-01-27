mod math;
use math::log2_up;

mod rmq_array;
#[allow(unused_imports)] // Importing RMQArray to expose it; I don't use it here
use rmq_array::RMQArray;
use rmq_array::{Idx, RMQArrayImpl, RMQArray_, Val};

mod inter_table;
use inter_table::InterTable;

mod power_table;
use power_table::{adjusted_index, TwoD};

/// Finds the left-most index with the smallest value in x.
/// The index returned is indexed from the start of x, so if x
/// is a sub-range of some y, x = y[i..j], you should add i
/// to get the index in y.
fn smallest_in_range(x: &[Val]) -> Idx {
    let (k, _) = x
        .iter()
        .enumerate()
        .fold((0, &std::u32::MAX as &u32), |(k1, v1), (k2, v2)| {
            if v1 <= v2 {
                (k1, v1)
            } else {
                (k2, v2)
            }
        });
    k
}

/// Implements RMQ by running through the [i,j) interval
/// and finding the index with the smallest value.
/// O(1) preprocessing and O(n) access.
pub struct BruteForceRMQImpl {
    values: Vec<u32>,
}

impl RMQArrayImpl for BruteForceRMQImpl {
    fn new(values: Vec<Val>) -> BruteForceRMQImpl {
        BruteForceRMQImpl { values }
    }

    fn len(&self) -> Idx {
        self.values.len()
    }
    fn val(&self, index: Idx) -> &Val {
        &self.values[index]
    }
    fn rmq(&self, i: Idx, j: Idx) -> Idx {
        // smallest_in_range() is relative to interval, so offset with start
        i + smallest_in_range(&self.values[i..j])
    }
}

pub type BruteForceRMQ = RMQArray_<BruteForceRMQImpl>;

/// Implements RMQ by table lookup. Has a complete table of all [i,j),
/// so uses O(n²) memory, takes O(n²) time preprocessing, but then
/// does RMQ in O(1).
pub struct FullTabulateRMQImpl {
    values: Vec<u32>,
    rmq: InterTable,
}
pub type FullTabulateRMQ = RMQArray_<FullTabulateRMQImpl>;

impl RMQArrayImpl for FullTabulateRMQImpl {
    fn new(values: Vec<u32>) -> FullTabulateRMQImpl {
        let mut rmq = InterTable::new(values.len());
        for i in 0..values.len() {
            rmq[(i, i + 1)] = i;
        }
        for i in 0..values.len() - 1 {
            for j in i + 2..values.len() + 1 {
                let k = rmq[(i, j - 1)];
                let v1 = values[k];
                let v2 = values[j - 1];
                rmq[(i, j)] = if v1 <= v2 { k } else { j - 1 }
            }
        }
        FullTabulateRMQImpl { values, rmq }
    }
    fn len(&self) -> Idx {
        self.values.len()
    }
    fn val(&self, index: Idx) -> &Val {
        &self.values[index]
    }
    fn rmq(&self, i: Idx, j: Idx) -> Idx {
        self.rmq[(i, j)]
    }
}

/// RMQ table that tabulates all [i,i+2^k] ranges (there are O(n log n)),
/// form which we can get the RMQ from the table by splitting [i,j) into
/// two, [i,2^k) and [j-2^k,j) (where k is the largest such k). We can get
/// the RMQ from those two intervals with a table lookup in O(1) and then
/// pick the one of those with the smallest value.
/// The result is O(n log n) preprocessing and O(1) lookup.
pub struct PowerRMQImpl {
    lcp: Vec<Val>,
    tbl: TwoD,
}

impl RMQArrayImpl for PowerRMQImpl {
    fn new(lcp: Vec<u32>) -> PowerRMQImpl {
        let n = lcp.len();
        let logn = log2_up(n);
        let mut tbl = TwoD::new(n);
        for i in 0..n {
            tbl[(i, 0)] = i;
        }
        // Dynamic programming construction of tables of increasing length.
        // We have O(log n) runs of the outer loop and O(n) of the inner,
        // so the total time is O(n log n).
        for k in 1..logn {
            for i in 0..(n - k) {
                let i1 = tbl[(i, k - 1)];
                let i2 = tbl[(i + k, k - 1)];
                let v1 = lcp[i1];
                let v2 = lcp[i2];
                tbl[(i, k)] = if v1 <= v2 { i1 } else { i2 }
            }
        }
        PowerRMQImpl { lcp, tbl }
    }
    fn len(&self) -> Idx {
        self.lcp.len()
    }
    fn val(&self, index: Idx) -> &Val {
        &self.lcp[index]
    }
    fn rmq(&self, i: Idx, j: Idx) -> Idx {
        // Work out k so [i,2^k) and [j-2^k,j) are overlapping (and are not overlapping)
        // anything outside of [i,j). Then use the table to get the index with the smallest
        // lcp in those intervals, and pick the smaller of the two (with the first index
        // in case of a tie). All in O(1).
        let (k, ii) = adjusted_index(i, j);
        let i1 = self.tbl[(i, k)];
        let i2 = self.tbl[(ii, k)];
        let v1 = self.lcp[i1];
        let v2 = self.lcp[i2];
        if v1 <= v2 {
            i1
        } else {
            i2
        }
    }
}

pub type PowerRMQ = RMQArray_<PowerRMQImpl>;

#[cfg(test)]
mod tests {
    use super::*;

    fn check_same_index<R: RMQArray>(rmqa: &R, vals: &Vec<u32>) {
        assert_eq!(rmqa.len(), vals.len());
        for (i, &v) in vals.iter().enumerate() {
            assert_eq!(v, rmqa[i]);
        }
    }

    fn check_min_in_interval<R: RMQArray>(rmqa: &R, i: Idx, j: Idx) {
        let k = rmqa.rmq(i, j);
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
        let v = vec![2, 1, 2, 5, 3, 6, 1, 3, 7, 4];
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

    fn check_rmq<R: RMQArray>() {
        // FIXME: figure out how to get a random vector here
        let v = vec![2, 1, 2, 5, 3, 6, 1, 3, 7, 4];
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
        // First a few checks of the Power specific table...
        // can we handle the diagonal (base case of the dynamic programming),
        // and can we handle the cases where we only look up in the table?
        let v = vec![2, 1, 2, 5, 3, 6, 1, 3, 7, 4];
        let rmqa = PowerRMQ::new(v.clone());
        println!("{:?}", &v);
        println!("{}", rmqa.tbl);

        // Checking diagonal
        for i in 0..v.len() {
            assert_eq!(i, rmqa.rmq(i, i + 1));
        }

        // Checking powers
        for i in 0..v.len() {
            for k in [0, 1, 2, 3] {
                let j = i + (1 << k);
                if j > v.len() {
                    continue;
                }
                let i1 = i + smallest_in_range(&v[i..j]);
                let i2 = rmqa.rmq(i, j);
                println!(
                    "[{},{}): {}, {}, [offset={}] {:?}",
                    i,
                    j,
                    i1,
                    i2,
                    i,
                    &v[i..j]
                );
                assert_eq!(i1, i2);
            }
        }

        // Full check
        check_rmq::<PowerRMQ>()
    }
}
