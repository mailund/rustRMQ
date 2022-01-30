mod math;
use math::{log2_up, log_table_size, round_down, round_up};

mod rmq_array;
#[allow(unused_imports)] // Importing RMQArray to expose it; I don't use it here
use rmq_array::RMQArray;
use rmq_array::{Idx, RMQArrayImpl, RMQArray_, Val};

mod inter_table;
use inter_table::InterTable;

mod power_table;
use power_table::{adjusted_index, TwoD};

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
    T::min(T::min(a, b), c)
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

/// A point is an index with the corresponding value
#[derive(Clone, Copy)]
pub struct Point(Idx, Val);

impl Point {
    #[inline]
    pub fn idx(&self) -> Idx {
        self.0
    }
    #[inline]
    pub fn val(&self) -> Val {
        self.1
    }
}

impl Min for Point {
    #[inline]
    fn min(p1: Point, p2: Point) -> Point {
        match p1.val() <= p2.val() {
            true => p1,
            false => p2,
        }
    }
}

/// Finds the left-most index with the smallest value in x.
/// Returns the index of the left-most minimal value and the
/// minimal value. If [i,j) is not a valid interval, you ged
/// None.
fn smallest_in_range(x: &[Val], i: Idx, j: Idx) -> Option<Point> {
    let y = &x[i..j];
    let min_val = y.iter().min()?;
    let pos = i + y.iter().position(|a| a == min_val)?;
    Some(Point(pos, *min_val))
}

/// Picks the index with the smallest value
#[inline]
fn select_index(x: &[Val], i: Idx, j: Idx) -> Idx {
    match x[i] <= x[j] {
        true => i,
        false => j,
    }
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
        smallest_in_range(&self.values, i, j).unwrap().idx()
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
                rmq[(i, j)] = select_index(&values, k, j - 1)
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
        let logn = log_table_size(n);
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
                tbl[(i, k)] = select_index(&lcp, i1, i2)
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
        select_index(&self.lcp, i1, i2)
    }
}

pub type PowerRMQ = RMQArray_<PowerRMQImpl>;

mod reduce {
    use super::{smallest_in_range, Idx, Point, Val};
    pub fn reduce_array(x: &Vec<Val>, block_size: usize) -> (Vec<Idx>, Vec<Val>) {
        let mut indices: Vec<Idx> = Vec::new();
        let mut values: Vec<Val> = Vec::new();
        let no_blocks = x.len() / block_size;
        for block in 0..no_blocks {
            let block_start = block * block_size;
            let block_end = block_start + block_size;
            let Point(pos, val) = smallest_in_range(&x, block_start, block_end).unwrap();
            indices.push(pos);
            values.push(val);
        }
        (indices, values)
    }

    #[cfg(test)]
    mod tests {
        use super::*;

        #[test]
        fn test_reduce() {
            let bs = 3;
            let v = vec![3, 2, 6, 1, 7, 3, 10, 1, 6, 2, 1, 7, 0, 2];
            let (idx, val) = reduce_array(&v, bs);
            for (i, &pos) in idx.iter().enumerate() {
                assert_eq!(v[pos], val[i]);
            }
        }
    }
}

use reduce::reduce_array;

pub struct ReducedPowerRMQImpl {
    lcp: Vec<Val>,
    block_size: usize,
    reduced_index: Vec<Idx>,
    tbl: PowerRMQImpl,
}

/// Reduces the input vector to one of length m = n/log(n); then preprocess
/// it with the PowerRMQ method (in time m log m = n/log(n) log(n/log(n)) = O(n)).
/// Lookup is now one constant time lookup in the PowerRMQ and two linear searches
/// in blocks of length log n. So, preprocessing O(n) and lookup O(log n).
impl RMQArrayImpl for ReducedPowerRMQImpl {
    fn new(lcp: Vec<Val>) -> ReducedPowerRMQImpl {
        let n = lcp.len();
        let block_size = log2_up(n);
        let (reduced_index, reduced_values) = reduce_array(&lcp, block_size);
        let tbl = PowerRMQImpl::new(reduced_values);
        ReducedPowerRMQImpl {
            lcp,
            block_size,
            reduced_index,
            tbl,
        }
    }
    fn len(&self) -> Idx {
        self.lcp.len()
    }
    fn val(&self, index: Idx) -> &Val {
        &self.lcp[index]
    }

    fn rmq(&self, i: Idx, j: Idx) -> Idx {
        let (bi, ii) = round_up(i, self.block_size);
        let (bj, jj) = round_down(j, self.block_size);
        if bi < bj {
            let p1 = smallest_in_range(&self.lcp, i, ii);
            let p2 = {
                let ri = self.tbl.rmq(bi, bj);
                Some(Point(self.reduced_index[ri], *self.tbl.val(ri)))
            };
            let p3 = smallest_in_range(&self.lcp, jj, j);
            min3(p1, p2, p3)
        } else {
            smallest_in_range(&self.lcp, i, j)
        }
        .unwrap()
        .idx()
    }
}

pub type ReducedPowerRMQ = RMQArray_<ReducedPowerRMQImpl>;

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
        let v = vec![2, 1, 2, 5, 3, 6, 1, 3, 7, 4, 1, 2, 4, 5, 6, 7];
        let rmqa = PowerRMQ::new(v.clone());
        println!("{:?}", &v);
        println!("{}", rmqa.tbl);
        println!("{}", rmqa.rmq(0, rmqa.len()));

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
                let i1 = smallest_in_range(&v, i, j).unwrap().idx();
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

    #[test]
    fn test_rmq_reduced_power() {
        // Full check
        check_rmq::<ReducedPowerRMQ>()
    }
}
