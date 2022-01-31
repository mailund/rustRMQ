use super::cmp::min;
use super::interval::{Idx, Point, Val};
use super::math::{log_table_size, Pow};
use super::power_table::{adjusted_index, TwoD};
use super::rmq_array::RMQArrayImpl;

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

        // When tbl is a TwoD table, interpret tbl[i,Pow(k)] as containing
        // values (at powers of two) in the range [i,i+2^k).
        let Pow(logn) = log_table_size(n);
        let mut tbl = TwoD::new(n);

        // Base case: intervals [i,i+1) = [i,i+2^0).
        for i in 0..n {
            tbl[(i, Pow(0))] = i;
        }

        // Dynamic programming construction of tables of increasing length.
        // We have O(log n) runs of the outer loop and O(n) of the inner,
        // so the total time is O(n log n).
        for k in 1..logn {
            for i in 0..(n - Pow(k - 1).value()) {
                // Interval [i,i+2^k) = [i,i+2^{k-1}) [i+2^{k-1},(i+2^{k-1})+2^{k-1})
                let left = Point::new(tbl[(i, Pow(k - 1))], &lcp);
                let right = Point::new(tbl[(i + Pow(k - 1).value(), Pow(k - 1))], &lcp);
                tbl[(i, Pow(k))] = min(left, right).idx();
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
    fn rmq(&self, i: Idx, j: Idx) -> Point {
        // Work out k so [i,2^k) and [j-2^k,j) are overlapping (and are not overlapping)
        // anything outside of [i,j). Then use the table to get the index with the smallest
        // lcp in those intervals, and pick the smaller of the two (with the first index
        // in case of a tie). All in O(1).
        let (p, ii) = adjusted_index(i, j);
        min(self.point(self.tbl[(i, p)]), self.point(self.tbl[(ii, p)]))
    }
}

#[cfg(test)]
mod tests {
    use super::super::interval::*;
    use super::*;
    #[test]
    fn test_rmq_power() {
        // First a few checks of the Power specific table...
        // can we handle the diagonal (base case of the dynamic programming),
        // and can we handle the cases where we only look up in the table?
        let v = vec![2, 1, 2, 5, 3, 6, 1, 3, 7, 4, 1, 2, 4, 5, 6, 7];
        let rmqa = PowerRMQImpl::new(v.clone());
        println!("{:?}", &v);
        println!("{}", rmqa.tbl);
        println!("{}", rmqa.rmq(0, rmqa.len()));

        // Checking diagonal
        for i in 0..v.len() {
            assert_eq!(i, rmqa.rmq(i, i + 1).idx());
        }

        // Checking powers
        for i in 0..v.len() {
            for k in [0, 1, 2, 3] {
                let j = i + (1 << k);
                if j > v.len() {
                    continue;
                }
                let i1 = smallest_in_range(&v, i, j).unwrap().idx();
                let i2 = rmqa.rmq(i, j).idx();
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
    }
}
