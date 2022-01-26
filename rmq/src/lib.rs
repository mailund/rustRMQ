/// The type we use for indexing into our arrays.
type Idx = usize;
/// The type our arrays hold. 32-bit are enough for genomic data.
type Val = u32;

/// Range Minimum Query interface.
pub trait RMQArray: std::ops::Index<Idx, Output = Val> {
    /// Create the necessary tables from a vector of values.
    fn new(values: Vec<Val>) -> Self;
    /// Get the length of the underlying array.
    fn len(&self) -> Idx;
    /// Get the value at a specific index.
    fn val(&self, index: Idx) -> &Val;
    /// Get the index of the first minimal value in the range [i,j).
    fn rmq(&self, i: Idx, j: Idx) -> Idx;
}

// Ugly hack for now! FIXME: Find a better solution to get a generic
// Index<> implementation.
macro_rules! adapt_index {
    ($t:ident) => {
        impl std::ops::Index<Idx> for $t {
            type Output = Val;
            fn index(&self, index: Idx) -> &Val {
                self.val(index)
            }
        }
    };
}

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
pub struct BruteForceRMQ {
    values: Vec<u32>,
}
adapt_index!(BruteForceRMQ); // Hack

impl RMQArray for BruteForceRMQ {
    fn new(values: Vec<Val>) -> BruteForceRMQ {
        BruteForceRMQ { values }
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

mod interval_table {
    use super::*;
    fn flat_idx(i: Idx, j: Idx, n: Idx) -> Idx {
        let k = n - i - 1;
        k * (k + 1) / 2 + j - i - 1
    }

    pub struct InterTable {
        n: Idx,
        table: Vec<Idx>,
    }

    impl InterTable {
        pub fn new(n: Idx) -> InterTable {
            let table: Vec<Idx> = vec![0; n * (n + 1) / 2];
            InterTable { n, table }
        }
    }

    impl std::ops::Index<(Idx, Idx)> for InterTable {
        type Output = Idx;
        fn index(&self, index: (Idx, Idx)) -> &Self::Output {
            let (i, j) = index;
            assert!(i < self.n);
            assert!(i < j && j <= self.n);
            &self.table[flat_idx(i, j, self.n)]
        }
    }

    impl std::ops::IndexMut<(Idx, Idx)> for InterTable {
        fn index_mut(&mut self, index: (Idx, Idx)) -> &mut Self::Output {
            let (i, j) = index;
            assert!(i < self.n);
            assert!(i < j && j <= self.n);
            &mut self.table[flat_idx(i, j, self.n)]
        }
    }
}

use interval_table::InterTable;

/// Implements RMQ by table lookup. Has a complete table of all [i,j),
/// so uses O(n²) memory, takes O(n²) time preprocessing, but then
/// does RMQ in O(1).
pub struct FullTabulateRMQ {
    values: Vec<u32>,
    rmq: InterTable,
}
adapt_index!(FullTabulateRMQ); // Hack

impl RMQArray for FullTabulateRMQ {
    fn new(values: Vec<u32>) -> FullTabulateRMQ {
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
        FullTabulateRMQ { values, rmq }
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

mod power_table {
    use super::*;
    /// Get k such that 2**k is j rounded down to the
    /// nearest power of 2.
    /// Can't use Idx here because I need Idx==u32, but
    /// I'll just use it as u32 and expect the type checker
    /// to complain if Idx ever changes type.
    /// We don't handle zero, because it is a special case
    /// and we don't need to round down zero.
    fn round_down_offset(j: Idx) -> Idx {
        assert!(j != 0);
        (32 - (j as u32).leading_zeros() - 1) as Idx
    }

    /// For n, get (rounded up) log(n)
    pub fn log(n: Idx) -> Idx {
        round_down_offset(n) + 1 // might go one too high; fix later.
    }

    /// From range [i,j), get values (k,j') where k is the offset
    /// into the TwoD table to look up the value for [i,j') value there,
    /// and where j' is the rounded down j, which is the index we should
    /// linear search for in lcp[j'..j].
    pub fn adjusted_index(i: Idx, j: Idx) -> (Idx, Idx) {
        let k = round_down_offset(j - i);
        (k, i + (1 << k))
    }

    /// A rather simple 2D array made from vectors of vectors.
    /// There are better solutions, but I can implement those later
    /// with the same interface.
    pub struct TwoD {
        table: Vec<Vec<Idx>>,
    }

    impl TwoD {
        pub fn new(n: Idx) -> TwoD {
            let table = vec![vec![0; log(n)]; n];
            TwoD { table }
        }
    }

    impl std::ops::Index<(Idx, Idx)> for TwoD {
        type Output = Idx;
        fn index(&self, index: (Idx, Idx)) -> &Self::Output {
            let (i, j) = index;
            &self.table[i][j]
        }
    }

    impl std::ops::IndexMut<(Idx, Idx)> for TwoD {
        fn index_mut(&mut self, index: (Idx, Idx)) -> &mut Self::Output {
            let (i, j) = index;
            &mut self.table[i][j]
        }
    }

    impl std::fmt::Display for TwoD {
        fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
            for row in &self.table {
                for val in row {
                    let _ = write!(f, "{} ", val);
                }
                let _ = write!(f, "\n");
            }
            Ok(())
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;

        #[test]
        fn test_round_down_offset() {
            assert_eq!(0, round_down_offset(1));
            assert_eq!(1, round_down_offset(2));
            assert_eq!(1, round_down_offset(3));
            assert_eq!(2, round_down_offset(4));
            assert_eq!(2, round_down_offset(5));
            assert_eq!(2, round_down_offset(6));
            assert_eq!(2, round_down_offset(7));
            assert_eq!(3, round_down_offset(8));
            assert_eq!(3, round_down_offset(9));
        }

        #[test]
        fn test_adjusted_index() {
            let (k, jj) = adjusted_index(0, 1);
            assert_eq!(k, 0);
            assert_eq!(jj, 1);

            let (k, jj) = adjusted_index(0, 2);
            assert_eq!(k, 1);
            assert_eq!(jj, 2);

            let (k, jj) = adjusted_index(0, 3);
            assert_eq!(k, 1);
            assert_eq!(jj, 2);

            let (k, jj) = adjusted_index(0, 4);
            assert_eq!(k, 2);
            assert_eq!(jj, 4);

            let (k, jj) = adjusted_index(0, 5);
            assert_eq!(k, 2);
            assert_eq!(jj, 4);

            let (k, jj) = adjusted_index(0, 6);
            assert_eq!(k, 2);
            assert_eq!(jj, 4);

            let (k, jj) = adjusted_index(0, 7);
            assert_eq!(k, 2);
            assert_eq!(jj, 4);

            let (k, jj) = adjusted_index(0, 8);
            assert_eq!(k, 3);
            assert_eq!(jj, 8);
        }

        #[test]
        fn test_2d() {
            let n = 5;
            let mut tbl = TwoD::new(n);
            println!("{}", tbl);

            for i in 0..n {
                for j in i + 1..n + 1 {
                    let (k, _) = adjusted_index(i, j);
                    assert_eq!(0, tbl[(i, k)]);
                }
            }

            for i in 0..n {
                for j in i + 1..n + 1 {
                    // This is just saving the largest matching j
                    // in the entry those js should go to
                    let (k, _) = adjusted_index(i, j);
                    println!("({},{}) to offset {}", i, j, k);
                    tbl[(i, k)] = j;
                }
            }
            println!("{}", tbl);
            for i in 0..n {
                for j in i + 1..n + 1 {
                    let (k, _) = adjusted_index(i, j);
                    println!("({},{}): {} <? {}", i, j, j, tbl[(i, k)]);
                    assert!(j <= tbl[(i, k)]);
                }
            }
        }
    }
}

use power_table::{adjusted_index, log, TwoD};

/// RMQ table that tabulates all [i,i+2^k] ranges (there are O(n log n)),
/// form which we can get the RMQ from the table by splitting [i,j) into
/// two, [i,2^k) and [2^k,j) (where k is the largest such k). The first
/// we can get from the tables in O(1) and the second we can get with a
/// linear search in O(log n), since the range [2^k,j) can't be more than
/// O(log n) long if k is the largest such k.
/// The result is O(n log n) preprocessing and O(log n) lookup.
struct PowerRMQ {
    lcp: Vec<Val>,
    tbl: TwoD,
}
adapt_index!(PowerRMQ); // Hack

impl RMQArray for PowerRMQ {
    fn new(lcp: Vec<u32>) -> PowerRMQ {
        let n = lcp.len();
        let logn = log(n);
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
        PowerRMQ { lcp, tbl }
    }
    fn len(&self) -> Idx {
        self.lcp.len()
    }
    fn val(&self, index: Idx) -> &Val {
        &self.lcp[index]
    }
    fn rmq(&self, i: Idx, j: Idx) -> Idx {
        let (k, jj) = adjusted_index(i, j);
        let i1 = self.tbl[(i, k)];
        if jj == j {
            // If we do not have any left-over bits, we shouldn't
            // mess with our result. We already have the result
            // we need from the table.
            return i1;
        }
        // Otherwise, we need to check if there is a smaller value
        // in the flanking region. It can't be more than log(n)
        // long, so a linear search is O(log n).
        let i2 = jj + smallest_in_range(&self.lcp[jj..j]);
        let v1 = self.lcp[i1];
        let v2 = self.lcp[i2];
        if v1 <= v2 {
            i1
        } else {
            i2
        }
    }
}

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
