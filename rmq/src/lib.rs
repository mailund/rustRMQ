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
        let (k, _) = self.values[i..j].iter().enumerate().fold(
            (0, &std::u32::MAX as &u32),
            |(k1, v1), (k2, v2)| if v1 <= v2 { (k1, v1) } else { (k2, v2) },
        );
        i + k // k is relative to interval, so offset with start
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
}
