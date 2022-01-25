/// The type we use for indexing into our arrays.
type Idx = usize;
/// The type our arrays hold. 32-bit are enough for genomic data.
type Val = u32;

/// Range Minimum Query interface.
pub trait RMQArray: std::ops::Index<Idx, Output = Val> {
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

impl BruteForceRMQ {
    pub fn new(values: Vec<u32>) -> BruteForceRMQ {
        BruteForceRMQ { values }
    }
}

impl RMQArray for BruteForceRMQ {
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

/// Implements RMQ by table lookup. Has a complete table of all [i,j),
/// so uses O(n²) memory, takes O(n²) time preprocessing, but then
/// does RMQ in O(1).
pub struct FullTabulateRMQ {
    values: Vec<u32>,
}
adapt_index!(FullTabulateRMQ); // Hack

impl FullTabulateRMQ {
    pub fn new(values: Vec<u32>) -> FullTabulateRMQ {
        FullTabulateRMQ { values }
    }
}

impl RMQArray for FullTabulateRMQ {
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
        let v = rmqa[k];
        for l in i..k {
            assert!(rmqa[l] > v);
        }
        for l in k + 1..j {
            assert!(rmqa[l] >= v)
        }
    }

    fn check_min<R: RMQArray>(rmqa: &R) {
        for i in 0..rmqa.len() - 1 {
            for j in i + 1..rmqa.len() {
                check_min_in_interval(rmqa, i, j)
            }
        }
    }

    #[test]
    fn test_index() {
        // FIXME: figure out how to get a random vector here
        let v = vec![2, 1, 2, 5, 3, 6, 1, 3, 7, 4];
        let rmqa = BruteForceRMQ::new(v.clone());
        check_same_index(&rmqa, &v);
    }

    #[test]
    fn test_rmq() {
        // FIXME: figure out how to get a random vector here
        let v = vec![2, 1, 2, 5, 3, 6, 1, 3, 7, 4];
        let rmqa = BruteForceRMQ::new(v.clone());
        check_min(&rmqa);
    }
}
