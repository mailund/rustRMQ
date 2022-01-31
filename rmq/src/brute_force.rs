use super::interval::{smallest_in_range, Idx, Point, Val};
use super::rmq_array::RMQArrayImpl;

/// Implements RMQ by running through the [i,j) interval
/// and finding the index with the smallest value.
/// O(1) preprocessing and O(n) access.
pub struct BruteForceRMQImpl {
    values: Vec<Val>,
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
    fn rmq(&self, i: Idx, j: Idx) -> Point {
        smallest_in_range(&self.values, i, j).unwrap()
    }
}
