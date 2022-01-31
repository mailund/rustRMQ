use super::cmp::min;
use super::inter_table::InterTable;
use super::interval::{Idx, Point, Val};
use super::rmq_array::RMQArrayImpl;

/// Implements RMQ by table lookup. Has a complete table of all [i,j),
/// so uses O(n²) memory, takes O(n²) time preprocessing, but then
/// does RMQ in O(1).
pub struct FullTabulateRMQImpl {
    values: Vec<u32>,
    rmq: InterTable,
}

impl RMQArrayImpl for FullTabulateRMQImpl {
    fn new(values: Vec<u32>) -> FullTabulateRMQImpl {
        let mut rmq = InterTable::new(values.len());
        for i in 0..values.len() {
            rmq[(i, i + 1)] = i;
        }
        for i in 0..values.len() - 1 {
            for j in i + 2..values.len() + 1 {
                // Dynamic programming:
                // Min val in [i,j) is either min in [i,j-1) or
                // [j-1,j) (i.e. index j-1).
                let left = Point::new(rmq[(i, j - 1)], &values);
                let current = Point::new(j - 1, &values);
                rmq[(i, j)] = min(left, current).idx()
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
    fn rmq(&self, i: Idx, j: Idx) -> Point {
        let k = self.rmq[(i, j)];
        Point(k, self.values[k])
    }
}
