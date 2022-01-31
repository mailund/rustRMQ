use super::cmp::min3;
use super::interval::{smallest_in_range, Idx, Point, Val};
use super::math::{log2_up, round_down, round_up, Pow};
use super::rmq_array::RMQArrayImpl;

/// Reduce an array x to the smallest value in each block (of size block_size)
/// and the index in the original array that this minimal value sits at.
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

/// Wrapper around a power-RMQ with a block index interface.
/// It adds a bit of type safety when we have reduced the array
/// into blocks, so we don't confuse indices into blocks with
/// indices into the original array.
struct BlockRMQ(super::PowerRMQ);
use super::math::BlockIdx;

impl BlockRMQ {
    fn new(lcp: Vec<u32>) -> BlockRMQ {
        BlockRMQ(super::PowerRMQ::new(lcp))
    }
    fn rmq(&self, bi: BlockIdx, bj: BlockIdx) -> Point {
        let BlockIdx(i) = bi;
        let BlockIdx(j) = bj;
        self.0.rmq(i, j)
    }
}

pub struct ReducedPowerRMQImpl {
    lcp: Vec<Val>,
    block_size: usize,
    reduced_index: Vec<Idx>,
    tbl: BlockRMQ,
}

impl ReducedPowerRMQImpl {
    /// Move a point from the reduced sequence to the original sequence
    fn adjust_point(&self, p: Point) -> Point {
        let Point(ri, v) = p;
        Point(self.reduced_index[ri], v)
    }
}

/// Reduces the input vector to one of length m = n/log(n); then preprocess
/// it with the PowerRMQ method (in time m log m = n/log(n) log(n/log(n)) = O(n)).
/// Lookup is now one constant time lookup in the PowerRMQ and two linear searches
/// in blocks of length log n. So, preprocessing O(n) and lookup O(log n).
impl RMQArrayImpl for ReducedPowerRMQImpl {
    fn new(lcp: Vec<Val>) -> ReducedPowerRMQImpl {
        let Pow(block_size) = log2_up(lcp.len());
        let (reduced_index, reduced_values) = reduce_array(&lcp, block_size);
        let tbl = BlockRMQ::new(reduced_values);
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

    fn rmq(&self, i: Idx, j: Idx) -> Point {
        let (bi, ii) = round_up(i, self.block_size);
        let (bj, jj) = round_down(j, self.block_size);
        if bi < bj {
            let p1 = smallest_in_range(&self.lcp, i, ii);
            let p2 = Some(self.adjust_point(self.tbl.rmq(bi, bj)));
            let p3 = smallest_in_range(&self.lcp, jj, j);
            min3(p1, p2, p3)
        } else {
            smallest_in_range(&self.lcp, i, j)
        }
        .unwrap()
    }
}
