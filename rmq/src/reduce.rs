use super::{smallest_in_range, Idx, Point, Val};

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
