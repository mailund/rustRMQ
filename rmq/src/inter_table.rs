use super::*;

#[inline]
fn flat_idx(i: Idx, j: Idx, n: Idx) -> Idx {
    let k = n - i - 1;
    k * (k + 1) / 2 + j - i - 1
}

/// Table for looking up at [i,j) (j > i) intervals.
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
