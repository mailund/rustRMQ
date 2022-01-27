use super::math::{log2_down, log2_up};
use super::Idx;

/// From range [i,j), get values (k,j-2^k) where k is the offset
/// into the TwoD table to look up the value for [i,i+2^k) and [j-2^k,j)
/// from which we can get the RMQ.
pub fn adjusted_index(i: Idx, j: Idx) -> (Idx, Idx) {
    let k = log2_down(j - i);
    (k, j - (1 << k))
}

/// A rather simple 2D array made from vectors of vectors.
/// There are better solutions, but I can implement those later
/// with the same interface.
pub struct TwoD {
    table: Vec<Vec<Idx>>,
}

impl TwoD {
    pub fn new(n: Idx) -> TwoD {
        let table = vec![vec![0; log2_up(n)]; n];
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
    fn test_log2_down() {
        assert_eq!(0, log2_down(1));
        assert_eq!(1, log2_down(2));
        assert_eq!(1, log2_down(3));
        assert_eq!(2, log2_down(4));
        assert_eq!(2, log2_down(5));
        assert_eq!(2, log2_down(6));
        assert_eq!(2, log2_down(7));
        assert_eq!(3, log2_down(8));
        assert_eq!(3, log2_down(9));
    }

    #[test]
    fn test_adjusted_index() {
        let (k, ii) = adjusted_index(0, 1);
        assert_eq!(k, 0);
        assert_eq!(ii, 0);

        let (k, ii) = adjusted_index(0, 2);
        assert_eq!(k, 1);
        assert_eq!(ii, 0);

        let (k, ii) = adjusted_index(0, 3);
        assert_eq!(k, 1);
        assert_eq!(ii, 1);

        let (k, ii) = adjusted_index(0, 4);
        assert_eq!(k, 2);
        assert_eq!(ii, 0);

        let (k, ii) = adjusted_index(0, 5);
        assert_eq!(k, 2);
        assert_eq!(ii, 1);

        let (k, ii) = adjusted_index(0, 6);
        assert_eq!(k, 2);
        assert_eq!(ii, 2);

        let (k, ii) = adjusted_index(0, 7);
        assert_eq!(k, 2);
        assert_eq!(ii, 3);

        let (k, ii) = adjusted_index(0, 8);
        assert_eq!(k, 3);
        assert_eq!(ii, 0);

        let (k, ii) = adjusted_index(1, 8);
        assert_eq!(k, 2);
        assert_eq!(ii, 4);

        let (k, ii) = adjusted_index(1, 9);
        assert_eq!(k, 3);
        assert_eq!(ii, 1);
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
