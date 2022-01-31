use super::interval::Idx;
use super::math::{log2_down, log_table_size, Pow};

/// From range [i,j), get values (k,j-2^k) where k is the offset
/// into the TwoD table to look up the value for [i,i+2^k) and [j-2^k,j)
/// from which we can get the RMQ.
pub fn adjusted_index(i: Idx, j: Idx) -> (Pow, Idx) {
    let k = log2_down(j - i);
    (k, j - k.value())
}

/// A rather simple 2D array made from vectors of vectors.
/// There are better solutions, but I can implement those later
/// with the same interface.
pub struct TwoD {
    table: Vec<Vec<Idx>>,
}

impl TwoD {
    pub fn new(n: Idx) -> TwoD {
        let Pow(logn) = log_table_size(n);
        let table = vec![vec![0; logn]; n];
        TwoD { table }
    }
}

impl std::ops::Index<(Idx, Pow)> for TwoD {
    type Output = Idx;
    fn index(&self, index: (Idx, Pow)) -> &Self::Output {
        match index {
            (i, Pow(k)) => &self.table[i][k],
        }
    }
}

impl std::ops::IndexMut<(Idx, Pow)> for TwoD {
    fn index_mut(&mut self, index: (Idx, Pow)) -> &mut Self::Output {
        match index {
            (i, Pow(k)) => &mut self.table[i][k],
        }
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
    fn test_adjusted_index() {
        // [0,0) is undefined; j must be larger than i
        // [0,1) => offset=1=2^0, k=0, ii=1-2^k=0
        let (Pow(k), ii) = adjusted_index(0, 1);
        assert_eq!(k, 0);
        assert_eq!(ii, 0);

        // [0,2) => offset=2=2^1, k=1, ii=2-2^1=0
        let (Pow(k), ii) = adjusted_index(0, 2);
        assert_eq!(k, 1);
        assert_eq!(ii, 0);

        // [0,3) => offset=2, k=1 -- second offset; then ii=1
        let (Pow(k), ii) = adjusted_index(0, 3);
        assert_eq!(k, 1);
        assert_eq!(ii, 1);

        // [0,4) => offset=4=2^2, k=2, ii=4-4=0
        let (Pow(k), ii) = adjusted_index(0, 4);
        assert_eq!(k, 2);
        assert_eq!(ii, 0);

        let (Pow(k), ii) = adjusted_index(0, 5);
        assert_eq!(k, 2);
        assert_eq!(ii, 1);

        let (Pow(k), ii) = adjusted_index(0, 6);
        assert_eq!(k, 2);
        assert_eq!(ii, 2);

        let (Pow(k), ii) = adjusted_index(0, 7);
        assert_eq!(k, 2);
        assert_eq!(ii, 3);

        let (Pow(k), ii) = adjusted_index(0, 8);
        assert_eq!(k, 3);
        assert_eq!(ii, 0);

        let (Pow(k), ii) = adjusted_index(1, 8);
        assert_eq!(k, 2);
        assert_eq!(ii, 4);

        let (Pow(k), ii) = adjusted_index(1, 9);
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
