mod brute_force;
mod cmp;
mod full_table;
mod inter_table;
mod interval;
mod math;
mod power_table;
mod powers;
mod reduce;
mod rmq_array;

#[allow(unused_imports)] // Importing RMQArray to expose it; I don't use it here
use rmq_array::RMQArray;
#[allow(unused_imports)] // Importing RMQArray to expose it; I don't use it here
use rmq_array::RMQArrayImpl;
use rmq_array::RMQArray_;

/// Implements RMQ by running through the [i,j) interval
/// and finding the index with the smallest value.
/// O(1) preprocessing and O(n) access.
pub type BruteForceRMQ = RMQArray_<brute_force::BruteForceRMQImpl>;

/// Implements RMQ by table lookup. Has a complete table of all [i,j),
/// so uses O(n²) memory, takes O(n²) time preprocessing, but then
/// does RMQ in O(1).
pub type FullTabulateRMQ = RMQArray_<full_table::FullTabulateRMQImpl>;

/// RMQ table that tabulates all [i,i+2^k] ranges (there are O(n log n)),
/// form which we can get the RMQ from the table by splitting [i,j) into
/// two, [i,2^k) and [j-2^k,j) (where k is the largest such k). We can get
/// the RMQ from those two intervals with a table lookup in O(1) and then
/// pick the one of those with the smallest value.
/// The result is O(n log n) preprocessing and O(1) lookup.
pub type PowerRMQ = RMQArray_<powers::PowerRMQImpl>;

/// Reduces the input vector to one of length m = n/log(n); then preprocess
/// it with the PowerRMQ method (in time m log m = n/log(n) log(n/log(n)) = O(n)).
/// Lookup is now one constant time lookup in the PowerRMQ and two linear searches
/// in blocks of length log n. So, preprocessing O(n) and lookup O(log n).
pub type ReducedPowerRMQ = RMQArray_<reduce::ReducedPowerRMQImpl>;

#[cfg(test)]
mod tests {
    use super::interval::{Idx, Point, Val};
    use super::*;

    fn check_same_index<R: RMQArray>(rmqa: &R, vals: &Vec<u32>) {
        assert_eq!(rmqa.len(), vals.len());
        for (i, &v) in vals.iter().enumerate() {
            assert_eq!(v, rmqa[i]);
        }
    }

    fn check_min_in_interval<R: RMQArray>(rmqa: &R, i: Idx, j: Idx) {
        let Point(k, _) = rmqa.rmq(i, j);
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
        // Not power of two
        let v = vec![2, 1, 2, 5, 3, 6, 1, 3, 7, 4];
        let rmqa = R::new(v.clone());
        check_same_index(&rmqa, &v);
        // Power of two
        let v = vec![2, 1, 2, 5, 3, 6, 1, 3, 7, 4, 2, 6, 3, 4, 7, 9];
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

    #[test]
    fn test_index_reduced_power() {
        check_index::<ReducedPowerRMQ>()
    }

    fn check_rmq<R: RMQArray>() {
        // Not power of two
        let v = vec![2, 1, 2, 5, 3, 6, 1, 3, 7, 4];
        let rmqa = R::new(v.clone());
        check_min(&rmqa);
        // Power of two
        let v = vec![2, 1, 2, 5, 3, 6, 1, 3, 7, 4, 2, 6, 3, 4, 7, 9];
        let rmqa = R::new(v.clone());
        check_min(&rmqa);
        // Not power of two
        let v = vec![2, 1, 2, 0, 2, 1, 3, 7, 4];
        let rmqa = R::new(v.clone());
        check_min(&rmqa);
        // Power of two
        let v = vec![2, 1, 2, 5, 3, 6, 1, 3];
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
        check_rmq::<PowerRMQ>()
    }

    #[test]
    fn test_rmq_reduced_power() {
        check_rmq::<ReducedPowerRMQ>()
    }

    // Not really RMQ but trying to use them on lcp-intervals
    fn next_lidx<RMQ: RMQArray>(lcp: &RMQ, i: Idx, j: Idx, l: Val) -> Option<Idx> {
        match (i, j) {
            (i, j) if i == j => None,        // empty interval at the end...
            (i, j) if j == i + 1 => Some(j), // last interval
            (i, j) if j > i + 1 => match lcp.rmq(i + 1, j) {
                Point(ii, ll) if l == ll => Some(ii),
                _ => Some(j), // We reached the end
            },
            _ => panic!("Interval with end point before start point"),
        }
    }
    fn children<RMQ: RMQArray>(lcp: &RMQ, i: Idx, j: Idx) -> (Val, Vec<(Idx, Idx)>) {
        let mut res: Vec<(Idx, Idx)> = Vec::new();
        let Point(mut prev, l) = lcp.rmq(i + 1, j);
        res.push((i, prev));
        while let Some(ii) = next_lidx(lcp, prev, j, l) {
            res.push((prev, ii));
            prev = ii;
        }
        (l, res)
    }
    fn widest(v: &Vec<(Idx, Idx)>) -> (Idx, Idx) {
        let mut widest = (0, 0);
        for (i, j) in v {
            if (j - i) > (widest.1 - widest.0) {
                widest = (*i, *j);
            }
        }
        widest
    }

    // Collecting the arguments in a struct...
    struct BTR<'a, RMQ: RMQArray> {
        sa: &'a Vec<Idx>,
        isa: &'a Vec<Idx>,
        lcp: &'a RMQ,
    }

    fn branch_tr_rec<RMQ: RMQArray>(btr: &BTR<RMQ>, i: Idx, j: Idx, res: &mut Vec<(Idx, Val)>) {
        let (l, sub_intervals) = children(btr.lcp, i, j);
        let offset = l as Idx;
        let (wi, wj) = widest(&sub_intervals);
        for (ii, jj) in sub_intervals {
            if (ii, jj) == (wi, wj) {
                // skip widest interval
                continue;
            }
            if jj - ii > 1 {
                // sub interval rather than leaf
                branch_tr_rec(btr, ii, jj, res);
            }
            for q in ii..jj {
                if btr.sa[q] + offset < btr.sa.len() {
                    let r = btr.isa[btr.sa[q] + offset];
                    if (i <= r && r < ii) || (jj <= r && r < j) {
                        res.push((btr.sa[q], 2 * l))
                    }
                }
                if offset < btr.sa[q] {
                    let r = btr.isa[btr.sa[q] - offset];
                    if wi <= r && r < wj {
                        res.push((btr.sa[r], 2 * l))
                    }
                }
            }
        }
    }
    fn branch_tr<RMQ: RMQArray>(btr: &BTR<RMQ>, i: Idx, j: Idx) -> Vec<(Idx, Val)> {
        let mut res: Vec<(Idx, Val)> = Vec::new();
        branch_tr_rec(btr, i, j, &mut res);
        res
    }

    #[test]
    fn mississippi() {
        let x = "mississippi$";
        let sa: Vec<Idx> = vec![11, 10, 7, 4, 1, 0, 9, 8, 6, 3, 5, 2];
        let mut isa: Vec<Idx> = vec![0; sa.len()];
        for (i, ii) in sa.iter().enumerate() {
            isa[*ii] = i;
        }
        let lcp = PowerRMQ::new(vec![0, 0, 1, 1, 4, 0, 0, 1, 0, 2, 1, 3]);
        let btr = BTR {
            sa: &sa,
            isa: &isa,
            lcp: &lcp,
        };
        for (i, l) in branch_tr(&btr, 0, lcp.len()).iter() {
            println!(
                "Branching tandem repeat {} at {}",
                &x[*i..*i + *l as Idx],
                i
            );
        }
        //assert!(false);
    }
}
