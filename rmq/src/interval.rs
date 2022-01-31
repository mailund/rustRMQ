use super::cmp::Min;

/// The type we use for indexing into our arrays.
pub type Idx = usize;
/// The type our arrays hold. 32-bit are enough for genomic data.
pub type Val = u32;

/// A point is an index with the corresponding value
#[derive(Debug, Clone, Copy)]
pub struct Point(pub Idx, pub Val);

impl Point {
    #[inline]
    pub fn idx(&self) -> Idx {
        self.0
    }
    #[inline]
    pub fn val(&self) -> Val {
        self.1
    }

    pub fn new(i: Idx, x: &[Val]) -> Point {
        Point(i, x[i])
    }
}

impl Min for Point {
    #[inline]
    fn min(p1: Point, p2: Point) -> Point {
        if p1.idx() > p2.idx() {
            // min is symmetric, giving preference to the smallest index,
            // so if p2 has the smallest index, we flip the points.
            return Self::min(p2, p1);
        }
        // Pick the smallest value, but in case of ties, pick p1.
        match p1.val() <= p2.val() {
            true => p1,
            false => p2,
        }
    }
}

impl std::fmt::Display for Point {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "Point({},{})", self.0, self.1)
    }
}

/// Finds the left-most index with the smallest value in x.
/// Returns the index of the left-most minimal value and the
/// minimal value. If [i,j) is not a valid interval, you ged
/// None.
pub fn smallest_in_range(x: &[Val], i: Idx, j: Idx) -> Option<Point> {
    let y = &x[i..j];
    let min_val = y.iter().min()?;
    let pos = i + y.iter().position(|a| a == min_val)?;
    Some(Point(pos, *min_val))
}
