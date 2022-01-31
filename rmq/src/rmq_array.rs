use super::interval::{Idx, Point, Val};

/// To make a RMQArray, you should implement this and then wrap the type
/// in RMQArray_<>.
pub trait RMQArrayImpl {
    /// Create the necessary tables from a vector of values.
    fn new(values: Vec<Val>) -> Self;
    /// Get the length of the underlying array.
    fn len(&self) -> Idx;
    /// Get the value at a specific index.
    fn val(&self, index: Idx) -> &Val;
    /// Get the point at the first minimal value in the range [i,j).
    fn rmq(&self, i: Idx, j: Idx) -> Point;

    /// Get an index as a value
    fn point(&self, i: Idx) -> Point {
        Point(i, *self.val(i))
    }
}

/// Range Minimum Query interface.
/// We are just telling the type system that RMQArray implements RMQArrayImpl
/// and Index<Idx>.
pub trait RMQArray: RMQArrayImpl + std::ops::Index<Idx, Output = Val> {}

/// Wrapper type to extend RMQArrayImpl.
pub struct RMQArray_<Impl: RMQArrayImpl>(Impl);

/// A delegate that pulls the interface from T into RMQArray_<T>.
impl<T> RMQArrayImpl for RMQArray_<T>
where
    T: RMQArrayImpl,
{
    /// Create the necessary tables from a vector of values.
    fn new(values: Vec<Val>) -> Self {
        RMQArray_::<T>(T::new(values))
    }
    /// Get the length of the underlying array.
    fn len(&self) -> Idx {
        self.0.len()
    }
    /// Get the value at a specific index.
    fn val(&self, i: Idx) -> &Val {
        self.0.val(i)
    }
    /// Get the point at the first minimal value in the range [i,j).
    fn rmq(&self, i: Idx, j: Idx) -> Point {
        self.0.rmq(i, j)
    }
}

/// Make sure that an RMQArray_<T> implements RMQArray.
/// There is nothing to implement, we get it all from RMQArrayImpl,
/// but we still need to tell the type system.
impl<T> RMQArray for RMQArray_<T> where T: RMQArrayImpl {}

/// Index for RMQArray_<>. It gives us the Index<> interface for all
/// RMQArrayImpl wrapped in RMQArray_<>.
impl<T> std::ops::Index<Idx> for RMQArray_<T>
where
    T: RMQArrayImpl,
{
    type Output = Val;
    fn index(&self, index: Idx) -> &Val {
        self.0.val(index)
    }
}

// Deref gives us access to the inner workings of T without
// going through the zero'th index in the wrapper.
impl<T> std::ops::Deref for RMQArray_<T>
where
    T: RMQArrayImpl,
{
    type Target = T;
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}
impl<T> std::ops::DerefMut for RMQArray_<T>
where
    T: RMQArrayImpl,
{
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}
