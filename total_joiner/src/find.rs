use std::ops::{Index, Range};

pub trait Length {
    fn len(&self) -> usize;
}

impl<T> Length for Vec<T> {
    fn len(&self) -> usize {
        self.len()
    }
}

impl<T> Length for [T] {
    fn len(&self) -> usize {
        SliceExt::len(self)
    }
}

pub trait Find<T: PartialEq+?Sized> {
    fn find(&self, to_find: &T) -> Option<usize>;
}

impl<T> Find<<T as Index<Range<usize>>>::Output> for T
    where
        T: Index<Range<usize>> + Length,
        <T as Index<Range<usize>>>::Output: PartialEq + Length,
{
    fn find(&self, to_find: &<T as Index<Range<usize>>>::Output) -> Option<usize> {
        for i in range(0, self.len() - to_find.len() + 1) {
            if to_find == &self[i..i+to_find.len()] {
                return Some(i);
            }
        }
        None
    }
}
