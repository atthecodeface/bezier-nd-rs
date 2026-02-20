use crate::Num;
use std::any::{Any, TypeId};
use std::collections::HashMap;
use std::sync::{OnceLock, RwLock};

struct TableEntry {
    type_id: TypeId,
    table: Box<dyn Any + Send + Sync>,
}

impl TableEntry {
    /// Create a [Table] whose contents is a `Vec<F>` from a slice of `[T]`
    fn new_f_table<F: Num, T: Copy, C: Fn(T) -> F>(convert: C, table: &[T]) -> Self {
        let vec_f: Vec<F> = table.iter().copied().map(convert).collect();
        Self {
            type_id: TypeId::of::<F>(),
            table: Box::new(vec_f),
        }
    }

    /// Invoke a function of `&[F]` on a table whose contents are known to be type `F`
    fn invoke_f<F: Num, R, FN: FnMut(&[F]) -> R>(&self, mut f: FN) -> R {
        assert!(self.type_id == TypeId::of::<F>());
        return f(self.table.downcast_ref::<Vec<F>>().unwrap());
    }

    /// Invoke a function of `&[F]` where `F==f64` on a table whose contents are `&'static [f64]`
    ///
    /// Note that `Any` requires the table to be *static*
    fn invoke_f_of_ts<F: Num, T: Num, R, FN: FnMut(&[F]) -> R>(
        mut f: FN,
        table: &&'static [T],
    ) -> R {
        assert!(TypeId::of::<F>() == TypeId::of::<T>());
        let table = (table as &dyn Any).downcast_ref::<&[F]>().unwrap();
        f(table)
    }
}

/// A static [ConstantsTable] using a [OnceLock]
pub struct StaticConstantsTable(OnceLock<SharedConstantsTable>);

/// A global [ConstantsTable] that can be used by multiple threads, to invoke functions
/// on tables of any type `F:Num` where there is a table provided of (e.g.) type `&'static [f64]`
pub static CONSTANTS_TABLE: StaticConstantsTable = StaticConstantsTable::new();

impl StaticConstantsTable {
    /// Create a new [ConstantsTable]
    ///
    /// This is normally done statically with:
    ///
    /// ```text
    /// static CONSTANTS_TABLE: StaticConstantsTable = StaticConstantsTable::new();
    /// ```
    pub const fn new() -> Self {
        Self(OnceLock::new())
    }

    /// Invoke a function `f` on a slice of `[F]` derived from a *static* slice of `[f64]`
    ///
    /// If `F` is `f64` then this is performed without any conversions
    ///
    /// If `F` is a different type (supporting `Num`, hence providing and `of_f64`
    /// conversion method) then an intermediate `Vec` must be allocated, and filled
    /// with the converted values. This conversion is performed *once* when the lazy constants table is used.
    pub fn use_constants_table<F: Num, R, FN: FnMut(&[F]) -> R>(
        &self,
        f: FN,
        table: &'static [f64],
    ) -> R {
        let f_type = TypeId::of::<F>();
        let f64_type = TypeId::of::<f64>();
        if f_type == f64_type {
            TableEntry::invoke_f_of_ts(f, &table)
        } else {
            let constants_table = self.0.get_or_init(SharedConstantsTable::new);
            constants_table.invoke_using_table(f, table)
        }
    }
}

/// A trait supported by constants tables that allow invocation of a mapping
/// function given a table of type `&[T]`
pub trait ConstantsTable<T: Num> {
    /// Invoke the mapping function `map` on the table, having converted the mapping
    /// first into something that permits a `&[F]`
    fn invoke_using_table<F: Num, R, M: FnMut(&[F]) -> R>(&self, map: M, table: &[T]) -> R;

    /// Attempt to cache a constants table after conversion, given a static version thereof
    ///
    /// Returns `false` if the cache addition could not be performed (a non-caching table, for example)
    fn cache_constants_table<F: Num>(&self, _table: &'static [T]) -> bool {
        false
    }

    /// Invoke a function `f` on a slice of `[F]` derived from a *static* slice of `[T]`
    ///
    /// If `F` is `T` then this is performed without any conversions
    ///
    /// If `F` is a different type then an intermediate `Vec` must be allocated, and filled
    /// with the converted values. This conversion may be cached in the table.
    fn use_constants_table<F: Num, R, M: FnMut(&[F]) -> R>(
        &self,
        map: M,
        table: &'static [T],
    ) -> R {
        let f_type = TypeId::of::<F>();
        let t_type = TypeId::of::<T>();
        if f_type == t_type {
            TableEntry::invoke_f_of_ts(map, &table)
        } else {
            self.invoke_using_table(map, table)
        }
    }
}

/// A table of constants tables (each of type `[F]` for some `F:Num`), accessed through a key
/// based on the [TypeId] of `F` and the *address* of the initial contents of the table
/// (which must be a `&'static[T]`)
#[derive(Default)]
pub struct SharedConstantsTable {
    table: RwLock<HashMap<(usize, TypeId), TableEntry>>,
}

impl SharedConstantsTable {
    /// Create a new [ConstantsTable]
    ///
    /// This cannot be `const` as [HashMap] does not have a `const` new method.
    pub fn new() -> Self {
        Self {
            table: RwLock::new(HashMap::new()),
        }
    }
}

macro_rules! ConstantsTableImpl {
    {$t:ty, $conv:ident} => {
        impl ConstantsTable<$t> for SharedConstantsTable {
            /// Invoke a function on a table in the set
            fn invoke_using_table<F: Num, R, M: FnMut(&[F]) -> R>(&self, f: M, table: &[$t]) -> R {
                let f_type = TypeId::of::<F>();
                let table_address: *const [$t] = table;
                let table_address = table_address.addr();
                let key = (table_address, f_type);
                eprintln!("Here");
                let rd = self.table.read().unwrap();
                let Some(table) = rd.get(&key) else {
                    drop(rd);
                    self.table
                        .write()
                        .unwrap()
                        .entry(key)
                        .or_insert_with(|| TableEntry::new_f_table(F :: $conv, table));
                    return self.table.read().unwrap().get(&key).unwrap().invoke_f(f);
                };
                table.invoke_f(f)
            }

            fn cache_constants_table<F: Num>(&self, table: &'static [$t]) -> bool {
                let f_type = TypeId::of::<F>();
                let table_address: *const [$t] = table;
                let table_address = table_address.addr();
                let key = (table_address, f_type);

                self.table
                    .write()
                    .unwrap()
                    .entry(key)
                    .or_insert_with(|| TableEntry::new_f_table(F :: $conv, table));
                true
            }
        }
    };
}
ConstantsTableImpl! {f64, of_f64}

#[cfg(test)]
fn assert_eq_value<F: Num>(index: usize, value: F) -> impl Fn(&[F]) -> bool {
    move |table| {
        assert_eq!(
            table[index], value,
            "Expected table index {index} to match value: {table:?}"
        );
        true
    }
}
#[test]
fn test_table_f64() {
    assert!(CONSTANTS_TABLE.use_constants_table(assert_eq_value(0, 5.0), &[5., 6., 7., 8.],));
}

#[test]
fn test_table_f32() {
    assert!(CONSTANTS_TABLE.use_constants_table(assert_eq_value(0, 5.0), &[5., 6., 7., 8.],));
}

#[test]
fn test_table_rational() {
    assert!(CONSTANTS_TABLE.use_constants_table(
        assert_eq_value::<crate::bignum::RationalN<2>>(0, (5, 1).into()),
        &[5., 6., 7., 8.],
    ));
}
