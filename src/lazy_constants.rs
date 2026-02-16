use crate::Num;
use std::any::TypeId;
use std::collections::HashMap;
use std::sync::{OnceLock, RwLock};

enum Table {
    StaticF64(&'static [f64]),
    FVec(TypeId, Vec<u8>),
}
static CONSTANTS_TABLE_HASHMAP: OnceLock<RwLock<HashMap<(usize, TypeId), Table>>> = OnceLock::new();

fn create_f_table<F: Num>(table: &'static [f64]) -> Table {
    let x: Vec<F> = table
        .iter()
        .map(|f| {
            let f: F = (*f as f32).into();
            f
        })
        .collect();
    Table::FVec(TypeId::of::<F>(), unsafe {
        std::mem::transmute::<Vec<F>, Vec<u8>>(x)
    })
}
fn invoke_f<F: Num, R, FN: FnMut(&[F]) -> R>(mut f: FN, table: &Table) -> R {
    match table {
        Table::StaticF64(table) => {
            assert!(TypeId::of::<F>() == TypeId::of::<f64>());
            f(unsafe { std::mem::transmute::<&[f64], &[F]>(table) })
        }
        Table::FVec(t, table) => {
            assert!(TypeId::of::<F>() == *t);
            return f(unsafe { std::mem::transmute::<&[u8], &[F]>(table) });
        }
    }
}

fn add_constants_table<F: Num>(table: &'static [f64]) {
    let table_address: *const [f64] = table;
    let table_address = table_address.addr();
    let f_type = TypeId::of::<F>();
    let key = (table_address, f_type);
    let map = CONSTANTS_TABLE_HASHMAP.get_or_init(|| RwLock::new(HashMap::new()));
    let mut map_write = map.write().unwrap();
    map_write
        .entry(key)
        .or_insert_with(|| create_f_table::<F>(table));
}
pub fn use_constants_table<F: Num, R, FN: FnMut(&[F]) -> R>(f: FN, table: &'static [f64]) -> R {
    let f_type = TypeId::of::<F>();
    let f64_type = TypeId::of::<f64>();
    if f_type == f64_type {
        invoke_f(f, &Table::StaticF64(table))
    } else {
        let table_address: *const [f64] = table;
        let table_address = table_address.addr();
        let key = &(table_address, f_type);
        let map = CONSTANTS_TABLE_HASHMAP.get_or_init(|| RwLock::new(HashMap::new()));
        let map_read = map.read().unwrap();
        if let Some(table) = map_read.get(key) {
            return invoke_f(f, table);
        }
        drop(map_read);
        add_constants_table::<F>(table);
        use_constants_table(f, table)
    }
}

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
    assert!(use_constants_table(
        assert_eq_value(0, 5.0),
        &[5., 6., 7., 8.],
    ));
}

#[test]
fn test_table_f32() {
    assert!(use_constants_table(
        assert_eq_value(0, 5.0),
        &[5., 6., 7., 8.],
    ));
}

#[test]
fn test_table_rational() {
    assert!(use_constants_table(
        assert_eq_value::<crate::bignum::RationalN<2>>(0, (5, 1).into()),
        &[5., 6., 7., 8.],
    ));
}
