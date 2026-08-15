//! C ABI. Prefix `lc_`. Caller owns every buffer.

#![deny(unsafe_op_in_unsafe_fn)]

use std::cell::RefCell;
use std::ffi::{c_char, c_int, CString};
use std::ptr;

use crate::{Cell, Error};

/// Periodic parallelepiped. Lattice vectors are a, b, c (same as vesin rows).
#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct lc_cell {
    /// Lattice vector a, x.
    pub ax: f64,
    /// Lattice vector a, y.
    pub ay: f64,
    /// Lattice vector a, z.
    pub az: f64,
    /// Lattice vector b, x.
    pub bx: f64,
    /// Lattice vector b, y.
    pub by: f64,
    /// Lattice vector b, z.
    pub bz: f64,
    /// Lattice vector c, x.
    pub cx: f64,
    /// Lattice vector c, y.
    pub cy: f64,
    /// Lattice vector c, z.
    pub cz: f64,
    /// Origin x.
    pub ox: f64,
    /// Origin y.
    pub oy: f64,
    /// Origin z.
    pub oz: f64,
}

thread_local! {
    static LAST_ERROR: RefCell<Option<CString>> = RefCell::new(None);
}

fn set_error(msg: &str) {
    LAST_ERROR.with(|slot| {
        *slot.borrow_mut() = CString::new(msg).ok();
    });
}

fn clear_error() {
    LAST_ERROR.with(|slot| {
        *slot.borrow_mut() = None;
    });
}

fn fail(err: Error) -> c_int {
    set_error(&err.to_string());
    1
}

fn fail_msg(msg: &str) -> c_int {
    set_error(msg);
    1
}

/// Thread-local last-error string from this thread's most recent `lc_*` call.
///
/// Returns a pointer to a NUL-terminated UTF-8 C string, or `NULL` if the
/// last call on this thread succeeded (or no call has failed yet).
/// The string is thread-local: distinct threads have independent slots.
/// The pointer is valid until the next `lc_*` call on this thread.
/// Do not free it.
#[no_mangle]
pub extern "C" fn lc_last_error() -> *const c_char {
    LAST_ERROR.with(|slot| {
        slot.borrow()
            .as_ref()
            .map(|s| s.as_ptr())
            .unwrap_or(ptr::null())
    })
}

/// Library version string. Process-static, NUL-terminated. Do not free.
#[no_mangle]
pub extern "C" fn lc_version() -> *const c_char {
    concat!(env!("CARGO_PKG_VERSION"), "\0").as_ptr() as *const c_char
}

/// k-nearest neighbours for `n` points.
///
/// # Arguments
///
/// - `xyz`: packed row-major `n` triples `(x, y, z)` (`n * 3` doubles).
/// - `n`: point count (`size_t`). `0` is [`Error::Empty`].
/// - `simbox`: periodic cell (lattice vectors a, b, c and origin).
/// - `k`: neighbours per source (`size_t`). `0` is [`Error::ZeroK`].
/// - `mask`: `NULL` includes every point; otherwise `n` ints, nonzero to
///   include that point as both source and candidate.
/// - `cell_hint`: target cell edge. Values `<= 0` select the default
///   (3.0 in the box units).
/// - `out_nn`: caller-owned output. Row-major `n * k` neighbour indices
///   (`int`). Neighbours of source `i` occupy `out_nn[i*k .. i*k+k]`,
///   nearest first. Missing slots are `-1`. Length is `n * k`.
///
/// # Null rules
///
/// `xyz`, `simbox`, and `out_nn` must be non-null when `n > 0` and
/// `k > 0`. `mask` may be `NULL`.
///
/// # Buffer layout
///
/// `out_nn` is row-major `n * k`. Unused / missing neighbour slots are
/// the sentinel `-1`.
///
/// # Return
///
/// `0` on success. Nonzero on failure. Read [`lc_last_error`] on the
/// same thread. The last-error string is thread-local: it is not shared
/// across threads, is valid until the next `lc_*` call on this thread,
/// and must not be freed.
///
/// # Safety
///
/// `xyz` is aligned for `f64` and readable for `n * 3` doubles.
/// `simbox` is aligned and points at one valid [`lc_cell`].
/// `out_nn` is aligned for `int` and writable for `n * k` ints.
/// `mask`, if non-null, is aligned for `int` and readable for `n` ints.
/// `n * 3` and `n * k` fit in `usize`.
#[no_mangle]
pub unsafe extern "C" fn lc_knearest(
    xyz: *const f64,
    n: usize,
    simbox: *const lc_cell,
    k: usize,
    mask: *const c_int,
    cell_hint: f64,
    out_nn: *mut c_int,
) -> c_int {
    if n == 0 {
        return fail(Error::Empty);
    }
    if k == 0 {
        return fail(Error::ZeroK);
    }
    if xyz.is_null() || simbox.is_null() || out_nn.is_null() {
        return fail_msg("null pointer");
    }
    let Some(need) = n.checked_mul(k) else {
        return fail(Error::BufferSize);
    };
    let max_xyz = (isize::MAX as usize) / 3;
    let max_out = (isize::MAX as usize) / std::mem::size_of::<c_int>();
    if n > max_xyz || need > max_out {
        return fail(Error::BufferSize);
    };
    // SAFETY: simbox is non-null, aligned, and points at one valid lc_cell.
    let box_c = unsafe { *simbox };
    let sim = match Cell::from_vectors(
        [box_c.ax, box_c.ay, box_c.az],
        [box_c.bx, box_c.by, box_c.bz],
        [box_c.cx, box_c.cy, box_c.cz],
        [box_c.ox, box_c.oy, box_c.oz],
    ) {
        Ok(b) => b,
        Err(e) => {
            set_error(&e.to_string());
            return 1;
        }
    };
    // SAFETY: xyz is non-null, 8-byte aligned, and readable for n packed
    // xyz triples (n * 3 f64s). [f64; 3] has the same alignment as f64.
    let pts: &[[f64; 3]] = unsafe { std::slice::from_raw_parts(xyz.cast::<[f64; 3]>(), n) };
    let mask_vec: Option<Vec<bool>> = if mask.is_null() {
        None
    } else {
        // SAFETY: mask is non-null, aligned for int, and readable for n ints.
        let raw = unsafe { std::slice::from_raw_parts(mask, n) };
        Some(raw.iter().map(|&v| v != 0).collect())
    };
    let hint = if cell_hint > 0.0 {
        Some(cell_hint)
    } else {
        None
    };
    // SAFETY: out_nn is non-null, aligned for int, and writable for
    // n*k ints. n*k fits in usize.
    let out = unsafe { std::slice::from_raw_parts_mut(out_nn, need) };
    if let Err(e) = crate::knearest_into(pts, &sim, k, mask_vec.as_deref(), hint, out) {
        set_error(&e.to_string());
        return 1;
    }
    clear_error();
    0
}
