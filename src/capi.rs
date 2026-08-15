//! C ABI. Prefix `lc_`. Caller owns every buffer.

use std::ffi::{c_char, c_int, CString};
use std::ptr;
use std::sync::Mutex;

use crate::Cell;

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

static LAST_ERROR: Mutex<Option<CString>> = Mutex::new(None);

fn set_error(msg: &str) {
    if let Ok(mut slot) = LAST_ERROR.lock() {
        *slot = CString::new(msg).ok();
    }
}

/// Process-static last error. Do not free. NULL if the last call succeeded.
#[no_mangle]
pub extern "C" fn lc_last_error() -> *const c_char {
    match LAST_ERROR.lock() {
        Ok(slot) => slot
            .as_ref()
            .map(|s| s.as_ptr())
            .unwrap_or(ptr::null()),
        Err(_) => ptr::null(),
    }
}

/// Library version string. Process-static. Do not free.
#[no_mangle]
pub extern "C" fn lc_version() -> *const c_char {
    concat!(env!("CARGO_PKG_VERSION"), "\0").as_ptr() as *const c_char
}

/// k-nearest neighbours for `n` points.
///
/// `xyz` is `n` packed xyz triples. `mask` is NULL (all points) or `n`
/// ints, nonzero to include. `cell_hint <= 0` selects the default edge.
/// `out_nn` is caller-owned, length `n * k`; unused slots are `-1`.
///
/// Returns 0 on success, nonzero on failure. Read [`lc_last_error`].
#[no_mangle]
pub unsafe extern "C" fn lc_knearest(
    xyz: *const f64,
    n: c_int,
    simbox: *const lc_cell,
    k: c_int,
    mask: *const c_int,
    cell_hint: f64,
    out_nn: *mut c_int,
) -> c_int {
    if xyz.is_null() || simbox.is_null() || out_nn.is_null() {
        set_error("null pointer");
        return 1;
    }
    if n <= 0 {
        set_error("no points");
        return 1;
    }
    if k <= 0 {
        set_error("k must be at least 1");
        return 1;
    }
    let n_us = n as usize;
    let k_us = k as usize;
    let box_c = *simbox;
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
    // Packed xyz is already [[f64; 3]; n] in memory. Copying it was
    // the extra tax on every d-SEAMS call.
    let pts: &[[f64; 3]] = std::slice::from_raw_parts(xyz.cast::<[f64; 3]>(), n_us);
    let mask_vec: Option<Vec<bool>> = if mask.is_null() {
        None
    } else {
        Some(
            std::slice::from_raw_parts(mask, n_us)
                .iter()
                .map(|&v| v != 0)
                .collect(),
        )
    };
    let hint = if cell_hint > 0.0 {
        Some(cell_hint)
    } else {
        None
    };
    let out = std::slice::from_raw_parts_mut(out_nn, n_us * k_us);
    if let Err(e) = crate::knearest_into(&pts, &sim, k_us, mask_vec.as_deref(), hint, out)
    {
        set_error(&e.to_string());
        return 1;
    }
    if let Ok(mut slot) = LAST_ERROR.lock() {
        *slot = None;
    }
    0
}
