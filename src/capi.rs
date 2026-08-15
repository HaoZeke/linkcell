//! C ABI. Prefix `lc_`. Caller owns every buffer.

use std::ffi::{c_char, c_int, CString};
use std::ptr;
use std::sync::Mutex;

use crate::{knearest, OrthoBox};

/// Orthorhombic periodic box. Lengths must be positive.
#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct lc_box {
    /// Edge length x.
    pub lx: f64,
    /// Edge length y.
    pub ly: f64,
    /// Edge length z.
    pub lz: f64,
    /// Origin x of the dump cell.
    pub xlo: f64,
    /// Origin y of the dump cell.
    pub ylo: f64,
    /// Origin z of the dump cell.
    pub zlo: f64,
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
    simbox: *const lc_box,
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
    let sim = match OrthoBox::with_origin(
        box_c.lx, box_c.ly, box_c.lz, box_c.xlo, box_c.ylo, box_c.zlo,
    ) {
        Ok(b) => b,
        Err(e) => {
            set_error(&e.to_string());
            return 1;
        }
    };
    let pts: Vec<[f64; 3]> = std::slice::from_raw_parts(xyz, n_us * 3)
        .chunks_exact(3)
        .map(|c| [c[0], c[1], c[2]])
        .collect();
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
    let rows = match knearest(
        &pts,
        &sim,
        k_us,
        mask_vec.as_deref(),
        hint,
    ) {
        Ok(r) => r,
        Err(e) => {
            set_error(&e.to_string());
            return 1;
        }
    };
    let out = std::slice::from_raw_parts_mut(out_nn, n_us * k_us);
    for slot in out.iter_mut() {
        *slot = -1;
    }
    for (i, row) in rows.iter().enumerate() {
        for (t, &j) in row.indices.iter().enumerate() {
            if t < k_us {
                out[i * k_us + t] = j as c_int;
            }
        }
    }
    if let Ok(mut slot) = LAST_ERROR.lock() {
        *slot = None;
    }
    0
}
