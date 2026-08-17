//! CUDA path: DLPack device pointers into `lc_gpu_*`.

use std::ffi::{c_char, c_int, c_void, CStr};
use std::ptr;

use dlpk::sys::{
    DLDataType, DLDataTypeCode, DLDevice, DLDeviceType, DLManagedTensorVersioned, DLPackVersion,
    DLTensor, DLPACK_FLAG_BITMASK_IS_COPIED, DLPACK_MAJOR_VERSION, DLPACK_MINOR_VERSION,
};
use linkcell::lc_cell;
use pyo3::exceptions::{PyRuntimeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyCapsule, PyTuple};

use crate::{capsule_of, numel, with_dltensor};

#[repr(C)]
struct lc_gpu_workspace {
    _private: [u8; 0],
}

unsafe extern "C" {
    fn lc_gpu_available() -> c_int;
    fn lc_gpu_last_error() -> *const c_char;
    fn lc_gpu_workspace_new() -> *mut lc_gpu_workspace;
    fn lc_gpu_workspace_free(ws: *mut lc_gpu_workspace);
    fn lc_gpu_knearest(
        ws: *mut lc_gpu_workspace,
        xyz: *const f64,
        n: usize,
        simbox: *const lc_cell,
        k: usize,
        mask: *const c_int,
        cell_hint: f64,
        out_nn: *mut c_int,
    ) -> c_int;
    fn lc_gpu_queue(ws: *mut lc_gpu_workspace) -> *mut c_void;
    fn lc_gpu_alloc(ptr: *mut *mut c_void, bytes: usize) -> c_int;
    fn lc_gpu_free(ptr: *mut c_void);
    fn lc_gpu_fill_i32(ptr: *mut c_void, value: c_int, n: usize) -> c_int;
}

fn gpu_err() -> PyErr {
    let p = unsafe { lc_gpu_last_error() };
    if p.is_null() {
        PyRuntimeError::new_err("gpu knearest failed")
    } else {
        let s = unsafe { CStr::from_ptr(p) }.to_string_lossy().into_owned();
        PyRuntimeError::new_err(s)
    }
}

pub fn available() -> bool {
    unsafe { lc_gpu_available() != 0 }
}

fn capsule_on_stream<'py>(
    obj: &Bound<'py, PyAny>,
    stream: *mut c_void,
) -> PyResult<Bound<'py, PyCapsule>> {
    if !stream.is_null() {
        let kwargs = pyo3::types::PyDict::new(obj.py());
        kwargs.set_item("stream", stream as usize)?;
        if let Ok(cap) = obj.call_method("__dlpack__", (), Some(&kwargs)) {
            if let Ok(c) = cap.cast_into::<PyCapsule>() {
                return Ok(c);
            }
        }
    }
    capsule_of(obj)
}

struct CudaView {
    ptr: *const u8,
    shape: Vec<i64>,
    device_id: i32,
}

fn take_cuda_f64<'py>(
    obj: &Bound<'py, PyAny>,
    stream: *mut c_void,
) -> PyResult<(Bound<'py, PyCapsule>, CudaView)> {
    let cap = capsule_on_stream(obj, stream)?;
    let view = with_dltensor(&cap, |t| {
        if t.device.device_type != DLDeviceType::kDLCUDA {
            return Err(PyValueError::new_err("expected kDLCUDA xyz"));
        }
        if t.dtype.code != DLDataTypeCode::kDLFloat || t.dtype.bits != 64 {
            return Err(PyValueError::new_err("CUDA xyz must be float64"));
        }
        let shape = if t.ndim <= 0 || t.shape.is_null() {
            return Err(PyValueError::new_err("null CUDA shape"));
        } else {
            unsafe { std::slice::from_raw_parts(t.shape, t.ndim as usize) }.to_vec()
        };
        if !t.strides.is_null() {
            let strides = unsafe { std::slice::from_raw_parts(t.strides, shape.len()) };
            let mut expect = 1i64;
            for (&dim, &st) in shape.iter().rev().zip(strides.iter().rev()) {
                if dim > 1 && st != expect {
                    return Err(PyValueError::new_err("CUDA xyz must be C-contiguous"));
                }
                if dim > 1 {
                    expect = expect.saturating_mul(dim);
                }
            }
        }
        let ptr = t.data.wrapping_add(t.byte_offset as usize).cast::<u8>();
        Ok(CudaView {
            ptr: ptr.cast_const(),
            shape,
            device_id: t.device.device_id,
        })
    })?;
    Ok((cap, view))
}

fn xyz_n(shape: &[i64]) -> PyResult<usize> {
    match shape {
        [n, 3] if *n > 0 => Ok(*n as usize),
        [n] if *n > 0 && *n % 3 == 0 => Ok((*n as usize) / 3),
        _ => Err(PyValueError::new_err(
            "xyz must be float64 shape (n, 3) or (n*3,)",
        )),
    }
}

struct CudaOwner {
    ptr: *mut c_void,
    shape: [i64; 2],
    strides: [i64; 2],
}

unsafe extern "C" fn free_dlpack_capsule(capsule: *mut pyo3::ffi::PyObject) {
    let ptr = pyo3::ffi::PyCapsule_GetPointer(capsule, c"dltensor_versioned".as_ptr());
    if ptr.is_null() {
        return;
    }
    delete_cuda_i32(ptr.cast());
}

unsafe extern "C" fn delete_cuda_i32(managed: *mut DLManagedTensorVersioned) {
    if managed.is_null() {
        return;
    }
    let ctx = (*managed).manager_ctx as *mut CudaOwner;
    if !ctx.is_null() {
        lc_gpu_free((*ctx).ptr);
        drop(Box::from_raw(ctx));
    }
    drop(Box::from_raw(managed));
}

/// CUDA int32 neighbour buffer. `torch.from_dlpack` / `cupy.from_dlpack`.
#[pyclass]
pub struct DlpackCuda {
    ptr: usize,
    n: usize,
    k: usize,
    device_id: i32,
    exported: std::sync::atomic::AtomicBool,
}

impl Drop for DlpackCuda {
    fn drop(&mut self) {
        if !self.exported.load(std::sync::atomic::Ordering::SeqCst) && self.ptr != 0 {
            unsafe { lc_gpu_free(self.ptr as *mut c_void) };
            self.ptr = 0;
        }
    }
}

#[pymethods]
impl DlpackCuda {
    #[pyo3(signature = (*, stream = None, max_version = None, dl_device = None, copy = None))]
    fn __dlpack__<'py>(
        &self,
        py: Python<'py>,
        stream: Option<Bound<'py, PyAny>>,
        max_version: Option<Bound<'py, PyAny>>,
        dl_device: Option<Bound<'py, PyAny>>,
        copy: Option<Bound<'py, PyAny>>,
    ) -> PyResult<Bound<'py, PyCapsule>> {
        let _ = (stream, max_version, copy);
        if self.exported.load(std::sync::atomic::Ordering::SeqCst) {
            return Err(PyValueError::new_err("dlpack capsule already used"));
        }
        if let Some(dev) = dl_device {
            let want = self.__dlpack_device__(py)?;
            if dev.ne(want.bind(py))? {
                return Err(PyValueError::new_err("unsupported dl_device"));
            }
        }
        let owner = Box::new(CudaOwner {
            ptr: self.ptr as *mut c_void,
            shape: [self.n as i64, self.k as i64],
            strides: [self.k as i64, 1],
        });
        let owner_ptr = Box::into_raw(owner);
        let managed = Box::new(DLManagedTensorVersioned {
            version: DLPackVersion {
                major: DLPACK_MAJOR_VERSION,
                minor: DLPACK_MINOR_VERSION,
            },
            manager_ctx: owner_ptr.cast(),
            deleter: Some(delete_cuda_i32),
            flags: DLPACK_FLAG_BITMASK_IS_COPIED,
            dl_tensor: DLTensor {
                data: self.ptr as *mut c_void,
                device: DLDevice {
                    device_type: DLDeviceType::kDLCUDA,
                    device_id: self.device_id,
                },
                ndim: 2,
                dtype: DLDataType {
                    code: DLDataTypeCode::kDLInt,
                    bits: 32,
                    lanes: 1,
                },
                shape: unsafe { &mut (*owner_ptr).shape }.as_mut_ptr(),
                strides: unsafe { &mut (*owner_ptr).strides }.as_mut_ptr(),
                byte_offset: 0,
            },
        });
        let raw = Box::into_raw(managed);
        self.exported
            .store(true, std::sync::atomic::Ordering::SeqCst);
        let ptr = std::ptr::NonNull::new(raw.cast())
            .ok_or_else(|| PyRuntimeError::new_err("null dlpack tensor"))?;
        unsafe {
            PyCapsule::new_with_pointer_and_destructor(
                py,
                ptr,
                c"dltensor_versioned",
                Some(free_dlpack_capsule),
            )
        }
    }

    fn __dlpack_device__<'py>(&self, py: Python<'py>) -> PyResult<Py<PyTuple>> {
        Ok(PyTuple::new(py, [2i32, self.device_id])?.unbind())
    }
}

pub fn knearest_cuda<'py>(
    py: Python<'py>,
    xyz: &Bound<'py, PyAny>,
    sim: &lc_cell,
    k: usize,
    mask: Option<&Bound<'py, PyAny>>,
    cell_hint: Option<f64>,
) -> PyResult<Py<DlpackCuda>> {
    if !available() {
        return Err(PyRuntimeError::new_err(
            "CUDA xyz but the driver or nvrtc is not loaded",
        ));
    }
    let ws = unsafe { lc_gpu_workspace_new() };
    if ws.is_null() {
        return Err(gpu_err());
    }
    struct Ws(*mut lc_gpu_workspace);
    impl Drop for Ws {
        fn drop(&mut self) {
            if !self.0.is_null() {
                unsafe { lc_gpu_workspace_free(self.0) };
            }
        }
    }
    let ws = Ws(ws);
    let stream = unsafe { lc_gpu_queue(ws.0) };
    let (xyz_cap, view) = take_cuda_f64(xyz, stream)?;
    let n = xyz_n(&view.shape)?;
    let need = n
        .checked_mul(k)
        .ok_or_else(|| PyValueError::new_err("n * k overflows"))?;
    let mask_cap;
    let mask_ptr = if let Some(m) = mask {
        let cap = capsule_on_stream(m, stream)?;
        let ptr = with_dltensor(&cap, |t| {
            if t.device.device_type != DLDeviceType::kDLCUDA {
                return Err(PyValueError::new_err(
                    "mask must be on the same CUDA device",
                ));
            }
            if numel(unsafe {
                if t.shape.is_null() {
                    &[]
                } else {
                    std::slice::from_raw_parts(t.shape, t.ndim.max(0) as usize)
                }
            })? != n
            {
                return Err(PyValueError::new_err("mask length must be n"));
            }
            Ok(t.data.wrapping_add(t.byte_offset as usize).cast::<c_int>())
        })?;
        mask_cap = Some(cap);
        ptr
    } else {
        mask_cap = None;
        ptr::null()
    };
    let _keep = (xyz_cap, mask_cap);
    let mut out: *mut c_void = ptr::null_mut();
    if unsafe { lc_gpu_alloc(&mut out, need * std::mem::size_of::<c_int>()) } != 0 {
        return Err(gpu_err());
    }
    if unsafe { lc_gpu_fill_i32(out, -1, need) } != 0 {
        unsafe { lc_gpu_free(out) };
        return Err(gpu_err());
    }
    let hint = cell_hint.unwrap_or(0.0);
    let ws_p = ws.0 as usize;
    let xyz_p = view.ptr as usize;
    let mask_p = mask_ptr as usize;
    let out_p = out as usize;
    let sim_c = *sim;
    let st = py.detach(|| unsafe {
        lc_gpu_knearest(
            ws_p as *mut lc_gpu_workspace,
            xyz_p as *const f64,
            n,
            &sim_c,
            k,
            mask_p as *const c_int,
            hint,
            out_p as *mut c_int,
        )
    });
    if st != 0 {
        unsafe { lc_gpu_free(out) };
        return Err(gpu_err());
    }
    Py::new(
        py,
        DlpackCuda {
            ptr: out as usize,
            n,
            k,
            device_id: view.device_id,
            exported: std::sync::atomic::AtomicBool::new(false),
        },
    )
}
