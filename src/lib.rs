//! Periodic linked-cell k-nearest neighbour search.
//!
//! vesin builds cutoff pair lists. nanoflann builds Euclidean KD-trees
//! without a minimum-image convention. This crate is the missing piece:
//! Allen and Tildesley's linked cells, a k-heap per source, expanding
//! shells until the k-th neighbour cannot lie outside the visited cube.
//!
//! The C ABI (`lc_*`) is the hourglass waist. C++ lives in
//! `include/linkcell.hpp` as a RAII header over that ABI.

#![deny(missing_docs)]

mod error;
mod knearest;
mod ortho;

pub use error::Error;
pub use knearest::{knearest, knearest_brute, Neighbors};
pub use ortho::OrthoBox;

#[cfg(feature = "capi")]
mod capi;
#[cfg(feature = "capi")]
pub use capi::*;
