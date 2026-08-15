//! Recoverable failures from a neighbour search.

use std::fmt;

/// Why a search could not run.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Error {
    /// `k == 0`.
    ZeroK,
    /// A box length is not strictly positive.
    BadBox,
    /// `xyz` is empty.
    Empty,
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Error::ZeroK => write!(f, "k must be at least 1"),
            Error::BadBox => write!(f, "box lengths must be positive"),
            Error::Empty => write!(f, "no points"),
        }
    }
}

impl std::error::Error for Error {}
