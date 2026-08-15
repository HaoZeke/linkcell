//! Periodic cell: three lattice vectors and a minimum-image convention.
//!
//! Orthorhombic boxes are the diagonal case. Triclinic (and any
//! parallelepiped) use the same fractional wrap. Lattice vectors are
//! stored as the columns of H, so r = H s. The C ABI and vesin pass
//! rows (a, b, c); constructors accept either.

use crate::Error;

/// Periodic parallelepiped.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Cell {
    /// Columns of H: `h[0]` is lattice vector a.
    h: [[f64; 3]; 3],
    /// Inverse of H. Fractional coordinates are `s = Hinv (r - origin)`.
    hinv: [[f64; 3]; 3],
    origin: [f64; 3],
    /// Perpendicular widths |a · n_a| etc., used to size the cell list.
    widths: [f64; 3],
}

impl Cell {
    /// Diagonal box with origin at zero.
    pub fn ortho(lx: f64, ly: f64, lz: f64) -> Result<Self, Error> {
        Self::from_vectors([lx, 0.0, 0.0], [0.0, ly, 0.0], [0.0, 0.0, lz], [0.0, 0.0, 0.0])
    }

    /// Diagonal box with an explicit dump-cell origin.
    pub fn ortho_origin(
        lx: f64,
        ly: f64,
        lz: f64,
        origin: [f64; 3],
    ) -> Result<Self, Error> {
        Self::from_vectors([lx, 0.0, 0.0], [0.0, ly, 0.0], [0.0, 0.0, lz], origin)
    }

    /// Parallelepiped from lattice vectors a, b, c and an origin.
    pub fn from_vectors(
        a: [f64; 3],
        b: [f64; 3],
        c: [f64; 3],
        origin: [f64; 3],
    ) -> Result<Self, Error> {
        let h = [a, b, c];
        let (hinv, det) = invert_columns(h).ok_or(Error::BadBox)?;
        if !det.is_finite() || det.abs() < 1e-18 {
            return Err(Error::BadBox);
        }
        let bc = cross(b, c);
        let ca = cross(c, a);
        let ab = cross(a, b);
        let wa = det.abs() / norm(bc);
        let wb = det.abs() / norm(ca);
        let wc = det.abs() / norm(ab);
        if !(wa > 0.0 && wb > 0.0 && wc > 0.0) {
            return Err(Error::BadBox);
        }
        Ok(Self {
            h,
            hinv,
            origin,
            widths: [wa, wb, wc],
        })
    }

    /// Lattice vector a (first column of H).
    pub fn a(&self) -> [f64; 3] {
        self.h[0]
    }

    /// Lattice vector b.
    pub fn b(&self) -> [f64; 3] {
        self.h[1]
    }

    /// Lattice vector c.
    pub fn c(&self) -> [f64; 3] {
        self.h[2]
    }

    /// Dump-cell origin.
    pub fn origin(&self) -> [f64; 3] {
        self.origin
    }

    /// Perpendicular widths of the three faces.
    pub fn widths(&self) -> [f64; 3] {
        self.widths
    }

    /// Fractional coordinates in [0, 1).
    #[inline]
    pub fn fractional(&self, r: [f64; 3]) -> [f64; 3] {
        let d = [
            r[0] - self.origin[0],
            r[1] - self.origin[1],
            r[2] - self.origin[2],
        ];
        let mut s = mul(self.hinv, d);
        for e in &mut s {
            *e -= e.floor();
            if *e >= 1.0 {
                *e = 0.0;
            }
        }
        s
    }

    /// Cartesian from fractional.
    #[inline]
    pub fn cartesian(&self, s: [f64; 3]) -> [f64; 3] {
        let r = mul(self.h, s);
        [
            r[0] + self.origin[0],
            r[1] + self.origin[1],
            r[2] + self.origin[2],
        ]
    }

    /// Squared minimum-image distance.
    #[inline]
    pub fn dist2(&self, p: [f64; 3], q: [f64; 3]) -> f64 {
        let dp = [q[0] - p[0], q[1] - p[1], q[2] - p[2]];
        let mut ds = mul(self.hinv, dp);
        for e in &mut ds {
            *e -= e.round();
        }
        let dr = mul(self.h, ds);
        dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2]
    }
}

/// H is stored by columns. `mul(h, s)` is H s.
#[inline]
fn mul(h: [[f64; 3]; 3], v: [f64; 3]) -> [f64; 3] {
    [
        h[0][0] * v[0] + h[1][0] * v[1] + h[2][0] * v[2],
        h[0][1] * v[0] + h[1][1] * v[1] + h[2][1] * v[2],
        h[0][2] * v[0] + h[1][2] * v[1] + h[2][2] * v[2],
    ]
}

fn cross(u: [f64; 3], v: [f64; 3]) -> [f64; 3] {
    [
        u[1] * v[2] - u[2] * v[1],
        u[2] * v[0] - u[0] * v[2],
        u[0] * v[1] - u[1] * v[0],
    ]
}

fn norm(v: [f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

fn invert_columns(h: [[f64; 3]; 3]) -> Option<([[f64; 3]; 3], f64)> {
    let a = h[0];
    let b = h[1];
    let c = h[2];
    let det = a[0] * (b[1] * c[2] - b[2] * c[1])
        - a[1] * (b[0] * c[2] - b[2] * c[0])
        + a[2] * (b[0] * c[1] - b[1] * c[0]);
    if !det.is_finite() || det.abs() < 1e-18 {
        return None;
    }
    let invdet = 1.0 / det;
    // Inverse of [a b c]: rows are (b x c, c x a, a x b) / det,
    // so columns of Hinv are those divided by det... 
    // Hinv_{ij} such that Hinv * H = I.
    // Cofactor transpose / det.
    let inv = [
        [
            (b[1] * c[2] - b[2] * c[1]) * invdet,
            (a[2] * c[1] - a[1] * c[2]) * invdet,
            (a[1] * b[2] - a[2] * b[1]) * invdet,
        ],
        [
            (b[2] * c[0] - b[0] * c[2]) * invdet,
            (a[0] * c[2] - a[2] * c[0]) * invdet,
            (a[2] * b[0] - a[0] * b[2]) * invdet,
        ],
        [
            (b[0] * c[1] - b[1] * c[0]) * invdet,
            (a[1] * c[0] - a[0] * c[1]) * invdet,
            (a[0] * b[1] - a[1] * b[0]) * invdet,
        ],
    ];
    // inv is stored as inv[col][row] matching H.
    Some((inv, det))
}
