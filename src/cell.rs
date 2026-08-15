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
    /// Axis-aligned diagonal box: MIC is three independent wraps.
    ortho: bool,
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
            ortho: is_axis_aligned(h),
        })
    }

    /// True when H is diagonal. Distances then skip the two matvecs.
    pub fn is_ortho(&self) -> bool {
        self.ortho
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
        if self.ortho {
            return [
                wrap01((r[0] - self.origin[0]) / self.widths[0]),
                wrap01((r[1] - self.origin[1]) / self.widths[1]),
                wrap01((r[2] - self.origin[2]) / self.widths[2]),
            ];
        }
        let d = [
            r[0] - self.origin[0],
            r[1] - self.origin[1],
            r[2] - self.origin[2],
        ];
        let mut s = mul(self.hinv, d);
        for e in &mut s {
            *e = wrap01(*e);
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

    /// Cartesian translation by integer lattice counts (na, nb, nc).
    /// vesin and LAMMPS add this shift once per neighbour cell so the
    /// pair loop is a plain subtract, not a minimum-image wrap.
    #[inline]
    pub fn lattice_shift(&self, na: i32, nb: i32, nc: i32) -> [f64; 3] {
        if self.ortho {
            [
                f64::from(na) * self.widths[0],
                f64::from(nb) * self.widths[1],
                f64::from(nc) * self.widths[2],
            ]
        } else {
            let a = self.h[0];
            let b = self.h[1];
            let c = self.h[2];
            let fa = f64::from(na);
            let fb = f64::from(nb);
            let fc = f64::from(nc);
            [
                fa * a[0] + fb * b[0] + fc * c[0],
                fa * a[1] + fb * b[1] + fc * c[1],
                fa * a[2] + fb * b[2] + fc * c[2],
            ]
        }
    }

    /// Squared Cartesian distance after applying a lattice shift to `q`.
    #[inline]
    pub fn dist2_shifted(&self, p: [f64; 3], q: [f64; 3], shift: [f64; 3]) -> f64 {
        let dx = q[0] + shift[0] - p[0];
        let dy = q[1] + shift[1] - p[1];
        let dz = q[2] + shift[2] - p[2];
        dx * dx + dy * dy + dz * dz
    }

    /// Squared minimum-image distance.
    #[inline]
    pub fn dist2(&self, p: [f64; 3], q: [f64; 3]) -> f64 {
        if self.ortho {
            return dist2_ortho(self.widths, p, q);
        }
        dist2_general(self.h, self.hinv, p, q)
    }
}

#[inline]
fn wrap01(mut s: f64) -> f64 {
    s -= s.floor();
    if s >= 1.0 {
        0.0
    } else {
        s
    }
}

#[inline]
fn dist2_ortho(l: [f64; 3], p: [f64; 3], q: [f64; 3]) -> f64 {
    let mut r2 = 0.0;
    for k in 0..3 {
        let mut d = (q[k] - p[k]).abs();
        d -= l[k] * (d / l[k]).round();
        r2 += d * d;
    }
    r2
}

#[inline]
fn dist2_general(h: [[f64; 3]; 3], hinv: [[f64; 3]; 3], p: [f64; 3], q: [f64; 3]) -> f64 {
    let dp = [q[0] - p[0], q[1] - p[1], q[2] - p[2]];
    let mut ds = mul(hinv, dp);
    for e in &mut ds {
        *e -= e.round();
    }
    let dr = mul(h, ds);
    dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2]
}

fn is_axis_aligned(h: [[f64; 3]; 3]) -> bool {
    let scale = (norm(h[0]) + norm(h[1]) + norm(h[2])).max(1.0);
    let tol = 1e-12 * scale;
    h[0][1].abs() <= tol
        && h[0][2].abs() <= tol
        && h[1][0].abs() <= tol
        && h[1][2].abs() <= tol
        && h[2][0].abs() <= tol
        && h[2][1].abs() <= tol
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ortho_flag_and_dist2_match_general() {
        let b = Cell::ortho(10.0, 11.0, 12.0).unwrap();
        assert!(b.is_ortho());
        let p = [0.2, 1.0, 11.5];
        let q = [9.7, 10.8, 0.4];
        let fast = b.dist2(p, q);
        let slow = dist2_general(b.h, b.hinv, p, q);
        assert!((fast - slow).abs() <= 1e-12 * (1.0 + fast.abs()));
        let near = [1.0, 2.0, 3.0];
        let far = [2.0, 3.0, 4.0];
        let zero = [0.0, 0.0, 0.0];
        assert!(
            (b.dist2_shifted(near, far, zero) - b.dist2(near, far)).abs()
                <= 1e-12
        );
        let left = [0.2, 0.0, 0.0];
        let right = [9.4, 0.0, 0.0];
        let mic = b.dist2(left, right);
        let via = b.dist2_shifted(left, right, b.lattice_shift(-1, 0, 0));
        assert!((via - mic).abs() <= 1e-12 * (1.0 + mic.abs()));
        assert!((mic - 0.64).abs() <= 1e-12);
        let sheared = Cell::from_vectors(
            [10.0, 0.0, 0.0],
            [5.0, 8.66, 0.0],
            [0.0, 0.0, 10.0],
            [0.0, 0.0, 0.0],
        )
        .unwrap();
        assert!(!sheared.is_ortho());
    }
}
