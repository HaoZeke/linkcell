//! Orthorhombic periodic box and minimum-image distance.

/// Orthorhombic cell, periodic in all three directions.
///
/// Coordinates are interpreted in the same units as the lengths. The origin
/// `(xlo, ylo, zlo)` is only used to bin points; distances use the lengths.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct OrthoBox {
    pub lx: f64,
    pub ly: f64,
    pub lz: f64,
    pub xlo: f64,
    pub ylo: f64,
    pub zlo: f64,
}

impl OrthoBox {
    pub fn new(lx: f64, ly: f64, lz: f64) -> Result<Self, super::Error> {
        Self::with_origin(lx, ly, lz, 0.0, 0.0, 0.0)
    }

    pub fn with_origin(
        lx: f64,
        ly: f64,
        lz: f64,
        xlo: f64,
        ylo: f64,
        zlo: f64,
    ) -> Result<Self, super::Error> {
        if !(lx > 0.0 && ly > 0.0 && lz > 0.0) {
            return Err(super::Error::BadBox);
        }
        Ok(Self {
            lx,
            ly,
            lz,
            xlo,
            ylo,
            zlo,
        })
    }

    #[inline]
    pub fn wrap(&self, x: f64, y: f64, z: f64) -> [f64; 3] {
        [
            wrap_one(x - self.xlo, self.lx),
            wrap_one(y - self.ylo, self.ly),
            wrap_one(z - self.zlo, self.lz),
        ]
    }

    /// Squared minimum-image distance, same convention as a folded dump box.
    #[inline]
    pub fn dist2(&self, a: [f64; 3], b: [f64; 3]) -> f64 {
        let dx = mic_one((a[0] - b[0]).abs(), self.lx);
        let dy = mic_one((a[1] - b[1]).abs(), self.ly);
        let dz = mic_one((a[2] - b[2]).abs(), self.lz);
        dx * dx + dy * dy + dz * dz
    }
}

#[inline]
fn wrap_one(x: f64, length: f64) -> f64 {
    let t = x / length;
    (t - t.floor()) * length
}

#[inline]
fn mic_one(d: f64, length: f64) -> f64 {
    d - length * (d / length).round()
}
