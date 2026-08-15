//! Linked-cell k-nearest search (Allen and Tildesley).

use crate::cell::Cell;
use crate::Error;

/// Bounded max-heap of (dist2, index). k is 4 in ice; the small
/// case stays on the stack so the pair loop does not allocate.
struct KHeap {
    d2: [f64; 16],
    idx: [usize; 16],
    extra_d2: Vec<f64>,
    extra_idx: Vec<usize>,
    n: usize,
    k: usize,
}

impl KHeap {
    fn new(k: usize) -> Self {
        let mut extra_d2 = Vec::new();
        let mut extra_idx = Vec::new();
        if k > 16 {
            extra_d2.resize(k, 0.0);
            extra_idx.resize(k, 0);
        }
        Self {
            d2: [0.0; 16],
            idx: [0; 16],
            extra_d2,
            extra_idx,
            n: 0,
            k,
        }
    }

    fn d2_at(&self, t: usize) -> f64 {
        if self.k <= 16 {
            self.d2[t]
        } else {
            self.extra_d2[t]
        }
    }

    fn set(&mut self, t: usize, d2: f64, j: usize) {
        if self.k <= 16 {
            self.d2[t] = d2;
            self.idx[t] = j;
        } else {
            self.extra_d2[t] = d2;
            self.extra_idx[t] = j;
        }
    }

    fn push(&mut self, d2: f64, j: usize) {
        if self.n < self.k {
            self.set(self.n, d2, j);
            self.n += 1;
            return;
        }
        let mut worst = 0;
        for t in 1..self.n {
            if self.d2_at(t) > self.d2_at(worst) {
                worst = t;
            }
        }
        if d2 < self.d2_at(worst) {
            self.set(worst, d2, j);
        }
    }

    fn full(&self) -> bool {
        self.n >= self.k
    }

    fn worst(&self) -> f64 {
        let mut w = self.d2_at(0);
        for t in 1..self.n {
            let v = self.d2_at(t);
            if v > w {
                w = v;
            }
        }
        w
    }

    fn finish(self) -> Vec<(f64, usize)> {
        let mut pairs = Vec::with_capacity(self.n);
        for t in 0..self.n {
            let j = if self.k <= 16 {
                self.idx[t]
            } else {
                self.extra_idx[t]
            };
            pairs.push((self.d2_at(t), j));
        }
        pairs.sort_by(|a, b| a.0.total_cmp(&b.0));
        pairs
    }
}

/// One source's k nearest neighbours, nearest first.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Neighbors {
    /// Candidate indices into the input point list.
    pub indices: Vec<usize>,
    /// Squared minimum-image distances, parallel to [`Self::indices`].
    pub dist2: Vec<f64>,
}

/// k-nearest neighbours of every point (or of the masked subset).
///
/// `mask[i] == false` drops point `i` from both sources and candidates.
/// `cell_hint` is the target cell edge; `None` uses 3.0 in the same units
/// as the box. Each row has `min(k, n_active - 1)` entries.
pub fn knearest(
    xyz: &[[f64; 3]],
    simbox: &Cell,
    k: usize,
    mask: Option<&[bool]>,
    cell_hint: Option<f64>,
) -> Result<Vec<Neighbors>, Error> {
    if k == 0 {
        return Err(Error::ZeroK);
    }
    if xyz.is_empty() {
        return Err(Error::Empty);
    }
    let n = xyz.len();
    let active: Vec<usize> = (0..n)
        .filter(|&i| mask.map(|m| m.get(i).copied().unwrap_or(false)).unwrap_or(true))
        .collect();
    if active.is_empty() {
        return Ok(vec![Neighbors::default(); n]);
    }

    let mut edge = cell_hint.unwrap_or(3.0);
    if !(edge > 0.0) {
        edge = 3.0;
    }
    let w = simbox.widths();
    edge = edge.min(w[0]).min(w[1]).min(w[2]);

    let nx = (w[0] / edge).floor().max(1.0) as i32;
    let ny = (w[1] / edge).floor().max(1.0) as i32;
    let nz = (w[2] / edge).floor().max(1.0) as i32;
    let ncell = (nx * ny * nz) as usize;
    let invx = f64::from(nx);
    let invy = f64::from(ny);
    let invz = f64::from(nz);
    let cell_min = (w[0] / f64::from(nx))
        .min(w[1] / f64::from(ny))
        .min(w[2] / f64::from(nz));

    let mut head = vec![-1isize; ncell];
    let mut next = vec![-1isize; n];

    let cell_of = |i: usize| -> (i32, i32, i32) {
        let s = simbox.fractional(xyz[i]);
        let ix = ((s[0] * invx) as i32).clamp(0, nx - 1);
        let iy = ((s[1] * invy) as i32).clamp(0, ny - 1);
        let iz = ((s[2] * invz) as i32).clamp(0, nz - 1);
        (ix, iy, iz)
    };
    let cell_index = |ix: i32, iy: i32, iz: i32| -> usize {
        let cx = ix.rem_euclid(nx);
        let cy = iy.rem_euclid(ny);
        let cz = iz.rem_euclid(nz);
        ((cz * ny + cy) * nx + cx) as usize
    };

    for &i in &active {
        let (ix, iy, iz) = cell_of(i);
        let c = cell_index(ix, iy, iz);
        next[i] = head[c];
        head[c] = i as isize;
    }

    let max_reach = nx.max(ny).max(nz) / 2 + 1;
    let mut seen = vec![0u32; ncell];
    let mut stamp: u32 = 1;
    let mut out = vec![Neighbors::default(); n];

    for &i in &active {
        let mut heap = KHeap::new(k);
        let (ix, iy, iz) = cell_of(i);
        stamp = stamp.wrapping_add(1);
        if stamp == 0 {
            seen.fill(0);
            stamp = 1;
        }
        let mut reach = 1i32;
        while reach <= max_reach {
            for dx in -reach..=reach {
                for dy in -reach..=reach {
                    for dz in -reach..=reach {
                        let shell = reach == 1
                            || dx.abs() == reach
                            || dy.abs() == reach
                            || dz.abs() == reach;
                        if !shell && reach > 1 {
                            continue;
                        }
                        let c = cell_index(ix + dx, iy + dy, iz + dz);
                        if seen[c] == stamp {
                            continue;
                        }
                        seen[c] = stamp;
                        let mut j = head[c];
                        while j >= 0 {
                            let ju = j as usize;
                            if ju != i {
                                heap.push(simbox.dist2(xyz[i], xyz[ju]), ju);
                            }
                            j = next[ju];
                        }
                    }
                }
            }
            if heap.full() {
                let bound = f64::from(reach) * cell_min;
                if heap.worst() <= bound * bound {
                    break;
                }
            }
            reach += 1;
        }
        let pairs = heap.finish();
        let row = &mut out[i];
        row.dist2 = pairs.iter().map(|p| p.0).collect();
        row.indices = pairs.iter().map(|p| p.1).collect();
    }
    Ok(out)
}

/// Brute-force k-nearest. Tests and small systems only.
pub fn knearest_brute(
    xyz: &[[f64; 3]],
    simbox: &Cell,
    k: usize,
    mask: Option<&[bool]>,
) -> Result<Vec<Neighbors>, Error> {
    if k == 0 {
        return Err(Error::ZeroK);
    }
    if xyz.is_empty() {
        return Err(Error::Empty);
    }
    let n = xyz.len();
    let active: Vec<usize> = (0..n)
        .filter(|&i| mask.map(|m| m.get(i).copied().unwrap_or(false)).unwrap_or(true))
        .collect();
    let mut out = vec![Neighbors::default(); n];
    for &i in &active {
        let mut pairs: Vec<(f64, usize)> = active
            .iter()
            .copied()
            .filter(|&j| j != i)
            .map(|j| (simbox.dist2(xyz[i], xyz[j]), j))
            .collect();
        pairs.sort_by(|a, b| a.0.total_cmp(&b.0));
        pairs.truncate(k);
        out[i].dist2 = pairs.iter().map(|p| p.0).collect();
        out[i].indices = pairs.iter().map(|p| p.1).collect();
    }
    Ok(out)
}
