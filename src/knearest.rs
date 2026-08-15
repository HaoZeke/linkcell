//! Linked-cell k-nearest search (Allen and Tildesley).

use crate::cell::Cell;
use crate::Error;

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
    let cell_index = |mut ix: i32, mut iy: i32, mut iz: i32| -> usize {
        ix = ix.rem_euclid(nx);
        iy = iy.rem_euclid(ny);
        iz = iz.rem_euclid(nz);
        ((iz * ny + iy) * nx + ix) as usize
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
        // Max-heap of size k on IEEE bits of dist2 (distances are finite).
        let mut heap: std::collections::BinaryHeap<(u64, usize)> =
            std::collections::BinaryHeap::new();
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
                                let d2 = simbox.dist2(xyz[i], xyz[ju]);
                                let key = d2.to_bits();
                                if heap.len() < k {
                                    heap.push((key, ju));
                                } else if let Some(&(top, _)) = heap.peek() {
                                    if key < top {
                                        heap.pop();
                                        heap.push((key, ju));
                                    }
                                }
                            }
                            j = next[ju];
                        }
                    }
                }
            }
            if heap.len() >= k {
                if let Some(&(top, _)) = heap.peek() {
                    let bound = f64::from(reach) * cell_min;
                    if f64::from_bits(top) <= bound * bound {
                        break;
                    }
                }
            }
            reach += 1;
        }
        let mut pairs: Vec<(f64, usize)> = heap
            .into_iter()
            .map(|(bits, j)| (f64::from_bits(bits), j))
            .collect();
        pairs.sort_by(|a, b| a.0.total_cmp(&b.0));
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
