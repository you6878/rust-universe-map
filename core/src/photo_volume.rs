//! A galaxy as a volume from three pictures, one seen along each of its axes (tools/galaxy_xyz.py,
//! map/galaxies/M31_xyz.json): no formula, the pixels alone.
//!
//! The pictures' own axes: U down the top view's rows, W along its columns (in the disc's plane), Y
//! up out of the disc; n x n pixels each, `pc` parsecs a pixel:
//!   y (seen along Y, from above): rows +U, columns +W;
//!   x (seen along U): columns +W, rows -Y;
//!   z (seen along W): columns +U, rows -Y.
//! In the galaxy's own frame as the game has it (fitted.json, `ParamGalaxy::to_box`: X its long
//! axis, Y its pole, Z its short axis, parsecs from its middle) a place is at W = X w.x + Z w.y,
//! U = X u.x + Z u.y: the top view keeps the photograph's orientation (it was made from it,
//! stretched face on). The frame file says where (tools/galaxy_xyz.py --photo found it: seen from the
//! Earth, the top view laid so is most like the photograph).
//! The stars of a voxel (a pixel each way) are the galaxy's stars times two shares: the top view's
//! light there over all of its light (where on the disc), and, in that column, the product of the
//! two edge views' light at that height over the same summed over the column (how high). Dust
//! across the mid-plane of an edge view hides stars, it does not take them away: its dip is filled
//! (`fill`). The stars of the voxels add up to the galaxy's stars.
//!
//! The light is the pictures' put back (the frame file's "light": 16-bit pictures of log10 of it,
//! tools/galaxy_xyz.py --profile): the pictures are made to be looked at, their light squeezed, the
//! middle as white as the arms. `load_as_drawn` keeps the pictures' own light, for comparing them
//! with the photograph (made to be looked at too).

use glam::DVec3;
use std::path::Path;

pub struct PhotoVolume {
    pub name: String,
    /// Pixels each way, and parsecs a pixel.
    pub n: usize,
    pub pc: f64,
    /// The stars in all.
    pub stars: f64,
    /// The galaxy's middle on the pictures (pixels: column, row of the top view), and the
    /// directions of W and U in its frame ((X, Z) components).
    pub middle: (f64, f64),
    pub w: (f64, f64),
    pub u: (f64, f64),
    /// How like the photograph the top view is, seen from the Earth (from the frame file).
    pub like_photograph: f64,
    /// Light (sRGB^2.2, the mean of the three) of each picture, row by row.
    top: Vec<f32>,
    along_x: Vec<f32>,
    along_z: Vec<f32>,
    /// For each column (X, Z) through the disc: the two edge views' product summed over the heights.
    norm: Vec<f32>,
    top_sum: f64,
    /// The top view's pixel colours.
    pub top_rgb: Vec<[u8; 3]>,
    /// For drawing heights: each column of the view along U added up (row by row, column-major:
    /// column k's from k * n), and each column of the view along W's highest.
    along_x_cumulative: Vec<f32>,
    along_z_top: Vec<f32>,
}

impl PartialEq for PhotoVolume {
    fn eq(&self, other: &Self) -> bool {
        self.name == other.name && self.n == other.n && self.pc == other.pc && self.middle == other.middle && self.stars == other.stars
    }
}

/// An edge view's column with its dust lane filled: a dip (lower than the light on both sides of
/// it) is raised. Not flat (to the lower of the two sides' highest, which made a slab of even
/// stars across the mid-plane): the sides' slopes go on into it, a parabola in log of the light (a
/// Gaussian across the plane) fitted to FLANK pixels each side, peaking in the middle; no higher
/// than PEAK_MOST the sides' light. Where the sides do not make such a hill the dip is filled flat.
const FLANK: usize = 6;
const PEAK_MOST: f32 = 2.0;
fn fill(column: &mut [f32]) {
    let n = column.len();
    let mut up = vec![0.0f32; n];
    let mut m = 0.0f32;
    for (k, &v) in column.iter().enumerate() {
        m = m.max(v);
        up[k] = m;
    }
    let mut flat = column.to_vec();
    m = 0.0;
    for k in (0..n).rev() {
        m = m.max(column[k]);
        flat[k] = column[k].max(up[k].min(m));
    }
    // each run of raised pixels
    let mut k = 0;
    while k < n {
        if flat[k] <= column[k] {
            k += 1;
            continue;
        }
        let a = k;
        while k < n && flat[k] > column[k] {
            k += 1;
        }
        let b = k - 1;
        let level = flat[a];
        let sides: Vec<usize> = (a.saturating_sub(FLANK)..a).chain(b + 1..(b + 1 + FLANK).min(n)).filter(|&j| column[j] > 0.0).collect();
        let hill = (sides.len() >= 5).then(|| hill_through(column, &sides, (a + b) as f64 / 2.0)).flatten();
        for j in a..=b {
            column[j] = match hill {
                Some((c0, c1, c2)) => {
                    let z = j as f64 - (a + b) as f64 / 2.0;
                    let v = (c0 + c1 * z + c2 * z * z).exp() as f32;
                    v.clamp(level, level * PEAK_MOST)
                }
                None => level,
            };
        }
    }
}

/// log of the light = c0 + c1 z + c2 z^2 (z pixels from `middle`) through the pixels `at`, least
/// squares; only a hill (c2 < 0) with its top between the first and last of them.
fn hill_through(column: &[f32], at: &[usize], middle: f64) -> Option<(f64, f64, f64)> {
    let mut a = [[0.0f64; 3]; 3];
    let mut r = [0.0f64; 3];
    for &j in at {
        let z = j as f64 - middle;
        let p = [1.0, z, z * z];
        let y = (column[j] as f64).ln();
        for i in 0..3 {
            for k in 0..3 {
                a[i][k] += p[i] * p[k];
            }
            r[i] += p[i] * y;
        }
    }
    // Cramer's rule
    let det = |m: &[[f64; 3]; 3]| {
        m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
    };
    let d = det(&a);
    if d.abs() < 1e-12 {
        return None;
    }
    let mut c = [0.0; 3];
    for (i, ci) in c.iter_mut().enumerate() {
        let mut m = a;
        for k in 0..3 {
            m[k][i] = r[k];
        }
        *ci = det(&m) / d;
    }
    let top = -c[1] / (2.0 * c[2]);
    let (lo, hi) = (at[0] as f64 - middle, at[at.len() - 1] as f64 - middle);
    (c[2] < 0.0 && top > lo && top < hi).then_some((c[0], c[1], c[2]))
}

#[cfg(test)]
mod fill_tests {
    use super::*;

    /// A disc seen edge on (a Gaussian across its plane) with a dust lane cut out of its middle:
    /// filled, its middle is again the brightest, as the disc is, and not a flat slab.
    #[test]
    fn a_dust_lane_is_filled_as_a_hill() {
        let disc: Vec<f32> = (0..41).map(|k| (-((k as f32 - 20.0) / 6.0).powi(2) / 2.0).exp()).collect();
        let mut seen = disc.clone();
        for v in &mut seen[17..=23] {
            *v *= 0.15;
        }
        fill(&mut seen);
        let flat_level = disc[16].min(disc[24]);
        println!("disc middle {:.3}, filled {:.3} (flat would be {:.3})", disc[20], seen[20], flat_level);
        assert!(seen[20] > seen[18] && seen[18] > seen[16] * 0.999, "peaks in the middle: {:?}", &seen[15..26]);
        assert!((seen[20] - disc[20]).abs() < 0.05, "as the disc: {} vs {}", seen[20], disc[20]);
        // the sides untouched
        assert_eq!(seen[10], disc[10]);
    }
}

/// A 16-bit picture of log10 of the light (0: none) back to the light.
fn true_light_of(file: &Path, n: usize, range: (f64, f64)) -> Result<Vec<f32>, String> {
    let img = image::open(file).map_err(|e| format!("{}: {e}", file.display()))?.to_luma16();
    if img.width() as usize != n || img.height() as usize != n {
        return Err(format!("{}: {}x{}, expected {n}x{n}", file.display(), img.width(), img.height()));
    }
    Ok(img.pixels().map(|p| if p.0[0] == 0 { 0.0 } else { 10f64.powf(range.0 + (p.0[0] as f64 - 1.0) / 65534.0 * (range.1 - range.0)) as f32 }).collect())
}

fn light_of(file: &Path, n: usize) -> Result<(Vec<f32>, Vec<[u8; 3]>), String> {
    let img = image::open(file).map_err(|e| format!("{}: {e}", file.display()))?.to_rgb8();
    if img.width() as usize != n || img.height() as usize != n {
        return Err(format!("{}: {}x{}, expected {n}x{n}", file.display(), img.width(), img.height()));
    }
    let lin = |c: u8| (c as f32 / 255.0).powf(2.2);
    Ok((img.pixels().map(|p| (lin(p.0[0]) + lin(p.0[1]) + lin(p.0[2])) / 3.0).collect(), img.pixels().map(|p| p.0).collect()))
}

/// A galaxy's three pictures, read only while they may be wanted: the map has a dozen such
/// galaxies, each tens of megabytes and a third of a second to read, and only those near the
/// camera are drawn star by star or as voxels (from afar a galaxy is its photograph). What is known
/// from its frame file alone (its name, how far its pictures reach) is kept always.
pub struct LazyVolume {
    pub file: std::path::PathBuf,
    pub name: String,
    /// How far its pictures reach from its middle (parsecs: along the disc, up), as
    /// `PhotoVolume::reach`.
    pub reach: (f64, f64),
    slot: std::sync::RwLock<Option<std::sync::Arc<PhotoVolume>>>,
    reading: std::sync::Mutex<()>,
}

impl PartialEq for LazyVolume {
    fn eq(&self, other: &Self) -> bool {
        self.file == other.file
    }
}

impl LazyVolume {
    /// From the frame file (tools/galaxy_xyz.py), its pictures not read yet.
    pub fn new(file: &Path) -> Result<LazyVolume, String> {
        let v: serde_json::Value = serde_json::from_str(&std::fs::read_to_string(file).map_err(|e| format!("{}: {e}", file.display()))?).map_err(|e| format!("{}: {e}", file.display()))?;
        let n = v["size"].as_u64().ok_or("no size")? as f64;
        let pc = v["frame"]["kpcPerPx"].as_f64().ok_or("no frame.kpcPerPx")? * 1000.0;
        let h = n / 2.0 * pc;
        Ok(LazyVolume {
            file: file.to_path_buf(),
            name: v["name"].as_str().unwrap_or("").to_owned(),
            reach: (h * std::f64::consts::SQRT_2, h),
            slot: Default::default(),
            reading: Default::default(),
        })
    }

    /// Its pictures if they are read.
    pub fn get(&self) -> Option<std::sync::Arc<PhotoVolume>> {
        self.slot.read().unwrap().clone()
    }

    /// Its pictures, read now if they are not (once, however many ask at the same time).
    pub fn load(&self) -> Result<std::sync::Arc<PhotoVolume>, String> {
        if let Some(v) = self.get() {
            return Ok(v);
        }
        let _one = self.reading.lock().unwrap();
        if let Some(v) = self.get() {
            return Ok(v);
        }
        let v = std::sync::Arc::new(PhotoVolume::load(&self.file)?);
        *self.slot.write().unwrap() = Some(v.clone());
        Ok(v)
    }

    /// Lets its pictures go (those still using them keep them until they are done).
    pub fn release(&self) {
        *self.slot.write().unwrap() = None;
    }
}

impl PhotoVolume {
    /// From the frame file tools/galaxy_xyz.py writes (its pictures next to it), with the light put
    /// back where the file has it.
    pub fn load(file: &Path) -> Result<PhotoVolume, String> {
        Self::load_with(file, true)
    }

    /// ... with the pictures' own light, as drawn.
    #[allow(dead_code)]
    pub fn load_as_drawn(file: &Path) -> Result<PhotoVolume, String> {
        Self::load_with(file, false)
    }

    fn load_with(file: &Path, put_back: bool) -> Result<PhotoVolume, String> {
        let v: serde_json::Value = serde_json::from_str(&std::fs::read_to_string(file).map_err(|e| format!("{}: {e}", file.display()))?).map_err(|e| format!("{}: {e}", file.display()))?;
        let n = v["size"].as_u64().ok_or("no size")? as usize;
        let frame = &v["frame"];
        let pc = frame["kpcPerPx"].as_f64().ok_or("no frame.kpcPerPx")? * 1000.0;
        let pair = |k: &str| (frame[k][0].as_f64().unwrap_or(0.0), frame[k][1].as_f64().unwrap_or(0.0));
        let dir = file.parent().unwrap_or(Path::new("."));
        let view = |k: &str| dir.join(v["views"][k].as_str().unwrap_or_default());
        let (mut top, top_rgb) = light_of(&view("y"), n)?;
        let (mut along_x, _) = light_of(&view("x"), n)?;
        let (mut along_z, _) = light_of(&view("z"), n)?;
        let light = &v["light"];
        if put_back && light.is_object() {
            let range = (light["log10Range"][0].as_f64().unwrap_or(0.0), light["log10Range"][1].as_f64().unwrap_or(1.0));
            let file = |k: &str| dir.join(light["files"][k].as_str().unwrap_or_default());
            top = true_light_of(&file("y"), n, range)?;
            along_x = true_light_of(&file("x"), n, range)?;
            along_z = true_light_of(&file("z"), n, range)?;
        }
        for (pic, _) in [(&mut along_x, 0), (&mut along_z, 1)] {
            for c in 0..n {
                let mut col: Vec<f32> = (0..n).map(|r| pic[r * n + c]).collect();
                fill(&mut col);
                for r in 0..n {
                    pic[r * n + c] = col[r];
                }
            }
        }
        // the columns' sums: for each (row X of the top view, column Z), over the heights
        let mut norm = vec![0.0f32; n * n];
        for y in 0..n {
            let (ax, az) = (&along_x[y * n..(y + 1) * n], &along_z[y * n..(y + 1) * n]);
            for i in 0..n {
                let zi = az[i];
                if zi <= 0.0 {
                    continue;
                }
                let row = &mut norm[i * n..(i + 1) * n];
                for k in 0..n {
                    row[k] += ax[k] * zi;
                }
            }
        }
        let top_sum = top.iter().map(|&l| l as f64).sum();
        let mid = (n as f64 - 1.0) / 2.0;
        let mut along_x_cumulative = vec![0.0f32; n * n];
        let mut along_z_top = vec![0.0f32; n];
        for c in 0..n {
            let mut acc = 0.0f32;
            for r in 0..n {
                acc += along_x[r * n + c];
                along_x_cumulative[c * n + r] = acc;
                along_z_top[c] = along_z_top[c].max(along_z[r * n + c]);
            }
        }
        Ok(PhotoVolume {
            name: v["name"].as_str().unwrap_or_default().to_owned(),
            n,
            pc,
            stars: v["stars"].as_f64().unwrap_or(1e12),
            middle: (frame["middlePx"]["col"].as_f64().unwrap_or(mid), frame["middlePx"]["row"].as_f64().unwrap_or(mid)),
            w: pair("wAxis"),
            u: pair("uAxis"),
            like_photograph: frame["likePhotograph"].as_f64().unwrap_or(0.0),
            top,
            along_x,
            along_z,
            norm,
            top_sum,
            top_rgb,
            along_x_cumulative,
            along_z_top,
        })
    }

    /// What `place` draws by: the top view's pixels' stars added up (only over columns that have
    /// light in the edge views). Made once for many draws.
    pub fn cumulative(&self) -> Vec<f64> {
        let mut sum = 0.0;
        self.top
            .iter()
            .zip(&self.norm)
            .map(|(&t, &c)| {
                if c > 0.0 {
                    sum += t as f64;
                }
                sum
            })
            .collect()
    }

    /// A place of one of its stars, drawn as the pictures have them (its own axes, parsecs from
    /// its middle), and the top view's pixel it is over. `table` is `cumulative()`. The height:
    /// drawn from the view along U's column there, kept as often as the view along W's row there
    /// is bright (their product, as `density` has it).
    pub fn place(&self, table: &[f64], rng: &mut crate::procedural::Rng) -> (DVec3, usize) {
        let n = self.n;
        let u = rng.next() * table.last().copied().unwrap_or(0.0);
        let pixel = table.partition_point(|&c| c < u).min(n * n - 1);
        let (i, k) = (pixel / n, pixel % n);
        let column = &self.along_x_cumulative[k * n..(k + 1) * n];
        let mut j = n / 2;
        for _ in 0..64 {
            let v = rng.next() as f32 * column[n - 1];
            j = column.partition_point(|&c| c < v).min(n - 1);
            if self.along_z_top[i] <= 0.0 || (rng.next() as f32) * self.along_z_top[i] <= self.along_z[j * n + i] {
                break;
            }
        }
        (self.place_of(i as f64 + rng.next() - 0.5, k as f64 + rng.next() - 0.5, j as f64 + rng.next() - 0.5), pixel)
    }

    /// How far out its pictures reach from its middle (parsecs): along the disc, and up.
    pub fn reach(&self) -> (f64, f64) {
        let h = self.n as f64 / 2.0 * self.pc;
        (h * std::f64::consts::SQRT_2, h)
    }

    /// The pixel nearest to a place on the pictures (pixels from their first, as a float).
    fn pixel(&self, p: f64) -> Option<usize> {
        let p = (p + 0.5).floor();
        (p >= 0.0 && p < self.n as f64).then_some(p as usize)
    }

    /// Where a place of the galaxy (its own axes, parsecs from its middle) is on the pictures:
    /// (row of the top view, its column, row of the edge views), in pixels.
    pub fn on_pictures(&self, q: DVec3) -> (f64, f64, f64) {
        let (w, u) = (q.x * self.w.0 + q.z * self.w.1, q.x * self.u.0 + q.z * self.u.1);
        (self.middle.1 + u / self.pc, self.middle.0 + w / self.pc, (self.n as f64 - 1.0) / 2.0 - q.y / self.pc)
    }

    /// ... and back: the place at a voxel (row, column of the top view, row of the edge views).
    pub fn place_of(&self, row: f64, col: f64, height_row: f64) -> DVec3 {
        let (u, w) = ((row - self.middle.1) * self.pc, (col - self.middle.0) * self.pc);
        DVec3::new(w * self.w.0 + u * self.u.0, ((self.n as f64 - 1.0) / 2.0 - height_row) * self.pc, w * self.w.1 + u * self.u.1)
    }

    /// Stars per cubic parsec at a place (its own axes, parsecs from the galaxy's middle).
    pub fn density(&self, q: DVec3) -> f64 {
        let (r, c, h) = self.on_pictures(q);
        let (Some(i), Some(k), Some(j)) = (self.pixel(r), self.pixel(c), self.pixel(h)) else { return 0.0 };
        let n = self.n;
        let col = self.norm[i * n + k];
        if col <= 0.0 {
            return 0.0;
        }
        let height = self.along_x[j * n + k] as f64 * self.along_z[j * n + i] as f64 / col as f64;
        self.stars * self.top[i * n + k] as f64 / self.top_sum * height / (self.pc * self.pc * self.pc)
    }

    /// Stars per square parsec of the face view at a place (its own axes; the height not used):
    /// what its column holds, however it is spread in height.
    pub fn column(&self, q: DVec3) -> f64 {
        let (r, c, _) = self.on_pictures(q);
        let (Some(i), Some(k)) = (self.pixel(r), self.pixel(c)) else { return 0.0 };
        if self.norm[i * self.n + k] <= 0.0 {
            return 0.0;
        }
        self.stars * self.top[i * self.n + k] as f64 / self.top_sum / (self.pc * self.pc)
    }

    /// The stars of all the voxels, added up (the galaxy's stars, but for columns with no light in
    /// the edge views).
    pub fn total(&self) -> f64 {
        let n = self.n;
        let mut sum = 0.0;
        for i in 0..n {
            for k in 0..n {
                if self.norm[i * n + k] > 0.0 {
                    sum += self.top[i * n + k] as f64;
                }
            }
        }
        self.stars * sum / self.top_sum
    }

    /// Its pictures' light (the game draws its disc from afar from the top one): the top view's,
    /// the two edge views' (dust lanes filled), the columns' sums, row by row, and the top view's
    /// light in all.
    #[allow(dead_code)]
    pub fn maps(&self) -> (&[f32], &[f32], &[f32], &[f32], f64) {
        (&self.top, &self.along_x, &self.along_z, &self.norm, self.top_sum)
    }

    /// Stars per cubic parsec in a voxel by its pixels (column and row of the top view, row of the
    /// edge views); none outside the pictures.
    pub fn density_at(&self, col: i64, row: i64, h: i64) -> f32 {
        let n = self.n as i64;
        if col < 0 || row < 0 || h < 0 || col >= n || row >= n || h >= n {
            return 0.0;
        }
        let (i, k, j, n) = (row as usize, col as usize, h as usize, self.n);
        let column = self.norm[i * n + k];
        if column <= 0.0 {
            return 0.0;
        }
        let height = self.along_x[j * n + k] as f64 * self.along_z[j * n + i] as f64 / column as f64;
        (self.stars * self.top[i * n + k] as f64 / self.top_sum * height / (self.pc * self.pc * self.pc)) as f32
    }

    /// A cube of `side`^3 voxels from `origin` (column, row, height row): stars per cubic parsec,
    /// height by height, row by row (the voxels around the ship, drawn finely: GalaxyLook).
    pub fn cube(&self, origin: (i64, i64, i64), side: usize) -> Vec<f32> {
        let s = side as i64;
        let mut out = Vec::with_capacity(side * side * side);
        for h in 0..s {
            for r in 0..s {
                for c in 0..s {
                    out.push(self.density_at(origin.0 + c, origin.1 + r, origin.2 + h));
                }
            }
        }
        out
    }

    /// A box of voxels from `origin` (column, row, height row), `size` each way: stars per cubic
    /// parsec, height by height, row by row (a slab of the cube round the ship that came into it:
    /// GalaxyLook._voxels_cube).
    pub fn cube_box(&self, origin: (i64, i64, i64), size: (usize, usize, usize)) -> Vec<f32> {
        let mut out = Vec::with_capacity(size.0 * size.1 * size.2);
        for h in 0..size.2 as i64 {
            for r in 0..size.1 as i64 {
                for c in 0..size.0 as i64 {
                    out.push(self.density_at(origin.0 + c, origin.1 + r, origin.2 + h));
                }
            }
        }
        out
    }

    /// Its voxels averaged over blocks of `f`^3 (the galaxy drawn coarsely, beyond the cube round
    /// the ship): stars per cubic parsec, block height by height, row by row, n / f each way; and
    /// each block column's colour in the top view (its pixels' light-weighted mean, linear,
    /// the brightest channel 255).
    pub fn coarse(&self, f: usize) -> (usize, Vec<f32>, Vec<[u8; 3]>) {
        let n = self.n;
        let m = n / f;
        // a column's share of the stars over its heights' sum: top / top_sum / norm
        let a: Vec<f64> = (0..n * n).map(|i| if self.norm[i] > 0.0 { self.top[i] as f64 / self.norm[i] as f64 } else { 0.0 }).collect();
        let mut out = vec![0.0f64; m * m * m];
        let mut row_sums = vec![0.0f64; m];
        for h in 0..n {
            let hb = h / f;
            let ax = &self.along_x[h * n..(h + 1) * n];
            let az = &self.along_z[h * n..(h + 1) * n];
            for i in 0..n {
                let zi = az[i] as f64;
                if zi <= 0.0 {
                    continue;
                }
                row_sums.iter_mut().for_each(|v| *v = 0.0);
                let ai = &a[i * n..(i + 1) * n];
                for k in 0..n {
                    row_sums[k / f] += ai[k] * ax[k] as f64;
                }
                let base = (hb * m + i / f) * m;
                for (kb, v) in row_sums.iter().enumerate() {
                    out[base + kb] += v * zi;
                }
            }
        }
        let scale = self.stars / self.top_sum / (self.pc * self.pc * self.pc) / (f * f * f) as f64;
        let rho = out.iter().map(|v| (v * scale) as f32).collect();
        let mut colour = vec![[0u8; 3]; m * m];
        for rb in 0..m {
            for kb in 0..m {
                let mut c = [0.0f64; 3];
                for i in rb * f..(rb + 1) * f {
                    for k in kb * f..(kb + 1) * f {
                        let w = self.top[i * n + k] as f64;
                        for ch in 0..3 {
                            c[ch] += w * (self.top_rgb[i * n + k][ch] as f64 / 255.0).powf(2.2);
                        }
                    }
                }
                let most = c[0].max(c[1]).max(c[2]);
                if most > 0.0 {
                    colour[rb * m + kb] = [0, 1, 2].map(|ch| (255.0 * c[ch] / most).round() as u8);
                }
            }
        }
        (m, rho, colour)
    }

    /// The light of each picture, for comparing.
    pub fn view(&self, axis: char) -> &[f32] {
        match axis {
            'x' => &self.along_x,
            'z' => &self.along_z,
            _ => &self.top,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::map::Map;

    fn file() -> std::path::PathBuf {
        Path::new(env!("CARGO_MANIFEST_DIR")).join("../map/galaxies/M31_xyz.json")
    }
    fn andromeda() -> Option<PhotoVolume> {
        PhotoVolume::load(&file()).ok()
    }
    fn andromeda_as_drawn() -> Option<PhotoVolume> {
        PhotoVolume::load_as_drawn(&file()).ok()
    }

    /// The correlation of the square roots of two lights (as the eye sees them, roughly), where
    /// either has light.
    fn likeness(a: &[f64], b: &[f64]) -> f64 {
        let (ta, tb) = (a.iter().cloned().fold(0.0, f64::max), b.iter().cloned().fold(0.0, f64::max));
        let pairs: Vec<(f64, f64)> = a.iter().zip(b).filter(|(x, y)| **x > 0.02 * ta || **y > 0.02 * tb).map(|(x, y)| (x.sqrt(), y.sqrt())).collect();
        let n = pairs.len() as f64;
        let (ma, mb) = (pairs.iter().map(|p| p.0).sum::<f64>() / n, pairs.iter().map(|p| p.1).sum::<f64>() / n);
        let (mut sab, mut saa, mut sbb) = (0.0, 0.0, 0.0);
        for (x, y) in &pairs {
            sab += (x - ma) * (y - mb);
            saa += (x - ma) * (x - ma);
            sbb += (y - mb) * (y - mb);
        }
        sab / (saa * sbb).sqrt()
    }

    /// Every galaxy made of three pictures loads, its voxels add up to its stars, and its coarse
    /// voxels hold them (cargo test every_xyz_galaxy -- --nocapture prints each).
    #[test]
    fn every_xyz_galaxy_loads() {
        let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("../map/galaxies");
        let mut files: Vec<_> = std::fs::read_dir(&dir).unwrap().filter_map(|e| e.ok()).map(|e| e.path()).filter(|p| p.to_string_lossy().ends_with("_xyz.json")).collect();
        files.sort();
        assert!(files.len() >= 13, "{} galaxies", files.len());
        for f in files {
            let t = std::time::Instant::now();
            let v = PhotoVolume::load(&f).unwrap_or_else(|e| panic!("{}: {e}", f.display()));
            let took = t.elapsed();
            let total = v.total();
            let (along, up) = v.reach();
            let (_, rho, _) = v.coarse(4);
            let block = (4.0 * v.pc).powi(3);
            let coarse: f64 = rho.iter().map(|&r| r as f64 * block).sum();
            eprintln!("{:24} {:6.1} pc/px, reach {:6.0} x {:5.0} pc, stars {:.3e} (of {:.3e}), coarse {:.4}, like the photograph {:.3}, loaded in {took:?}",
                v.name, v.pc, along, up, total, v.stars, coarse / total, v.like_photograph);
            assert!((total / v.stars - 1.0).abs() < 0.03, "{}: {total:e} of {:e}", v.name, v.stars); // less the columns with no light in the edge views
            assert!((coarse / total - 1.0).abs() < 1e-3, "{}", v.name);
            assert!(up > 0.0 && along > up, "{}: reach {along} x {up}", v.name);
        }
    }

    /// The coarse voxels (8^3 averaged) hold the same stars as the fine ones, block by block, and
    /// a cube round a place is the fine voxels there.
    #[test]
    fn its_coarse_voxels_hold_its_stars() {
        let Some(v) = andromeda() else { return };
        let t = std::time::Instant::now();
        let (m, rho, colour) = v.coarse(8);
        let took = t.elapsed();
        let block = (8.0 * v.pc).powi(3);
        let total: f64 = rho.iter().map(|&r| r as f64 * block).sum();
        println!("coarse {m}^3 in {took:?}: stars {total:.4e} of {:.4e} ({:.4})", v.total(), total / v.total());
        assert!((total / v.total() - 1.0).abs() < 1e-3);
        // a block against its fine voxels added up
        let (cb, rb, hb) = (m / 2 + 3, m / 2 - 2, m / 2);
        let mut fine = 0.0;
        for h in 0..8 {
            for r in 0..8 {
                for c in 0..8 {
                    fine += v.density_at((cb * 8 + c) as i64, (rb * 8 + r) as i64, (hb * 8 + h) as i64) as f64;
                }
            }
        }
        let coarse = rho[(hb * m + rb) * m + cb] as f64 * 512.0;
        println!("a block: fine {fine:.5e}, coarse {coarse:.5e}");
        assert!((fine / coarse - 1.0).abs() < 1e-3);
        assert!(colour.iter().any(|c| c[0] == 255 || c[2] == 255));
        let cube = v.cube((500, 505, 500), 26);
        assert_eq!(cube.len(), 26 * 26 * 26);
        assert_eq!(cube[(3 * 26 + 2) * 26 + 1], v.density_at(501, 507, 503));
    }

    #[test]
    fn the_voxels_add_up_to_its_stars() {
        let Some(g) = andromeda() else { return };
        let t = g.total();
        eprintln!("{}: {:e} stars in all, {:.1} pc a voxel", g.name, t, g.pc);
        assert!((t / g.stars - 1.0).abs() < 0.01, "{t:e}");
        // added up voxel by voxel, too (every 2nd voxel each way, times 8)
        let mut sum = 0.0;
        for i in (0..g.n).step_by(2) {
            for k in (0..g.n).step_by(2) {
                for j in (0..g.n).step_by(2) {
                    let q = g.place_of(i as f64 + 0.5, k as f64 + 0.5, j as f64 + 0.5);
                    sum += g.density(q) * 8.0 * g.pc.powi(3);
                }
            }
        }
        eprintln!("  voxel by voxel: {sum:e}");
        assert!((sum / g.stars - 1.0).abs() < 0.05, "{sum:e}");
    }

    /// Seen along each axis, the volume gives back each picture (the top one exactly; the edge
    /// views as the product of the two can).
    #[test]
    fn seen_along_each_axis_it_is_its_picture() {
        let Some(g) = andromeda_as_drawn() else { return };
        let n = g.n;
        let mut seen = [vec![0.0f64; n * n], vec![0.0f64; n * n], vec![0.0f64; n * n]];
        for i in (0..n).step_by(2) {
            for k in (0..n).step_by(2) {
                for j in 0..n {
                    let q = g.place_of(i as f64, k as f64, j as f64);
                    let d = g.density(q);
                    if d > 0.0 {
                        seen[0][i * n + k] += d;
                        seen[1][j * n + k] += d;
                        seen[2][j * n + i] += d;
                    }
                }
            }
        }
        for (s, axis) in seen.iter().zip(['y', 'x', 'z']) {
            let pic: Vec<f64> = g.view(axis).iter().map(|&l| l as f64).collect();
            // on the pixels the sums reached (every 2nd row and column of the top view, every 2nd
            // column of the edge views)
            let (a, b): (Vec<f64>, Vec<f64>) = s.iter().zip(&pic).enumerate().filter(|(p, _)| if axis == 'y' { (p / n) % 2 == 0 && (p % n) % 2 == 0 } else { (p % n) % 2 == 0 }).map(|(_, (x, y))| (*x, *y)).unzip();
            let like = likeness(&a, &b);
            eprintln!("  seen along {axis}: like its picture {like:.3}");
            assert!(like > if axis == 'y' { 0.98 } else { 0.75 }, "{axis}: {like}");
        }
    }

    /// The decisive one: put in space where the game has Andromeda (its middle, its axes), seen
    /// from the Earth the volume is the photograph, pixel for pixel (the user's rule: 0.9 at
    /// least), as its frame file says the top view alone is. Beside it what the game draws now:
    /// its look from afar (its glow points) and its stars (the grid's).
    #[test]
    fn from_the_earth_it_is_the_photograph() {
        let Some(g) = andromeda_as_drawn() else { return };
        let Ok(map) = Map::load(None, None) else { return };
        let Some(look) = map.galaxy_looks.iter().find(|l| l.name == "Andromeda Galaxy") else { return };
        let file = std::path::PathBuf::from(look.photo.replace("$data", &map.data.to_string_lossy()));
        let Ok(photo) = image::open(&file) else { return };
        let photo = photo.to_rgb32f();
        let (w, hgt) = (photo.dimensions().0 as usize, photo.dimensions().1 as usize);
        let lum: Vec<f64> = photo.pixels().map(|p| p.0.iter().map(|&c| (c as f64).powf(2.2)).sum::<f64>() / 3.0).collect();
        // the sky's own light off (the median of the edges, as galaxy_views.py)
        let mut edge: Vec<f64> = (0..w).flat_map(|x| [lum[x], lum[(hgt - 1) * w + x]]).chain((0..hgt).flat_map(|y| [lum[y * w], lum[y * w + w - 1]])).collect();
        edge.sort_by(f64::total_cmp);
        let sky = edge[edge.len() / 2] * 1.5;
        let to_eq = crate::render::galactic_to_equatorial();
        let centre = to_eq * look.centre;
        let sight = centre.normalize();
        let north = (DVec3::Y - sight * sight.y).normalize();
        let east = north.cross(sight);
        let rad_per_px = look.photo_kpc_per_px * 1000.0 / centre.length();
        let (cx, cy) = (look.photo_middle[0], look.photo_middle[1]);
        let block = 4usize;
        let (bw, bh) = (w / block, hgt / block);
        let mut pic = vec![0.0f64; bw * bh];
        for y in 0..bh * block {
            for x in 0..bw * block {
                pic[(y / block) * bw + x / block] += (lum[y * w + x] - sky).max(0.0);
            }
        }
        // a place of the galaxy (its own axes, pc) on the photograph's blocks
        let on_photo = |q: DVec3| -> Option<usize> {
            let d = (to_eq * (look.centre + look.to_box.transpose() * q)).normalize();
            let (col, row) = (cx - d.dot(east) / rad_per_px, cy - d.dot(north) / rad_per_px);
            (col >= 0.0 && row >= 0.0 && (col as usize) < bw * block && (row as usize) < bh * block).then(|| (row as usize / block) * bw + col as usize / block)
        };
        // the volume's voxels (every 2nd each way) on the photograph
        let mut seen = vec![0.0f64; bw * bh];
        for i in (0..g.n).step_by(2) {
            for k in (0..g.n).step_by(2) {
                for j in (0..g.n).step_by(2) {
                    let q = g.place_of(i as f64, k as f64, j as f64);
                    let d = g.density(q);
                    if d > 0.0 {
                        if let Some(b) = on_photo(q) {
                            seen[b] += d;
                        }
                    }
                }
            }
        }
        let volume = likeness(&seen, &pic);
        eprintln!("from the Earth, the volume like the photograph {volume:.3} (the top view alone {:.3}); {:.1} pc a pixel", g.like_photograph, g.pc);
        // what the game draws now: its look's glow points, and the star field's grid
        let points = |places: &[DVec3], centre: DVec3, to_box: glam::DMat3| -> Vec<f64> {
            let mut img = vec![0.0f64; bw * bh];
            for q in places {
                let d = (to_eq * (centre + to_box.transpose() * *q)).normalize();
                let (col, row) = (cx - d.dot(east) / rad_per_px, cy - d.dot(north) / rad_per_px);
                if col >= 0.0 && row >= 0.0 && (col as usize) < bw * block && (row as usize) < bh * block {
                    img[(row as usize / block) * bw + col as usize / block] += 1.0;
                }
            }
            img
        };
        let (glow, _, _) = look.glow(2_000_000, 3);
        let game_look = likeness(&points(&glow, look.centre, look.to_box), &pic);
        let grid = map.star_field.as_ref().and_then(|f| f.galaxies.iter().find(|g| g.name == "Andromeda")).map(|grid| {
            let (places, _, _) = grid.glow(2_000_000, 11);
            likeness(&points(&places, grid.centre, grid.to_box), &pic)
        });
        eprintln!("the game now: its look from afar {game_look:.3}; its stars (the grid) {}", grid.map_or("-".into(), |l| format!("{l:.3}")));
        assert!(volume >= 0.9, "{volume}");
    }

    /// Against the stars the game makes now (procedural.rs: Andromeda's grid, made from the
    /// photograph): how many in all, and how dense at the same places of the galaxy.
    #[test]
    fn against_the_games_stars() {
        let Some(g) = andromeda() else { return };
        let Ok(map) = Map::load(None, None) else { return };
        let Some(field) = map.star_field.as_ref() else { return };
        let Some(look) = map.galaxy_looks.iter().find(|l| l.name == "Andromeda Galaxy") else { return };
        let Some(grid) = field.galaxies.iter().find(|g| g.name == "Andromeda") else { return };
        let per = field.stars_per_system();
        eprintln!("stars in all: the pictures {:.3e}; the game {:.3e} systems x {per:.2} = {:.3e} stars", g.total(), grid.systems, grid.systems * per);
        // the same places: the galaxy's own axes (fitted.json's, which the pictures are laid on),
        // to galactic, to the grid's
        eprintln!("  place (kpc: X long, Y up, Z short)   pictures      game (systems x {per:.2})  per pc^3");
        for (x, y, z) in [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (3.0, 0.0, 0.0), (8.0, 0.0, 0.0), (8.0, 0.5, 0.0), (0.0, 0.0, 8.0), (15.0, 0.0, 0.0), (25.0, 0.0, 0.0)] {
            let q = DVec3::new(x, y, z) * 1000.0;
            let galactic = look.centre + look.to_box.transpose() * q;
            eprintln!("  ({x:>4}, {y:>3}, {z:>4})                    {:>10.3}   {:>10.3}", g.density(q), grid.density(galactic) * per);
        }
        // where the two galaxies' middles are (the look's and the grid's), and how their axes differ
        let apart = (look.centre - grid.centre).length();
        let tilt = (look.to_box.row(1).dot(grid.to_box.row(1))).clamp(-1.0, 1.0).acos().to_degrees();
        eprintln!("  the look's middle and the grid's {apart:.0} pc apart; their poles {tilt:.1} deg apart");
    }

    /// With the light put back the middle is as much brighter than the disc as measured (Courteau
    /// et al. 2011: about 7 magnitudes from the middle to 8 kpc out, face on), and the stars are
    /// where they should be: dense in the bulge, about as near the Sun's (0.1 a cubic parsec) some
    /// way out in the disc.
    #[test]
    fn its_middle_is_as_bright_as_measured() {
        let (Some(g), Some(drawn)) = (andromeda(), andromeda_as_drawn()) else { return };
        let r = |g: &PhotoVolume, kpc: f64| g.density(DVec3::new(kpc * 1000.0, 0.0, 0.0));
        let (put_back, as_drawn) = (r(&g, 0.0) / r(&g, 8.0), r(&drawn, 0.0) / r(&drawn, 8.0));
        eprintln!("the middle against 8 kpc out, stars per pc^3: {put_back:.0} times (as drawn {as_drawn:.1} times)");
        for kpc in [0.0, 0.3, 1.0, 3.0, 5.0, 8.0, 12.0, 20.0, 25.0] {
            eprintln!("  {kpc:>4} kpc: {:>10.4} per pc^3 (as drawn {:.4})", r(&g, kpc), r(&drawn, kpc));
        }
        assert!(put_back > 100.0, "{put_back}");
        let t = g.total();
        assert!((t / g.stars - 1.0).abs() < 0.01, "{t:e}");
    }
}
#[cfg(test)]
mod in_the_game {
    use super::*;

    /// On the way in from the Earth's side (13 degrees above its disc, the trip's way), the stars
    /// by distance once the field is all made, against what its density says: how many there
    /// are near at all (blue within 10 pc, white 100, red 1000).
    #[test]
    fn by_distance_on_the_way_in() {
        let Ok(map) = crate::map::Map::load(None, None) else { return };
        let f = map.star_field.as_ref().unwrap();
        let Some(g) = f.galaxies.iter().find(|g| g.photo.is_some()) else { return };
        let to_eq = crate::render::galactic_to_equatorial();
        let sight = (to_eq * g.centre).normalize();
        for ly in [23000.0, 14400.0, 5000.0, 1500.0] {
            let p = g.centre - to_eq.transpose() * sight * ly / 3.2616;
            let rho = f.density(p);
            let mut field = crate::procedural::Field::new(f.clone());
            let mut rounds = 0;
            while field.update(p, 400) && rounds < 300 {
                rounds += 1;
            }
            let mut bands = [0usize; 4];
            for (q, _) in field.stars() {
                let d = q.distance(p);
                bands[if d < 10.0 { 0 } else if d < 100.0 { 1 } else if d < 1000.0 { 2 } else { 3 }] += 1;
            }
            let ball = |r: f64| rho * 4.0 / 3.0 * std::f64::consts::PI * r * r * r;
            eprintln!("{ly} ly in: {rho:.4} systems/pc^3 here; within 10 pc the density says {:.0}, made {}; 10-100 pc made {}; 100-1000 pc made {}; beyond {}", ball(10.0), bands[0], bands[1], bands[2], bands[3]);
        }
    }

    /// Stars brighter than a few apparent magnitudes near its middle: as many as the model says?
    #[test]
    fn bright_stars_near_its_middle() {
        // how many stars brighter than an apparent magnitude there should be in all the sky near
        // its middle (its density there, the model's share of stars brighter than each absolute
        // magnitude, out to where none would be), against what the field gives
        let Ok(map) = crate::map::Map::load(None, None) else { return };
        let f = map.star_field.as_ref().unwrap();
        let Some(g) = f.galaxies.iter().find(|g| g.photo.is_some()) else { return };
        let to_eq = crate::render::galactic_to_equatorial();
        let sight = (to_eq * g.centre).normalize();
        eprintln!("share of stars brighter than an absolute magnitude: {}", [-4.0, -1.5, 0.0, 1.4, 4.8, 10.0].iter().map(|&m| format!("{m}: {:.3e}", f.brighter_than(m))).collect::<Vec<_>>().join(", "));
        for kpc in [0.4, 0.0003] {
            let p = g.centre - to_eq.transpose() * sight * kpc * 1000.0;
            f.keep_near(p);
            let rho = f.density(p);
            // shells round the place, each with its density averaged over directions; stars as
            // bright as the model makes them (none brighter than BRIGHTEST, -7)
            let dirs: Vec<DVec3> = (0..64).map(|i| {
                let z = 1.0 - (i as f64 + 0.5) / 32.0;
                let a = i as f64 * 2.399963;
                let r = (1.0 - z * z).sqrt();
                DVec3::new(r * a.cos(), z, r * a.sin())
            }).collect();
            let shells: Vec<(f64, f64, f64)> = {
                let mut v = Vec::new();
                let mut r = 0.001f64;
                while r < 1500.0 {
                    let dr = r * 0.02;
                    let mean = dirs.iter().map(|d| f.density(p + *d * r)).sum::<f64>() / dirs.len() as f64;
                    v.push((r, dr, mean));
                    r += dr;
                }
                v
            };
            let expect = |m: f64| {
                shells.iter().map(|&(r, dr, d)| {
                    let mag = m - 5.0 * (r / 10.0).log10();
                    let share = if mag <= -7.0 { 0.0 } else { f.brighter_than(mag) - f.brighter_than(-7.0) };
                    d * 4.0 * std::f64::consts::PI * r * r * dr * share
                }).sum::<f64>()
            };
            let mut field = crate::procedural::Field::new(f.clone());
            let mut rounds = 0;
            while field.update(p, 400) && rounds < 200 {
                rounds += 1;
            }
            let seen = |m: f64| field.stars().filter(|(q, s)| s.abs_mag as f64 + 5.0 * (q.distance(p).max(1e-6) / 10.0).log10() < m).count();
            eprintln!("{kpc} kpc from its middle, {rho:.1} systems/pc^3: brighter than -4: {:.0} (field {}), -1.5: {:.0} ({}), 0: {:.0} ({}), 3: {:.0} ({}), 6: {:.0} ({})",
                expect(-4.0), seen(-4.0), expect(-1.5), seen(-1.5), expect(0.0), seen(0.0), expect(3.0), seen(3.0), expect(6.0), seen(6.0));
            // the stars seen brightest are all there (within chance): to magnitude 3
            for m in [-1.5, 0.0, 3.0] {
                let (want, got) = (expect(m), seen(m) as f64);
                assert!(got > want * 0.8 - 3.0 * want.sqrt() && got < want * 1.2 + 3.0 * want.sqrt(), "{kpc} kpc, brighter than {m}: {got} of {want:.0}");
            }
        }
    }

    /// In the game's star field (procedural.rs: the map's Andromeda spread by its three pictures):
    /// flying in from the Earth's side toward its middle, the stars around the ship, as the field
    /// keeps them (all levels of cubes, within reach of the camera), and how long they take.
    #[test]
    fn flying_in_the_stars_are_there() {
        let Ok(map) = crate::map::Map::load(None, None) else { return };
        let f = map.star_field.as_ref().unwrap();
        let Some(g) = f.galaxies.iter().find(|g| g.photo.is_some()) else { return };
        let to_eq = crate::render::galactic_to_equatorial();
        let sight = (to_eq * g.centre).normalize();
        let mut counts = Vec::new();
        for kpc in [12.0, 2.0, 0.4, 0.05] {
            let p = g.centre - to_eq.transpose() * sight * kpc * 1000.0;
            let mut field = crate::procedural::Field::new(f.clone());
            let t0 = std::time::Instant::now();
            let mut rounds = 0;
            while field.update(p, 400) && rounds < 200 {
                rounds += 1;
            }
            let n = field.stars().count();
            eprintln!("{kpc:>5} kpc from the middle (the Earth's side): {:.4} systems per pc^3, {n} stars around, in {rounds} rounds, {} ms", g.density(p), t0.elapsed().as_millis());
            counts.push(n);
        }
        assert!(counts[1] > 100_000 && counts[2] > counts[1], "{counts:?}");
    }
}
