//! Galaxies as formulas: a few dozen numbers each, fitted to a photograph by tools/galaxy_fit.py.
//!
//! A galaxy here is a round bulge (a Sersic profile), an exponential disc with spiral arms and
//! dust lanes along them, and a coarse map (32 x 32) of what the photograph has that the formula
//! has not. From these numbers alone the star formula (procedural.rs) gets how many stars there
//! are at any place, so that a catalog of many galaxies costs a few hundred bytes each and no star
//! is stored: as with the Milky Way, only the stars around the camera are made.
//!
//! With its photograph at hand (the usual case), the galaxy is made from the photograph itself:
//! each pixel is put back where it lies on the disc (the disc's tilt and the direction of its long
//! axis say where), and its light is how many stars are there. The photograph is the one of the
//! data folder, so this costs nothing to store but those few numbers. The formula is for the
//! galaxies without a photograph.
//!
//! Lengths along the galaxy's own axes (the long axis of its disc as seen from the Earth, its
//! pole, the short axis), as in the voxel grids of photo_galaxy.py.

use std::sync::Arc;

use glam::{DMat3, DVec3};
use serde_json::Value;

use crate::procedural::Rng;

#[derive(Clone, Debug, PartialEq)]
pub struct ParamGalaxy {
    pub name: String,
    pub names: Vec<String>,
    /// The photograph and where it has the galaxy: its middle (pixels) and scale (kpc a pixel).
    pub photo: String,
    pub photo_middle: [f64; 2],
    pub photo_kpc_per_px: f64,
    /// The disc's long and short axes in the photograph (to the right, down).
    pub photo_long: [f64; 2],
    pub photo_short: [f64; 2],
    /// Middle (galactic, parsecs) and the rows that turn galactic axes into its own.
    pub centre: DVec3,
    pub to_box: DMat3,
    pub systems: f64,
    pub abs_mag: f64,
    pub bulge_share: f64,
    pub bulge_radius: f64,
    pub sersic: f64,
    pub bulge_flat: f64,
    pub disc_scale: f64,
    pub disc_height: f64,
    pub disc_edge: f64,
    pub arm_count: f64,
    /// m over the tangent of the pitch, signed: how the arms wind with the log of the radius.
    pub arm_wind: f64,
    pub arm_phase: f64,
    pub arm_strength: f64,
    pub arm_sharpness: i32,
    pub arm_start: f64,
    pub dust: f64,
    pub dust_offset: f64,
    pub colour_middle: [f32; 3],
    pub colour_disc: [f32; 3],
    pub colour_arms: [f32; 3],
    pub detail_size: usize,
    pub detail_half: f64,
    pub detail: Vec<f32>,
    /// Filled in on loading: what the parts add up to, for systems per cubic parsec.
    disc_total: f64,
    bulge_total: f64,
    /// The disc's pattern (arms, lanes, map) at most, for drawing stars by it.
    pattern_max: f64,
    /// The bulge's radius for a share of its stars, from 0 to 1 in 1024 steps.
    bulge_radii: Vec<f64>,
    /// The photograph, read pixel by pixel: where it is, the stars are as it has them. Or the
    /// map of the disc from above that tools/galaxy_views.py made, which is the photograph where
    /// the photograph can tell.
    pub photo_map: Option<Arc<PhotoMap>>,
    /// How thick the disc is along its radius (parsecs, parsecs), from galaxy_views.py.
    pub thickness: Vec<(f64, f64)>,
    /// The side matrix: at each radius, how the light is spread in height (galaxy_views.py).
    pub side: Option<SideMap>,
    /// The middle a little thicker, as a convex lens, and some light of its own there.
    pub lens: Option<Lens>,
    /// Pictures of the map from above as the game draws its disc (rows along the long axis), for
    /// a galaxy with no photograph (the Milky Way): from afar (half its three pictures' size)
    /// and near (four times it), whole paths; empty if it has none (faceMap "far", "near").
    pub face_far: String,
    pub face_near: String,
}

/// The light of a galaxy seen edge on, as a matrix: for each of `rows` radii out to `radius`,
/// how much of it is at each of `heights` heights from -`height` to +`height` (parsecs).
#[derive(Clone, Debug, PartialEq)]
pub struct SideMap {
    pub rows: usize,
    pub heights: usize,
    pub radius: f64,
    pub height: f64,
    /// Each row's share at each height (each row adds up to 1), and the same added up.
    share: Vec<f32>,
    cumulative: Vec<f32>,
}

impl SideMap {
    fn from_json(v: &Value) -> Option<SideMap> {
        let (rows, heights) = (v["radii"].as_u64()? as usize, v["heights"].as_u64()? as usize);
        let values: Vec<f32> = v["values"].as_array()?.iter().map(|x| x.as_f64().unwrap_or(0.0) as f32).collect();
        if values.len() != rows * heights || rows == 0 || heights == 0 {
            return None;
        }
        let (mut share, mut cumulative) = (Vec::with_capacity(values.len()), Vec::with_capacity(values.len()));
        for row in values.chunks(heights) {
            let sum: f32 = row.iter().sum::<f32>().max(1e-12);
            let mut acc = 0.0;
            for &x in row {
                share.push(x / sum);
                acc += x / sum;
                cumulative.push(acc);
            }
        }
        Some(SideMap { rows, heights, radius: v["radiusKpc"].as_f64()? * 1000.0, height: v["heightKpc"].as_f64()? * 1000.0, share, cumulative })
    }

    fn row(&self, r: f64) -> usize {
        ((r / self.radius * self.rows as f64) as usize).min(self.rows - 1)
    }

    /// The two rows around a radius (their middles) and how far it is from the first to the
    /// second: the matrix read smoothly along the radius, not in steps.
    fn rows_around(&self, r: f64) -> (usize, usize, f64) {
        let f = (r / self.radius * self.rows as f64 - 0.5).clamp(0.0, (self.rows - 1) as f64);
        let k0 = f.floor() as usize;
        (k0, (k0 + 1).min(self.rows - 1), f - k0 as f64)
    }

    /// Of the light at a radius, how much there is per parsec of height at a height (the two
    /// rows around the radius blended).
    fn per_pc(&self, r: f64, y: f64) -> f64 {
        let f = (y + self.height) / (2.0 * self.height) * self.heights as f64;
        if f < 0.0 || f >= self.heights as f64 {
            return 0.0;
        }
        let (k0, k1, t) = self.rows_around(r);
        let j = f as usize;
        let share = self.share[k0 * self.heights + j] as f64 * (1.0 - t) + self.share[k1 * self.heights + j] as f64 * t;
        share / (2.0 * self.height / self.heights as f64)
    }

    /// A height at a radius, for a draw `u` (0..1) and another for within its step. Of the two
    /// rows around the radius, one is drawn as much as the radius is near it: on average the
    /// blend of the two.
    fn height_of(&self, r: f64, u: f64, within: f64) -> f64 {
        let (k0, k1, t) = self.rows_around(r);
        // the draw for the row, from the draw for the height's low digits
        let pick = if (u * 4096.0).fract() < t { k1 } else { k0 };
        let row = &self.cumulative[pick * self.heights..(pick + 1) * self.heights];
        let j = row.partition_point(|&c| (c as f64) < u).min(self.heights - 1);
        -self.height + (j as f64 + within) * 2.0 * self.height / self.heights as f64
    }
}

/// A convex lens in the middle of a disc (galaxy_views.py "lens"; lengths in parsecs).
#[derive(Clone, Debug, PartialEq)]
pub struct Lens {
    pub radius: f64,
    pub strength: f64,
    pub extra_light: f64,
    pub extra_radius: f64,
    pub extra_flat: f64,
    /// The disc's radius, which sets how far the lens pushes.
    pub disc_radius: f64,
}

impl Lens {
    fn dome(&self, r: f64) -> f64 {
        (1.0 - (r / self.radius).powi(2)).max(0.0).sqrt()
    }

    fn reach(&self, dome: f64) -> f64 {
        (0.12 * dome + 0.02) * self.disc_radius
    }

    /// A height pushed apart: near the plane by 1 + strength dome, the far wings by no more than
    /// strength dome h (tanh: no edge anywhere).
    fn push(&self, r: f64, y: f64) -> f64 {
        let dome = self.dome(r);
        let h = self.reach(dome);
        y + self.strength * dome * h * (y / h).tanh()
    }

    /// The height a pushed height came from, and how much the push stretched it there.
    fn unpush(&self, r: f64, pushed: f64) -> (f64, f64) {
        let dome = self.dome(r);
        if dome <= 0.0 {
            return (pushed, 1.0);
        }
        let h = self.reach(dome);
        let k = self.strength * dome;
        let mut y = pushed / (1.0 + k);
        for _ in 0..8 {
            let t = (y / h).tanh();
            let f = y + k * h * t - pushed;
            let df = 1.0 + k * (1.0 - t * t);
            y -= f / df;
        }
        let t = (y / h).tanh();
        (y, 1.0 + k * (1.0 - t * t))
    }
}

/// A photograph of a galaxy as the stars of its disc: the light of each pixel, put back where it
/// lies on the disc.
#[derive(Clone, PartialEq)]
pub struct PhotoMap {
    pub width: usize,
    pub height: usize,
    /// The light of each pixel, linear, the sky's taken off; nothing beyond the disc.
    pub light: Vec<f32>,
    /// The colour of each pixel over its light, in 64ths (3 bytes a pixel, not 12).
    pub colour: Vec<[u8; 3]>,
    /// Where the disc's middle is (pixels), kpc a pixel, its long and short axes in the picture
    /// (to the right, down), the cosine of its tilt, how far out it reaches (parsecs).
    pub middle: [f64; 2],
    pub pc_per_px: f64,
    pub long: [f64; 2],
    pub short: [f64; 2],
    pub cos_tilt: f64,
    pub reach: f64,
    /// The light added up pixel by pixel, for drawing stars by it.
    cumulative: Vec<f64>,
}

impl std::fmt::Debug for PhotoMap {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "PhotoMap({}x{}, {:.1} pc a pixel, tilt cos {:.2})", self.width, self.height, self.pc_per_px, self.cos_tilt)
    }
}

impl PhotoMap {
    fn load(file: &std::path::Path, meta: &Value) -> Result<PhotoMap, String> {
        let img = image::open(file).map_err(|e| format!("{}: {e}", file.display()))?.to_rgb32f();
        let (w, h) = (img.width() as usize, img.height() as usize);
        let pair = |key: &str| -> [f64; 2] {
            let a = meta[key].as_array().cloned().unwrap_or_default();
            [a.first().and_then(Value::as_f64).unwrap_or(0.0), a.get(1).and_then(Value::as_f64).unwrap_or(0.0)]
        };
        let mut map = PhotoMap {
            width: w,
            height: h,
            light: Vec::with_capacity(w * h),
            colour: Vec::with_capacity(w * h),
            middle: pair("middlePx"),
            pc_per_px: meta["kpcPerPx"].as_f64().ok_or(format!("{}: no kpcPerPx", file.display()))? * 1000.0,
            long: pair("longAxisPx"),
            short: pair("shortAxisPx"),
            cos_tilt: meta["cosTilt"].as_f64().ok_or(format!("{}: no cosTilt", file.display()))?.max(0.05),
            reach: meta["reachKpc"].as_f64().unwrap_or(1e9) * 1000.0,
            cumulative: Vec::new(),
        };
        let linear: Vec<[f32; 3]> = img.pixels().map(|p| p.0.map(|c| c.powf(2.2))).collect();
        // the sky: the light along the picture's edges
        let mut edge: Vec<f32> = (0..w * h).filter(|i| i % w == 0 || i % w == w - 1 || i / w == 0 || i / w == h - 1).map(|i| linear[i].iter().sum::<f32>() / 3.0).collect();
        edge.sort_by(f32::total_cmp);
        let sky = edge[edge.len() / 2] * 1.5;
        for (i, c) in linear.iter().enumerate() {
            let (x, z) = map.disc_of(i);
            let l = (c.iter().sum::<f32>() / 3.0 - sky).max(0.0);
            map.light.push(if x.hypot(z) < map.reach { l } else { 0.0 });
            let mean = (c.iter().sum::<f32>() / 3.0).max(1e-6);
            map.colour.push(c.map(|v| ((v / mean).min(3.98) * 64.0) as u8));
        }
        let mut sum = 0.0;
        map.cumulative = map.light.iter().map(|&l| {
            sum += l as f64;
            sum
        }).collect();
        Ok(map)
    }

    /// Where on the disc (parsecs along its long and short axes) the middle of a pixel is.
    fn disc_of(&self, i: usize) -> (f64, f64) {
        let (dx, dy) = ((i % self.width) as f64 - self.middle[0], (i / self.width) as f64 - self.middle[1]);
        ((dx * self.long[0] + dy * self.long[1]) * self.pc_per_px, (dx * self.short[0] + dy * self.short[1]) * self.pc_per_px / self.cos_tilt)
    }

    /// The light at a place of the disc (parsecs), the pixels blended.
    fn light_at(&self, x: f64, z: f64) -> f64 {
        let col = self.middle[0] + (x * self.long[0] + z * self.cos_tilt * self.short[0]) / self.pc_per_px;
        let row = self.middle[1] + (x * self.long[1] + z * self.cos_tilt * self.short[1]) / self.pc_per_px;
        if col < 0.0 || row < 0.0 || col > (self.width - 1) as f64 || row > (self.height - 1) as f64 {
            return 0.0;
        }
        let (c0, r0) = (col.floor() as usize, row.floor() as usize);
        let (c1, r1) = ((c0 + 1).min(self.width - 1), (r0 + 1).min(self.height - 1));
        let (tc, tr) = (col - c0 as f64, row - r0 as f64);
        let at = |c: usize, r: usize| self.light[r * self.width + c] as f64;
        at(c0, r0) * (1.0 - tc) * (1.0 - tr) + at(c1, r0) * tc * (1.0 - tr) + at(c0, r1) * (1.0 - tc) * tr + at(c1, r1) * tc * tr
    }

    /// The colour at a place of the disc (parsecs), the nearest pixel's.
    fn colour_at(&self, x: f64, z: f64) -> [f32; 3] {
        let col = self.middle[0] + (x * self.long[0] + z * self.cos_tilt * self.short[0]) / self.pc_per_px;
        let row = self.middle[1] + (x * self.long[1] + z * self.cos_tilt * self.short[1]) / self.pc_per_px;
        if col < 0.0 || row < 0.0 || col > (self.width - 1) as f64 || row > (self.height - 1) as f64 {
            return [1.0; 3];
        }
        self.colour[row.round() as usize * self.width + col.round() as usize].map(|c| c as f32 / 64.0)
    }

    /// Bytes kept for it.
    pub fn bytes(&self) -> usize {
        self.light.len() * 4 + self.colour.len() * 3 + self.cumulative.len() * 8
    }

    fn total(&self) -> f64 {
        self.cumulative.last().copied().unwrap_or(0.0)
    }

    /// How much of the disc one pixel covers, seen from above (square parsecs).
    fn pixel_area(&self) -> f64 {
        self.pc_per_px * self.pc_per_px / self.cos_tilt
    }
}

fn sersic_b(n: f64) -> f64 {
    2.0 * n - 1.0 / 3.0 + 0.009876 / n
}

fn smoothstep(lo: f64, hi: f64, v: f64) -> f64 {
    let t = ((v - lo) / (hi - lo)).clamp(0.0, 1.0);
    t * t * (3.0 - 2.0 * t)
}

impl ParamGalaxy {
    /// One galaxy of tools/galaxy_fit.py's file.
    /// `data` is the folder `$data` stands for; other photographs are relative to `base`.
    pub fn from_json(v: &Value, data: &std::path::Path, base: &std::path::Path) -> Result<ParamGalaxy, String> {
        let name = v["name"].as_str().ok_or("a galaxy has no name")?.to_owned();
        let n = |obj: &Value, key: &str| obj[key].as_f64().ok_or(format!("{name}: no {key}"));
        let v3 = |a: &Value| -> Result<DVec3, String> {
            let a = a.as_array().filter(|a| a.len() == 3).ok_or(format!("{name}: not three numbers"))?;
            Ok(DVec3::new(a[0].as_f64().unwrap_or(0.0), a[1].as_f64().unwrap_or(0.0), a[2].as_f64().unwrap_or(0.0)))
        };
        let c3 = |a: &Value| -> [f32; 3] {
            let a = a.as_array().cloned().unwrap_or_default();
            [0, 1, 2].map(|k| a.get(k).and_then(Value::as_f64).unwrap_or(1.0) as f32)
        };
        let axes = v["axes"].as_array().filter(|a| a.len() == 3).ok_or(format!("{name}: no axes"))?;
        let (a0, a1, a2) = (v3(&axes[0])?, v3(&axes[1])?, v3(&axes[2])?);
        let (b, d, arms, dust, meta, detail) = (&v["bulge"], &v["disc"], &v["arms"], &v["dust"], &v["photoMetadata"], &v["detail"]);
        let count = n(arms, "count")?;
        let pitch = n(arms, "pitchDeg")?;
        let middle = meta["middlePx"].as_array().cloned().unwrap_or_default();
        let mut g = ParamGalaxy {
            names: v["names"].as_array().map(|a| a.iter().filter_map(|s| s.as_str().map(String::from)).collect()).unwrap_or_else(|| vec![name.clone()]),
            photo: v["photo"].as_str().unwrap_or("").to_owned(),
            photo_middle: [middle.first().and_then(Value::as_f64).unwrap_or(0.0), middle.get(1).and_then(Value::as_f64).unwrap_or(0.0)],
            photo_kpc_per_px: meta["kpcPerPx"].as_f64().unwrap_or(0.0),
            photo_long: [meta["longAxisPx"][0].as_f64().unwrap_or(1.0), meta["longAxisPx"][1].as_f64().unwrap_or(0.0)],
            photo_short: [meta["shortAxisPx"][0].as_f64().unwrap_or(0.0), meta["shortAxisPx"][1].as_f64().unwrap_or(1.0)],
            centre: v3(&v["centrePc"])?,
            to_box: DMat3::from_cols(a0, a1, a2).transpose(),
            systems: n(v, "systems")?,
            abs_mag: v["absMag"].as_f64().unwrap_or(-20.0),
            bulge_share: n(b, "share")?.clamp(0.0, 1.0),
            bulge_radius: n(b, "radiusKpc")?.max(0.01) * 1000.0,
            sersic: n(b, "sersic")?.clamp(0.3, 8.0),
            bulge_flat: b["flattening"].as_f64().unwrap_or(0.7),
            disc_scale: n(d, "scaleKpc")? * 1000.0,
            disc_height: d["heightKpc"].as_f64().unwrap_or(0.3) * 1000.0,
            disc_edge: n(d, "edgeKpc")? * 1000.0,
            arm_count: count,
            arm_wind: count / pitch.abs().to_radians().tan() * pitch.signum(),
            arm_phase: n(arms, "phaseDeg")?.to_radians(),
            arm_strength: n(arms, "strength")?,
            arm_sharpness: arms["sharpness"].as_i64().unwrap_or(2) as i32,
            arm_start: n(arms, "startKpc")? * 1000.0,
            dust: n(dust, "strength")?,
            dust_offset: n(dust, "offsetDeg")?.to_radians(),
            colour_middle: c3(&v["colours"]["middle"]),
            colour_disc: c3(&v["colours"]["disc"]),
            colour_arms: c3(&v["colours"]["arms"]),
            detail_size: detail["size"].as_u64().unwrap_or(0) as usize,
            detail_half: detail["halfKpc"].as_f64().unwrap_or(1.0) * 1000.0,
            detail: detail["values"].as_array().map(|a| a.iter().map(|x| (x.as_f64().unwrap_or(0.0) / 16.0).exp() as f32).collect()).unwrap_or_default(),
            name,
            disc_total: 1.0,
            bulge_total: 1.0,
            pattern_max: 1.0,
            bulge_radii: Vec::new(),
            photo_map: None,
            thickness: Vec::new(),
            side: None,
            lens: None,
            face_far: String::new(),
            face_near: String::new(),
        };
        if g.detail.len() != g.detail_size * g.detail_size {
            g.detail_size = 0;
        }
        g.prepare();
        // the photograph, if it is there and the numbers to read it are
        let file = if g.photo.starts_with("$data") { std::path::PathBuf::from(g.photo.replace("$data", &data.to_string_lossy())) } else { base.join(&g.photo) };
        // a photograph of our own (map/photos, relative to the file): kept as a whole path, as the
        // game reads it from anywhere (GalaxyLook's disc from afar)
        if !g.photo.is_empty() && !g.photo.starts_with("$data") {
            g.photo = file.to_string_lossy().into_owned();
        }
        let face = &v["faceMap"];
        for (key, to) in [("far", &mut g.face_far), ("near", &mut g.face_near)] {
            if let Some(f) = face[key].as_str() {
                *to = base.join(f).to_string_lossy().into_owned();
            }
        }
        if let (Some(f), Some(half), Some(n)) = (face["file"].as_str(), face["halfKpc"].as_f64(), face["size"].as_f64()) {
            // the map from above is read as a photograph taken from straight above it: rows along
            // the long axis, columns along the short one
            // ("rowsAlongShort": the other way round, as galaxy_xyz.py's top view is: the Milky
            // Way's face map is its MW_xyz_y.png itself)
            let mid = (n - 1.0) / 2.0;
            let (long, short) = if face["rowsAlongShort"].as_bool().unwrap_or(false) { ([1.0, 0.0], [0.0, 1.0]) } else { ([0.0, 1.0], [1.0, 0.0]) };
            let as_photo = serde_json::json!({"middlePx": [mid, mid], "kpcPerPx": 2.0 * half / n, "longAxisPx": long, "shortAxisPx": short, "cosTilt": 1.0, "reachKpc": half});
            g.photo_map = Some(Arc::new(PhotoMap::load(&base.join(f), &as_photo)?));
        } else if meta["cosTilt"].is_number() && file.exists() {
            g.photo_map = Some(Arc::new(PhotoMap::load(&file, meta)?));
        }
        g.side = SideMap::from_json(&v["sideMap"]);
        let l = &v["lens"];
        if let (Some(r), Some(sm)) = (l["radiusKpc"].as_f64(), &g.side) {
            g.lens = Some(Lens {
                radius: r * 1000.0,
                strength: l["strength"].as_f64().unwrap_or(0.0),
                extra_light: l["extraLight"].as_f64().unwrap_or(0.0).clamp(0.0, 0.9),
                extra_radius: l["extraRadiusKpc"].as_f64().unwrap_or(1.0) * 1000.0,
                extra_flat: l["extraFlat"].as_f64().unwrap_or(0.3),
                disc_radius: sm.radius,
            });
        }
        let t = &v["thickness"];
        if let (Some(r), Some(h)) = (t["radiusKpc"].as_array(), t["heightKpc"].as_array()) {
            g.thickness = r.iter().zip(h).filter_map(|(r, h)| Some((r.as_f64()? * 1000.0, h.as_f64()? * 1000.0))).collect();
        }
        Ok(g)
    }

    /// Adds the parts up once: the disc over its plane (with its pattern), the bulge over its
    /// radius; and the radii of the bulge by share, for drawing its stars.
    fn prepare(&mut self) {
        let steps = 400;
        let half = self.disc_edge;
        let cell = 2.0 * half / steps as f64;
        let (mut total, mut top) = (0.0, 0.0f64);
        for i in 0..steps {
            for k in 0..steps {
                let (x, z) = (-half + (i as f64 + 0.5) * cell, -half + (k as f64 + 0.5) * cell);
                let r = x.hypot(z);
                let p = self.pattern(x, z);
                top = top.max(p);
                total += self.disc_profile(r) * p * cell * cell;
            }
        }
        // through the disc: sech^2(y / z0) adds up to 2 z0
        self.disc_total = total * 2.0 * self.disc_height;
        self.pattern_max = top.max(1e-9);
        let reach = (8.0 * self.bulge_radius).max(1.0);
        let n = 4096;
        let mut cumulative = Vec::with_capacity(n);
        let mut sum = 0.0;
        for i in 0..n {
            let r = (i as f64 + 0.5) / n as f64 * reach;
            sum += self.bulge_profile(r) * 4.0 * std::f64::consts::PI * r * r * (reach / n as f64);
            cumulative.push(sum);
        }
        self.bulge_total = sum * self.bulge_flat;
        self.bulge_radii = (0..1024)
            .map(|k| {
                let u = (k as f64 + 0.5) / 1024.0 * sum;
                (cumulative.partition_point(|&c| c < u) as f64 + 0.5) / n as f64 * reach
            })
            .collect();
    }

    fn disc_profile(&self, r: f64) -> f64 {
        (-r / self.disc_scale).exp() * smoothstep(self.disc_edge, self.disc_edge * 0.85, r)
    }

    fn bulge_profile(&self, r: f64) -> f64 {
        (-sersic_b(self.sersic) * ((r / self.bulge_radius).powf(1.0 / self.sersic) - 1.0)).exp()
    }

    /// The arms, the dust lanes and the map, across the plane (1 on average where there is none).
    pub fn pattern(&self, x: f64, z: f64) -> f64 {
        let r = x.hypot(z).max(1.0);
        let psi = self.arm_count * z.atan2(x) - self.arm_wind * (r / 1000.0).ln() - self.arm_phase;
        let on = smoothstep(self.arm_start * 0.6, self.arm_start * 1.4, r);
        let k = self.arm_sharpness.max(1);
        // ((1 + cos)/2)^k is on average (2k choose k) / 4^k
        let mean: f64 = (1..=k).map(|j| (k + j) as f64 / j as f64).product::<f64>() / 4f64.powi(k);
        let arms = (1.0 + self.arm_strength * on * (((1.0 + psi.cos()) / 2.0).powi(k) / mean - 1.0)).max(0.0);
        let lane = 1.0 - self.dust * on * ((1.0 + (psi - self.dust_offset).cos()) / 2.0).powi(8);
        arms * lane * self.detail_at(x, z)
    }

    fn detail_at(&self, x: f64, z: f64) -> f64 {
        let n = self.detail_size;
        if n == 0 {
            return 1.0;
        }
        let f = |v: f64| (v + self.detail_half) / (2.0 * self.detail_half) * n as f64 - 0.5;
        let (fx, fz) = (f(x).clamp(0.0, n as f64 - 1.0), f(z).clamp(0.0, n as f64 - 1.0));
        let (x0, z0) = (fx.floor() as usize, fz.floor() as usize);
        let (x1, z1) = ((x0 + 1).min(n - 1), (z0 + 1).min(n - 1));
        let (tx, tz) = (fx - x0 as f64, fz - z0 as f64);
        let at = |i: usize, k: usize| self.detail[i * n + k] as f64;
        at(x0, z0) * (1.0 - tx) * (1.0 - tz) + at(x1, z0) * tx * (1.0 - tz) + at(x0, z1) * (1.0 - tx) * tz + at(x1, z1) * tx * tz
    }

    /// How thick the disc is at a radius: thicker in the middle, where the bulge is (a photograph
    /// cannot tell the two apart).
    fn height_at(&self, r: f64, reach: f64) -> f64 {
        if let (Some(first), Some(last)) = (self.thickness.first(), self.thickness.last()) {
            if r <= first.0 {
                return first.1;
            }
            if r >= last.0 {
                return last.1;
            }
            let k = self.thickness.partition_point(|t| t.0 < r);
            let (a, b) = (self.thickness[k - 1], self.thickness[k]);
            return a.1 + (b.1 - a.1) * (r - a.0) / (b.0 - a.0);
        }
        self.disc_height * (1.0 + 1.5 * (-(r / (0.12 * reach)).powi(2)).exp())
    }

    /// Systems per cubic parsec at a place (galactic, parsecs).
    pub fn density(&self, p: DVec3) -> f64 {
        let q = self.to_box * (p - self.centre);
        if let (Some(m), Some(side)) = (&self.photo_map, &self.side) {
            // the face matrix (where on the disc) times the side matrix (how high), the lens
            // pushing the middle apart, and the lens's own light
            let r = q.x.hypot(q.z);
            let mut density = 0.0;
            if let Some(l) = &self.lens {
                let s = l.extra_radius;
                let d2 = (q.x * q.x + q.z * q.z + (q.y / l.extra_flat).powi(2)) / (s * s);
                density += self.systems * l.extra_light * (-d2).exp() / (std::f64::consts::PI.powf(1.5) * s * s * s * l.extra_flat);
            }
            if r <= m.reach {
                let (y, stretch) = self.lens.as_ref().map_or((q.y, 1.0), |l| l.unpush(r, q.y));
                let disc_share = 1.0 - self.lens.as_ref().map_or(0.0, |l| l.extra_light);
                density += self.systems * disc_share * m.light_at(q.x, q.z) / m.total().max(1e-30) / m.pixel_area() * side.per_pc(r, y) / stretch;
            }
            return density;
        }
        if let Some(m) = &self.photo_map {
            let r = q.x.hypot(q.z);
            if r > m.reach {
                return 0.0;
            }
            let hz = self.height_at(r, m.reach);
            let sech = 1.0 / (q.y / hz).cosh();
            // the pixel's share of the stars, over the area it covers, through the disc
            return self.systems * m.light_at(q.x, q.z) / m.total().max(1e-30) / m.pixel_area() * sech * sech / (2.0 * hz);
        }
        let r = q.x.hypot(q.z);
        let reach = self.disc_edge.max(8.0 * self.bulge_radius);
        if r > reach || q.y.abs() > reach {
            return 0.0;
        }
        let sech = 1.0 / (q.y / self.disc_height).cosh();
        let disc = self.disc_profile(r) * self.pattern(q.x, q.z) * sech * sech / self.disc_total;
        let rb = (q.x * q.x + q.z * q.z + (q.y / self.bulge_flat).powi(2)).sqrt();
        let bulge = if self.bulge_share > 0.0 { self.bulge_profile(rb) / self.bulge_total } else { 0.0 };
        self.systems * ((1.0 - self.bulge_share) * disc + self.bulge_share * bulge)
    }

    /// Whether a cube (corner, side) could hold any of its stars.
    pub fn touches(&self, corner: DVec3, side: f64) -> bool {
        let mid = corner + DVec3::splat(side * 0.5);
        let reach = self.photo_map.as_ref().map_or(self.disc_edge.max(8.0 * self.bulge_radius) * 1.2, |m| m.reach);
        mid.distance(self.centre) < reach + side
    }

    /// How the light lies in height, for drawing the galaxy as a volume: `rows` radii from the
    /// middle out to the returned radius by `heights` heights over twice the returned height
    /// (parsecs), each value the light per parsec of height at that radius and height (each row
    /// adds up, over its heights, to 1). From the side matrix where it has one, else the disc's
    /// thickness (a sech²).
    pub fn height_profile(&self, rows: usize, heights: usize) -> (Vec<f32>, f64, f64) {
        let reach = self.reach();
        let (radius, half) = match &self.side {
            Some(side) => (side.radius, side.height),
            None => (reach, (4.0 * self.height_at(0.0, reach)).max(4.0 * self.disc_height)),
        };
        let dy = 2.0 * half / heights as f64;
        let mut out = Vec::with_capacity(rows * heights);
        for i in 0..rows {
            let r = (i as f64 + 0.5) / rows as f64 * radius;
            let row: Vec<f64> = (0..heights)
                .map(|j| {
                    let y = -half + (j as f64 + 0.5) * dy;
                    match &self.side {
                        Some(side) => side.per_pc(r, y),
                        None => {
                            let hz = self.height_at(r, reach);
                            let sech = 1.0 / (y / hz).cosh();
                            sech * sech / (2.0 * hz)
                        }
                    }
                })
                .collect();
            let sum = row.iter().sum::<f64>() * dy;
            out.extend(row.iter().map(|v| (v / sum.max(1e-30)) as f32));
        }
        (out, radius, half)
    }

    /// How far out the galaxy reaches (parsecs).
    pub fn reach(&self) -> f64 {
        self.photo_map.as_ref().map_or(self.disc_edge, |m| m.reach)
    }

    /// The galaxy as seen from straight above, for drawing it from afar as one textured disc:
    /// `n` x `n` pixels over its reach (rows along its long axis, columns along its short one),
    /// RGB, as light^(1/2.2) with the brightest 0.1% white. The light of its stars is the
    /// photograph's (or the face-on map's) where it has one, else the formula's.
    pub fn face_image(&self, n: usize) -> Vec<u8> {
        let reach = self.reach();
        let mut rgb = vec![[0.0f32; 3]; n * n];
        for i in 0..n {
            for k in 0..n {
                let (x, z) = ((i as f64 + 0.5) / n as f64 * 2.0 - 1.0, (k as f64 + 0.5) / n as f64 * 2.0 - 1.0);
                let (x, z) = (x * reach, z * reach);
                if x.hypot(z) > reach {
                    continue;
                }
                rgb[i * n + k] = match &self.photo_map {
                    Some(m) => {
                        let l = m.light_at(x, z) as f32;
                        m.colour_at(x, z).map(|c| c * l)
                    }
                    None => {
                        let l = (self.disc_profile(x.hypot(z)) * self.pattern(x, z)) as f32;
                        self.colour_disc.map(|c| c * l)
                    }
                };
            }
        }
        let mut lums: Vec<f32> = rgb.iter().map(|c| (c[0] + c[1] + c[2]) / 3.0).filter(|&l| l > 0.0).collect();
        lums.sort_by(f32::total_cmp);
        let top = lums.get((lums.len() as f64 * 0.999) as usize).copied().unwrap_or(1.0).max(1e-12);
        rgb.iter().flat_map(|c| c.map(|v| ((v / top).clamp(0.0, 1.0).powf(1.0 / 2.2) * 255.0 + 0.5) as u8)).collect()
    }

    /// The glow of the galaxy as `count` points that share its light, drawn by the formula:
    /// places along its own axes (parsecs from its middle), colours, absolute magnitude of each.
    pub fn glow(&self, count: usize, seed: u64) -> (Vec<DVec3>, Vec<[f32; 3]>, f64) {
        let mut rng = Rng::seeded(seed ^ 0x6669_7474);
        let (mut places, mut colours) = (Vec::with_capacity(count), Vec::with_capacity(count));
        let mix = |a: [f32; 3], b: [f32; 3], t: f32| [0, 1, 2].map(|k| a[k] + (b[k] - a[k]) * t);
        if let Some(m) = &self.photo_map {
            // pixel by pixel, as much as its light; anywhere within the pixel; through the disc
            let total = m.total();
            let middle = m.colour_at(0.0, 0.0);
            while places.len() < count && total > 0.0 {
                // the lens's own light: a flattened round Gaussian, the colour of the middle
                if let Some(l) = self.lens.as_ref().filter(|l| rng.next() < l.extra_light) {
                    let mut normal = || (-2.0 * rng.next().max(1e-12).ln()).sqrt() * (std::f64::consts::TAU * rng.next()).cos();
                    let s = l.extra_radius / 2f64.sqrt();
                    places.push(DVec3::new(normal() * s, normal() * s * l.extra_flat, normal() * s));
                    colours.push(middle);
                    continue;
                }
                let u = rng.next() * total;
                let i = m.cumulative.partition_point(|&c| c < u).min(m.light.len() - 1);
                let (dx, dy) = ((i % m.width) as f64 - m.middle[0] + rng.next() - 0.5, (i / m.width) as f64 - m.middle[1] + rng.next() - 0.5);
                let (x, z) = ((dx * m.long[0] + dy * m.long[1]) * m.pc_per_px, (dx * m.short[0] + dy * m.short[1]) * m.pc_per_px / m.cos_tilt);
                let r = x.hypot(z);
                let y = match &self.side {
                    // how high: from the side matrix's row at this radius, pushed by the lens
                    Some(side) => {
                        let y = side.height_of(r, rng.next(), rng.next());
                        self.lens.as_ref().map_or(y, |l| l.push(r, y))
                    }
                    None => {
                        let u = rng.next().clamp(1e-9, 1.0 - 1e-9);
                        self.height_at(r, m.reach) / 2.0 * (u / (1.0 - u)).ln()
                    }
                };
                places.push(DVec3::new(x, y, z));
                colours.push(m.colour[i].map(|c| c as f32 / 64.0));
            }
            return (places, colours, self.abs_mag + 2.5 * (count.max(1) as f64).log10());
        }
        while places.len() < count {
            if rng.next() < self.bulge_share {
                let r = self.bulge_radii[((rng.next() * 1024.0) as usize).min(1023)];
                let (cos_t, phi) = (rng.next() * 2.0 - 1.0, rng.next() * std::f64::consts::TAU);
                let s = (1.0 - cos_t * cos_t).sqrt();
                places.push(DVec3::new(r * s * phi.cos(), r * cos_t * self.bulge_flat, r * s * phi.sin()));
                colours.push(self.colour_middle);
                continue;
            }
            // the disc: the radius from r e^(-r/h) (two draws), the angle evenly; then kept as
            // the pattern says
            let r = -self.disc_scale * (rng.next().max(1e-12) * rng.next().max(1e-12)).ln();
            if r > self.disc_edge || rng.next() > smoothstep(self.disc_edge, self.disc_edge * 0.85, r) {
                continue;
            }
            let a = rng.next() * std::f64::consts::TAU;
            let (x, z) = (r * a.cos(), r * a.sin());
            let p = self.pattern(x, z);
            if rng.next() * self.pattern_max > p {
                continue;
            }
            // sech^2 through the disc: a logistic draw
            let u = rng.next().clamp(1e-9, 1.0 - 1e-9);
            let y = self.disc_height / 2.0 * (u / (1.0 - u)).ln();
            places.push(DVec3::new(x, y, z));
            let arm = ((p - 1.0) as f32).clamp(0.0, 1.0);
            let middle = (1.0 - (r / (2.0 * self.disc_scale)) as f32).clamp(0.0, 1.0);
            colours.push(mix(mix(self.colour_disc, self.colour_arms, arm), self.colour_middle, middle * 0.5));
        }
        (places, colours, self.abs_mag + 2.5 * (count.max(1) as f64).log10())
    }
}

/// The galaxies of a file of tools/galaxy_fit.py.
/// `data` is the folder that `$data` in the photographs' paths stands for.
pub fn load(path: &std::path::Path, data: &std::path::Path) -> Result<Vec<ParamGalaxy>, String> {
    let text = std::fs::read_to_string(path).map_err(|e| format!("{}: {e}", path.display()))?;
    let v: Value = serde_json::from_str(&text).map_err(|e| format!("{}: {e}", path.display()))?;
    v["galaxies"].as_array().ok_or(format!("{}: no \"galaxies\"", path.display()))?.iter().map(|g| ParamGalaxy::from_json(g, data, path.parent().unwrap_or(std::path::Path::new(".")))).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Each galaxy of map/galaxies/fitted.json, drawn by its formula and seen from the Earth on
    /// its own photograph's pixels: it must look like the photograph (the correlation of the
    /// square roots of the light, in blocks of a 128th of the picture), about as well as the
    /// fitting tool said.
    #[test]
    fn fitted_galaxies_look_like_their_photographs() {
        let file = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../map/galaxies/fitted.json");
        let data = std::path::PathBuf::from(std::env::var("HOME").unwrap_or_default()).join(".gaiasky/data");
        let Ok(galaxies) = load(&file, &data) else { return };
        let to_eq = crate::render::galactic_to_equatorial();
        let text: Value = serde_json::from_str(&std::fs::read_to_string(&file).unwrap()).unwrap();
        for (g, entry) in galaxies.iter().zip(text["galaxies"].as_array().unwrap()) {
            let photo_file = if g.photo.starts_with("$data") { std::path::PathBuf::from(g.photo.replace("$data", &data.to_string_lossy())) } else { file.parent().unwrap().join(&g.photo) };
            let Ok(photo) = image::open(&photo_file) else { continue };
            let photo = photo.to_rgb32f();
            let (w, h) = photo.dimensions();
            let blocks = 128usize;
            let bin = |col: f64, row: f64| -> Option<usize> {
                if col < 0.0 || row < 0.0 || col >= w as f64 || row >= h as f64 {
                    return None;
                }
                Some(((row / h as f64 * blocks as f64) as usize).min(blocks - 1) * blocks + ((col / w as f64 * blocks as f64) as usize).min(blocks - 1))
            };
            let mut pic = vec![0.0f64; blocks * blocks];
            for (x, y, p) in photo.enumerate_pixels() {
                pic[bin(x as f64, y as f64).unwrap()] += p.0.iter().map(|&c| (c as f64).powf(2.2)).sum::<f64>() / 3.0;
            }
            let centre = to_eq * g.centre;
            let sight = centre.normalize();
            let north = (DVec3::Y - sight * sight.y).normalize();
            let east = north.cross(sight);
            let rad_per_px = g.photo_kpc_per_px * 1000.0 / centre.length();
            let (places, _, _) = g.glow(1_000_000, 5);
            let mut img = vec![0.0f64; blocks * blocks];
            for q in &places {
                let d = (to_eq * (g.centre + g.to_box.transpose() * *q)).normalize();
                if let Some(b) = bin(g.photo_middle[0] - d.dot(east) / rad_per_px, g.photo_middle[1] - d.dot(north) / rad_per_px) {
                    img[b] += 1.0;
                }
            }
            if let Ok(dir) = std::env::var("DUMP") {
                // the stars in colour, on the photograph's pixels, 512 across
                let px: usize = std::env::var("DUMP_PX").ok().and_then(|v| v.parse().ok()).unwrap_or(512);
                let (dw, dh) = (px, px * h as usize / w as usize);
                let (places, colours, _) = g.glow(std::env::var("DUMP_STARS").ok().and_then(|v| v.parse().ok()).unwrap_or(4_000_000), 9);
                let mut acc = vec![[0.0f64; 3]; dw * dh];
                for (q, c) in places.iter().zip(&colours) {
                    let d = (to_eq * (g.centre + g.to_box.transpose() * *q)).normalize();
                    let (col, row) = (g.photo_middle[0] - d.dot(east) / rad_per_px, g.photo_middle[1] - d.dot(north) / rad_per_px);
                    if col >= 0.0 && row >= 0.0 && col < w as f64 && row < h as f64 {
                        let k = (row / h as f64 * dh as f64) as usize * dw + (col / w as f64 * dw as f64) as usize;
                        for j in 0..3 {
                            acc[k][j] += c[j] as f64;
                        }
                    }
                }
                let mut sorted: Vec<f64> = acc.iter().map(|a| a.iter().sum::<f64>() / 3.0).collect();
                sorted.sort_by(f64::total_cmp);
                let top = sorted[(sorted.len() as f64 * 0.999) as usize].max(1e-9);
                let bytes: Vec<u8> = acc.iter().flat_map(|a| a.map(|v| ((v / top).clamp(0.0, 1.0).powf(1.0 / 2.2) * 255.0) as u8)).collect();
                image::RgbImage::from_raw(dw as u32, dh as u32, bytes).unwrap().save(format!("{dir}/stars_{}.png", g.names.last().unwrap().replace(' ', "_"))).unwrap();
                // from above (along the pole) and from the side (along the short axis)
                let reach = g.photo_map.as_ref().map_or(g.disc_edge, |m| m.reach);
                // from above, tilted by 60 and 30 degrees, and from the side: looking down at the
                // disc from that elevation (screen right along the long axis)
                // and from above turned as the photograph is (its long axis where the photograph
                // has it), over the range DUMP_ABOVE_KPC if given: to set beside a face-on picture
                // made from the photograph
                let span: f64 = std::env::var("DUMP_ABOVE_KPC").ok().and_then(|v| v.parse::<f64>().ok()).map_or(reach, |k| k * 1000.0);
                for (view, elevation) in [("above", 90.0f64), ("abovephoto", 91.0), ("tilt60", 60.0), ("tilt30", 30.0), ("side", 0.0)] {
                    let e = elevation.to_radians();
                    let reach = if elevation == 91.0 { span } else { reach };
                    let mut acc = vec![[0.0f64; 3]; px * px];
                    for (q, c) in places.iter().zip(&colours) {
                        let (u, v) = if elevation == 91.0 {
                            (q.x * g.photo_long[0] + q.z * g.photo_short[0], q.x * g.photo_long[1] + q.z * g.photo_short[1])
                        } else if elevation == 90.0 {
                            (q.z, q.x)
                        } else {
                            (q.x, -(q.y * e.cos() - q.z * e.sin()))
                        };
                        let (col, row) = ((u / reach * 0.5 + 0.5) * px as f64, (v / reach * 0.5 + 0.5) * px as f64);
                        if col >= 0.0 && row >= 0.0 && col < px as f64 && row < px as f64 {
                            let k = row as usize * px + col as usize;
                            for j in 0..3 {
                                acc[k][j] += c[j] as f64;
                            }
                        }
                    }
                    let mut sorted: Vec<f64> = acc.iter().map(|a| a.iter().sum::<f64>() / 3.0).collect();
                    sorted.sort_by(f64::total_cmp);
                    let top = sorted[(sorted.len() as f64 * 0.999) as usize].max(1e-9);
                    let bytes: Vec<u8> = acc.iter().flat_map(|a| a.map(|v| ((v / top).clamp(0.0, 1.0).powf(1.0 / 2.2) * 255.0) as u8)).collect();
                    image::RgbImage::from_raw(px as u32, px as u32, bytes).unwrap().save(format!("{dir}/{view}_{}.png", g.names.last().unwrap().replace(' ', "_"))).unwrap();
                }
            }
            let top = pic.iter().cloned().fold(0.0, f64::max);
            let (a, b): (Vec<f64>, Vec<f64>) = img.iter().zip(&pic).filter(|(i, p)| **i > 0.0 || **p > 0.02 * top).map(|(i, p)| (i.sqrt(), p.sqrt())).unzip();
            let (ma, mb) = (a.iter().sum::<f64>() / a.len() as f64, b.iter().sum::<f64>() / b.len() as f64);
            let cov: f64 = a.iter().zip(&b).map(|(x, y)| (x - ma) * (y - mb)).sum();
            let like = cov / (a.iter().map(|x| (x - ma).powi(2)).sum::<f64>() * b.iter().map(|y| (y - mb).powi(2)).sum::<f64>()).sqrt();
            let said = entry["fit"]["likePhotograph"].as_f64().unwrap_or(0.0);
            let bytes = serde_json::to_string(entry).unwrap().len();
            eprintln!("{}: like its photograph {like:.3} ({}; the formula alone {said:.3}), {bytes} bytes of numbers, {:.2e} systems", g.name, if g.photo_map.is_some() { "read from the photograph" } else { "the formula" }, g.systems);
            // from the face matrix (a picture, not the photograph) the galaxy only lies as the
            // photograph has it; from the photograph itself it looks like it
            // a face-on picture made from the photograph (stretched face on, upscaled) must look
            // like it (0.9 or more: below that it is not used); one drawn freely is only reported
            let from_photo = entry["faceMap"]["fromPhotograph"].as_bool().unwrap_or(false);
            if from_photo || g.side.is_none() {
                assert!(like > 0.9, "{}: {like:.3}", g.name);
            } else if like < 0.9 {
                eprintln!("  {}: drawn freely, {like:.3}: not to be used", g.name);
            }
        }
    }

    /// From afar a galaxy is one disc with its picture from above: that picture is what its stars
    /// look like from straight above (the light of the glow points along its pole, against the
    /// picture, both in blocks), and it costs little to keep.
    #[test]
    fn the_picture_from_above_is_the_stars_from_above() {
        let galaxies = all_galaxies();
        let mut worst = (1.0f64, String::new());
        for g in &galaxies {
            let n = 128;
            let face = g.face_image(n);
            assert_eq!(face.len(), n * n * 3);
            let pic: Vec<f64> = face.chunks(3).map(|c| (c.iter().map(|&v| (v as f64 / 255.0).powf(2.2)).sum::<f64>() / 3.0)).collect();
            let reach = g.reach();
            let (places, _, _) = g.glow(1_000_000, 3);
            let mut img = vec![0.0f64; n * n];
            for q in &places {
                let (i, k) = (((q.x / reach * 0.5 + 0.5) * n as f64) as i64, ((q.z / reach * 0.5 + 0.5) * n as f64) as i64);
                if (0..n as i64).contains(&i) && (0..n as i64).contains(&k) {
                    img[i as usize * n + k as usize] += 1.0;
                }
            }
            let (a, b): (Vec<f64>, Vec<f64>) = img.iter().zip(&pic).filter(|(i, p)| **i > 0.0 || **p > 0.0).map(|(i, p)| (i.sqrt(), p.sqrt())).unzip();
            let (ma, mb) = (a.iter().sum::<f64>() / a.len() as f64, b.iter().sum::<f64>() / b.len() as f64);
            let cov: f64 = a.iter().zip(&b).map(|(x, y)| (x - ma) * (y - mb)).sum();
            let like = cov / (a.iter().map(|x| (x - ma).powi(2)).sum::<f64>() * b.iter().map(|y| (y - mb).powi(2)).sum::<f64>()).sqrt();
            let kept = g.photo_map.as_ref().map_or(0, |m| m.bytes());
            if like < worst.0 {
                worst = (like, g.name.clone());
            }
            let _ = kept;
            assert!(like > 0.9, "{}: {like:.3}", g.name);
        }
        eprintln!("{} galaxies: from above like their stars, the least {:.3} ({})", galaxies.len(), worst.0, worst.1);
    }

    /// Both files of galaxies: with a photograph, and imagined.
    fn all_galaxies() -> Vec<ParamGalaxy> {
        let dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../map/galaxies");
        let data = std::path::PathBuf::from(std::env::var("HOME").unwrap_or_default()).join(".gaiasky/data");
        ["fitted.json", "imagined.json"].iter().flat_map(|f| load(&dir.join(f), &data).unwrap_or_default()).collect()
    }

    /// From the side, the stars are as thick as the side matrix has them at each radius (the
    /// spread in height, which the lens widens in the middle), and the middle is the thickest.
    #[test]
    fn from_the_side_as_thick_as_the_side_matrix() {
        let galaxies = all_galaxies();
        let mut checked = 0;
        for g in galaxies.iter().filter(|g| g.side.is_some()) {
            let side = g.side.as_ref().unwrap();
            let (places, _, _) = g.glow(400_000, 4);
            let rms_of = |r0: f64, r1: f64| {
                let ys: Vec<f64> = places.iter().filter(|q| (r0..r1).contains(&q.x.hypot(q.z))).map(|q| q.y).collect();
                (ys.iter().map(|y| y * y).sum::<f64>() / ys.len().max(1) as f64).sqrt()
            };
            let matrix_rms = |r: f64| {
                let k = side.row(r);
                let step = 2.0 * side.height / side.heights as f64;
                (0..side.heights).map(|j| side.share[k * side.heights + j] as f64 * (-side.height + (j as f64 + 0.5) * step).powi(2)).sum::<f64>().sqrt()
            };
            let (r_out, r_in) = (side.radius * 0.7, side.radius * 0.1);
            let (outer, inner) = (rms_of(r_out * 0.9, r_out * 1.1), rms_of(0.0, r_in * 2.0));
            let want = matrix_rms(r_out);
            eprintln!("{}: thickness (rms) {:.0} pc at {:.1} kpc (the side matrix {:.0}), {:.0} pc in the middle", g.name, outer, r_out / 1000.0, want, inner);
            assert!((outer / want - 1.0).abs() < 0.35, "{}: {outer:.0} against {want:.0}", g.name);
            assert!(inner > outer * 0.9, "{}: the middle {inner:.0} pc, further out {outer:.0} pc", g.name);
            checked += 1;
        }
        eprintln!("{checked} galaxies checked from the side");
    }

    /// A galaxy as numbers is made of stars all the same: around a camera inside it, only those
    /// near are kept, a few hundred thousand of its hundreds of billions.
    #[test]
    fn a_fitted_galaxy_is_made_of_stars_around_the_camera() {
        let file = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../map/galaxies/fitted.json");
        let data = std::path::PathBuf::from(std::env::var("HOME").unwrap_or_default()).join(".gaiasky/data");
        let Ok(galaxies) = load(&file, &data) else { return };
        // the galaxies made from their matrices (those still read from their photograph directly
        // are the old way, which a coarse sum misses)
        for g in galaxies.iter().filter(|g| g.side.is_some()) {
            let model = crate::procedural::Model { local_density: 0.0, bulge_density: 0.0, fitted: vec![g.clone()], ..Default::default() };
            // the systems add up: the density summed over the galaxy
            let reach = g.photo_map.as_ref().map_or(g.disc_edge * 1.05, |m| m.reach);
            let (steps, ysteps) = (120, 60);
            let tall = g.side.as_ref().map_or(12.0 * g.disc_height.max(g.bulge_radius * g.bulge_flat), |s| s.height * 1.6);
            let (cell, ycell) = (2.0 * reach / steps as f64, 2.0 * tall / ysteps as f64);
            let mut sum = 0.0;
            for i in 0..steps {
                for k in 0..steps {
                    for j in 0..ysteps {
                        let q = DVec3::new(-reach + (i as f64 + 0.5) * cell, -3.0 * ysteps as f64 * ycell / 6.0 + (j as f64 + 0.5) * ycell, -reach + (k as f64 + 0.5) * cell);
                        sum += g.density(g.centre + g.to_box.transpose() * q) * cell * cell * ycell;
                    }
                }
            }
            assert!((sum / g.systems - 1.0).abs() < 0.1, "{}: {sum:.3e} of {:.3e} systems", g.name, g.systems);
            let camera = g.centre + g.to_box.transpose() * DVec3::new(g.disc_scale, 100.0, 0.0);
            let mut field = crate::procedural::Field::new(model);
            while field.update(camera, 10_000) {}
            let kept = field.stars().count();
            eprintln!("{}: {kept} stars kept around a camera inside it", g.name);
            assert!(kept > 1000 && kept < crate::procedural::MAX_STARS, "{}: {kept}", g.name);
        }
    }
}
