//! A galaxy made in the game from its three pictures (an image model's, godot/scripts/galaxy_builder.gd):
//! tools/galaxy_xyz.py's plain path in Rust (no Python where the game runs), and what the game
//! draws it with:
//!   - the views put on one grid: each side view scaled so that the width holding 99% of its light
//!     is the top view's along the same axis, its middle (the peak of its light) moved onto the top
//!     view's, its mid-plane (a dust lane, where there is one) onto the middle row; then
//!     `<name>_xyz_{y,x,z}.png` (size x size) and `<name>_xyz.json` (photo_volume.rs);
//!   - its disc from afar and near (galaxy_look.gd, memory galaxy-picture-sizes: a galaxy the game
//!     makes is 512 from afar, 2048 near): the top view, rows along its long axis as
//!     Universe.galaxy_face gives it, its light to the 99.9% brightest white;
//!   - its entry (`galaxy.json`): "look", the same kind as map/galaxies/fitted.json's, and "field",
//!     as map/universe.json starField.galaxies' (Universe.add_built_galaxy).
//! The views' axes as galaxy_xyz.py's: U down the top view's rows, W along its columns; x seen along
//! U (columns +W, rows -Y), z seen along W (columns +U, rows -Y).

use glam::{DMat3, DVec3};
use image::{imageops::FilterType, RgbImage};

type Picture = Vec<[f64; 3]>;

fn load(img: &image::DynamicImage, n: u32) -> Picture {
    let rgb = image::imageops::resize(&img.to_rgb8(), n, n, FilterType::Lanczos3);
    rgb.pixels().map(|p| p.0.map(|c| c as f64 / 255.0)).collect()
}

fn light(p: &Picture) -> Vec<f64> {
    p.iter().map(|c| c.iter().map(|v| v.powf(2.2)).sum::<f64>() / 3.0).collect()
}

/// Sums over rows (one per column) and over columns (one per row).
fn sums(l: &[f64], n: usize) -> (Vec<f64>, Vec<f64>) {
    let mut cols = vec![0.0; n];
    let mut rows = vec![0.0; n];
    for y in 0..n {
        for x in 0..n {
            cols[x] += l[y * n + x];
            rows[y] += l[y * n + x];
        }
    }
    (cols, rows)
}

/// A gaussian over a profile, sigma in samples, its ends reflected (scipy's gaussian_filter1d).
fn smooth(p: &[f64], sigma: f64) -> Vec<f64> {
    let r = (4.0 * sigma).ceil() as i64;
    let w: Vec<f64> = (-r..=r).map(|k| (-(k * k) as f64 / (2.0 * sigma * sigma)).exp()).collect();
    let total: f64 = w.iter().sum();
    let n = p.len() as i64;
    let at = |i: i64| -> f64 {
        let mut i = i;
        while i < 0 || i >= n {
            i = if i < 0 { -i - 1 } else { 2 * n - i - 1 };
        }
        p[i as usize]
    };
    (0..n).map(|i| (-r..=r).map(|k| w[(k + r) as usize] * at(i + k)).sum::<f64>() / total).collect()
}

/// Where the 0.5% and 99.5% points of a light profile are, and its middle: the peak of the profile
/// smoothed over 2% of the picture, to a fraction of a pixel (pixels).
fn extent(profile: &[f64]) -> (f64, f64, f64) {
    let total: f64 = profile.iter().sum::<f64>().max(1e-30);
    let p: Vec<f64> = profile.iter().map(|v| v / total).collect();
    let mut c = 0.0;
    let (mut lo, mut hi) = (None, None);
    for (i, v) in p.iter().enumerate() {
        c += v;
        if lo.is_none() && c >= 0.005 {
            lo = Some(i);
        }
        if hi.is_none() && c >= 0.995 {
            hi = Some(i);
        }
    }
    let s = smooth(&p, p.len() as f64 * 0.02);
    let k = s.iter().enumerate().max_by(|a, b| a.1.total_cmp(b.1)).map_or(0, |m| m.0);
    let mut mid = k as f64;
    if k > 0 && k + 1 < s.len() {
        let (a, b, d) = (s[k - 1], s[k], s[k + 1]);
        if a - 2.0 * b + d < 0.0 {
            mid = k as f64 + 0.5 * (a - d) / (a - 2.0 * b + d);
        }
    }
    (lo.unwrap_or(0) as f64, hi.unwrap_or(p.len() - 1) as f64, mid)
}

/// The picture scaled by `scale` about its middle, then moved by (dx, dy) pixels; read between
/// pixels, black outside.
fn place(img: &Picture, n: usize, scale: f64, dx: f64, dy: f64) -> Picture {
    let m = (n as f64 - 1.0) / 2.0;
    let read = |x: i64, y: i64| -> [f64; 3] {
        if x < 0 || y < 0 || x >= n as i64 || y >= n as i64 {
            [0.0; 3]
        } else {
            img[y as usize * n + x as usize]
        }
    };
    let mut out = vec![[0.0; 3]; n * n];
    for y in 0..n {
        for x in 0..n {
            let sx = (x as f64 - m - dx) / scale + m;
            let sy = (y as f64 - m - dy) / scale + m;
            let (x0, y0) = (sx.floor(), sy.floor());
            let (fx, fy) = (sx - x0, sy - y0);
            let (x0, y0) = (x0 as i64, y0 as i64);
            let (a, b, c, d) = (read(x0, y0), read(x0 + 1, y0), read(x0, y0 + 1), read(x0 + 1, y0 + 1));
            out[y * n + x] = [0, 1, 2].map(|k| (a[k] * (1.0 - fx) + b[k] * fx) * (1.0 - fy) + (c[k] * (1.0 - fx) + d[k] * fx) * fy);
        }
    }
    out
}

/// The row of a side view's dust lane: in each column of the disc's middle half, the darkest row
/// between two brighter ones near the brightest row; their median. None if there is none.
fn dust_lane(l: &[f64], n: usize) -> Option<f64> {
    let (cols, _) = sums(l, n);
    let top = cols.iter().cloned().fold(0.0, f64::max);
    let lit: Vec<usize> = (0..n).filter(|&c| cols[c] > 0.05 * top).collect();
    if lit.is_empty() {
        return None;
    }
    let q = |f: f64| lit[((lit.len() - 1) as f64 * f).round() as usize];
    let (lo, hi) = (q(0.2), q(0.8));
    let m = (0..n).max_by(|&a, &b| {
        let s = |r: usize| (lo..hi).map(|c| l[r * n + c]).sum::<f64>();
        s(a).total_cmp(&s(b))
    })?;
    let w = 12usize;
    if m < w || m + w >= n {
        return None;
    }
    let mut rows = Vec::new();
    for c in (lo..hi).step_by(2) {
        let p: Vec<f64> = (m - w..=m + w).map(|r| l[r * n + c]).collect();
        let k = (2..p.len() - 2).min_by(|&a, &b| p[a].total_cmp(&p[b]))?;
        let before = p[..k].iter().cloned().fold(0.0, f64::max);
        let after = p[k + 1..].iter().cloned().fold(0.0, f64::max);
        if p[k] < 0.5 * before.min(after) {
            rows.push((m - w + k) as f64);
        }
    }
    if rows.len() <= 10 {
        return None;
    }
    rows.sort_by(f64::total_cmp);
    Some(rows[rows.len() / 2])
}

fn to_rgb(p: &Picture, n: usize) -> RgbImage {
    RgbImage::from_fn(n as u32, n as u32, |x, y| image::Rgb(p[y as usize * n + x as usize].map(|v| (v.clamp(0.0, 1.0) * 255.0 + 0.5) as u8)))
}

/// The top view as the game's disc draws it (Universe.galaxy_face): `n` x `n`, rows along the
/// galaxy's long axis (the top view's columns: transposed), nothing outside the circle, its light
/// to the 99.9% brightest white.
fn disc_picture(top: &image::DynamicImage, n: u32) -> RgbImage {
    let rgb = image::imageops::resize(&top.to_rgb8(), n, n, FilterType::Lanczos3);
    let nn = n as usize;
    let mut lin = vec![[0.0f64; 3]; nn * nn];
    for (x, y, p) in rgb.enumerate_pixels() {
        let (u, v) = ((x as f64 + 0.5) / n as f64 * 2.0 - 1.0, (y as f64 + 0.5) / n as f64 * 2.0 - 1.0);
        if u.hypot(v) <= 1.0 {
            // transposed: row x, column y
            lin[x as usize * nn + y as usize] = p.0.map(|c| (c as f64 / 255.0).powf(2.2));
        }
    }
    let mut lums: Vec<f64> = lin.iter().map(|c| (c[0] + c[1] + c[2]) / 3.0).filter(|&l| l > 0.0).collect();
    lums.sort_by(f64::total_cmp);
    let top_l = lums.get((lums.len() as f64 * 0.999) as usize).copied().unwrap_or(1.0).max(1e-12);
    RgbImage::from_fn(n, n, |x, y| image::Rgb(lin[y as usize * nn + x as usize].map(|v| ((v / top_l).clamp(0.0, 1.0).powf(1.0 / 2.2) * 255.0 + 0.5) as u8)))
}

/// What a galaxy is to be made of: its name (and the file names' stem), where it is (galactic,
/// parsecs from the Sun), its axes (galactic: long, pole, short), its size (half the side of its
/// pictures, kpc), its stars and absolute magnitude.
pub struct Spec {
    pub name: String,
    pub stem: String,
    pub centre_pc: DVec3,
    pub axes: DMat3,
    pub half_kpc: f64,
    pub stars: f64,
    pub abs_mag: f64,
    pub about: String,
}

/// Makes everything in `dir` from its `y.png` (the top view), `x.png` (seen with its longest axis
/// left to right) and `z.png` (its intermediate one): the image model's pictures as
/// tools/galaxy_gpt.py asks for them. Returns what it did, a line.
pub fn make(dir: &std::path::Path, spec: &Spec, n: u32) -> Result<String, String> {
    let open = |f: &str| image::open(dir.join(f)).map_err(|e| format!("{}: {e}", dir.join(f).display()));
    let top_img = open("y.png")?;
    let nn = n as usize;
    let face = load(&top_img, n);
    let (fcols, frows) = sums(&light(&face), nn);
    let (fx, fy) = (extent(&fcols), extent(&frows));
    let m = (nn as f64 - 1.0) / 2.0;
    // x is seen along U: its columns are W, the top view's columns; z along W: its columns are U,
    // the top view's rows. The image model's "x" (its longest axis left to right) is the one whose
    // width is the top view's longer side.
    let (a, b) = (open("x.png")?, open("z.png")?);
    let wider_rows = (fy.1 - fy.0) > (fx.1 - fx.0);
    let (along_u, along_w) = if wider_rows { (b, a) } else { (a, b) };
    let mut views = vec![("y", face.clone())];
    let mut said = Vec::new();
    for (name, img, target) in [("x", along_u, fx), ("z", along_w, fy)] {
        let p = load(&img, n);
        let (c, r) = sums(&light(&p), nn);
        let (ex, ey) = (extent(&c), extent(&r));
        let scale = (target.1 - target.0) / (ex.1 - ex.0).max(1.0);
        let dx = target.2 - (m + (ex.2 - m) * scale);
        let dy = m - (m + (ey.2 - m) * scale);
        let mut placed = place(&p, nn, scale, dx, dy);
        let lane = dust_lane(&light(&placed), nn);
        if let Some(row) = lane {
            placed = place(&placed, nn, 1.0, 0.0, m - row);
        }
        said.push(format!("{name}: width {:.0} -> {:.0} px, lane {}", ex.1 - ex.0, target.1 - target.0, lane.map_or("none".to_string(), |r| format!("{r:.0}"))));
        views.push((name, placed));
    }
    let stem = format!("{}_xyz", spec.stem);
    for (v, p) in &views {
        to_rgb(p, nn).save(dir.join(format!("{stem}_{v}.png"))).map_err(|e| e.to_string())?;
    }
    let kpc_px = 2.0 * spec.half_kpc / n as f64;
    let xyz = serde_json::json!({
        "name": spec.name, "size": n, "stars": spec.stars,
        "views": {"y": format!("{stem}_y.png"), "x": format!("{stem}_x.png"), "z": format!("{stem}_z.png")},
        "frame": {"kpcPerPx": kpc_px, "middlePx": {"col": m, "row": m}, "wAxis": [1.0, 0.0], "uAxis": [0.0, 1.0], "likePhotograph": null},
        "bulgePx": {"col": fx.2, "row": fy.2},
        "about": "made in the game from an image model's three views (viewer/src/xyz_make.rs, as tools/galaxy_xyz.py): y seen from above, rows +U, columns +W; x seen along U, columns +W, rows -Y; z seen along W, columns +U, rows -Y. Light is sRGB^2.2.",
    });
    std::fs::write(dir.join(format!("{stem}.json")), serde_json::to_string_pretty(&xyz).unwrap()).map_err(|e| e.to_string())?;
    // its disc from afar and near
    disc_picture(&top_img, 512).save(dir.join("face_far.jpg")).map_err(|e| e.to_string())?;
    disc_picture(&top_img, top_img.width().max(n).min(2048)).save(dir.join("face_near.jpg")).map_err(|e| e.to_string())?;
    // its entry: the look (fitted.json's kind, its formula numbers plain: the pictures are what is
    // drawn) and the field's
    let ax = |v: DVec3| [v.x, v.y, v.z];
    let half = spec.half_kpc;
    let look = serde_json::json!({
        "name": spec.name, "names": [spec.name], "photo": "", "made": true, "about": spec.about,
        "centrePc": ax(spec.centre_pc), "axes": [ax(spec.axes.x_axis), ax(spec.axes.y_axis), ax(spec.axes.z_axis)],
        "absMag": spec.abs_mag, "systems": spec.stars,
        "bulge": {"share": 0.0, "radiusKpc": 0.05 * half, "sersic": 4.0, "flattening": 0.7},
        "disc": {"scaleKpc": 0.2 * half, "heightKpc": 0.02 * half, "edgeKpc": 0.8 * half},
        "arms": {"count": 2, "pitchDeg": 20.0, "phaseDeg": 0.0, "strength": 0.0, "sharpness": 2, "startKpc": 0.1 * half},
        "dust": {"strength": 0.0, "offsetDeg": 0.0},
        "colours": {"middle": [1, 1, 1], "disc": [1, 1, 1], "arms": [1, 1, 1]},
        "faceMap": {"file": format!("{stem}_y.png"), "halfKpc": half, "size": n, "rowsAlongShort": true, "far": "face_far.jpg", "near": "face_near.jpg"},
    });
    let field = serde_json::json!({"name": spec.name, "photo": format!("{stem}.json"), "look": spec.name, "systems": spec.stars, "absMag": spec.abs_mag});
    std::fs::write(dir.join("galaxy.json"), serde_json::to_string_pretty(&serde_json::json!({"look": look, "field": field})).unwrap()).map_err(|e| e.to_string())?;
    Ok(said.join("; "))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A made-up galaxy of three gaussian pictures: the side views, drawn twice as wide and off
    /// the middle, come out as wide as the top view and on its middle; the files are all there
    /// and readable as the game reads them (photo_volume.rs, fitted.rs).
    #[test]
    fn three_pictures_make_a_galaxy() {
        let dir = std::env::temp_dir().join(format!("xyz_make_test_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let blob = |n: u32, sx: f64, sy: f64, cx: f64, cy: f64| {
            image::DynamicImage::ImageRgb8(RgbImage::from_fn(n, n, |x, y| {
                let (u, v) = ((x as f64 - cx) / sx, (y as f64 - cy) / sy);
                let l = (-(u * u + v * v) / 2.0).exp();
                image::Rgb([(l * 250.0) as u8, (l * 220.0) as u8, (l * 180.0) as u8])
            }))
        };
        // top: 2048, wider along the rows (U) than the columns (W)
        blob(2048, 200.0, 300.0, 1024.0, 1024.0).save(dir.join("y.png")).unwrap();
        // x: the longest axis left to right, twice the size and off the middle
        blob(1024, 300.0, 60.0, 560.0, 500.0).save(dir.join("x.png")).unwrap();
        blob(1024, 200.0, 60.0, 470.0, 530.0).save(dir.join("z.png")).unwrap();
        let spec = Spec { name: "Test Galaxy".into(), stem: "test".into(), centre_pc: DVec3::new(0.0, 0.0, 1e8), axes: DMat3::IDENTITY, half_kpc: 50.0, stars: 1e11, abs_mag: -21.0, about: "a test".into() };
        make(&dir, &spec, 256).unwrap();
        let read = |f: &str| image::open(dir.join(f)).unwrap();
        let n = 256usize;
        let l = |f: &str| light(&load(&read(f), n as u32));
        let (yc, yr) = sums(&l("test_xyz_y.png"), n);
        let (xc, _) = sums(&l("test_xyz_x.png"), n);
        let (zc, _) = sums(&l("test_xyz_z.png"), n);
        let (ey_c, ey_r, ex, ez) = (extent(&yc), extent(&yr), extent(&xc), extent(&zc));
        // x shows W (the top view's columns), z shows U (its rows)
        assert!(((ex.1 - ex.0) / (ey_c.1 - ey_c.0) - 1.0).abs() < 0.08, "x {ex:?} top cols {ey_c:?}");
        assert!(((ez.1 - ez.0) / (ey_r.1 - ey_r.0) - 1.0).abs() < 0.08, "z {ez:?} top rows {ey_r:?}");
        assert!((ex.2 - ey_c.2).abs() < 2.0 && (ez.2 - ey_r.2).abs() < 2.0, "middles {ex:?} {ez:?}");
        assert_eq!(read("face_far.jpg").width(), 512);
        assert_eq!(read("face_near.jpg").width(), 2048);
        // as the game reads them
        let v = crate::photo_volume::PhotoVolume::load(&dir.join("test_xyz.json")).unwrap();
        assert_eq!(v.n, n);
        let entry: serde_json::Value = serde_json::from_str(&std::fs::read_to_string(dir.join("galaxy.json")).unwrap()).unwrap();
        let g = crate::fitted::ParamGalaxy::from_json(&entry["look"], &dir, &dir).unwrap();
        assert_eq!(g.name, "Test Galaxy");
        assert!(g.photo_map.is_some() && g.face_far.ends_with("face_far.jpg"));
        let _ = std::fs::remove_dir_all(&dir);
    }
}
