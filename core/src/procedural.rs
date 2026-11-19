//! Stars of the whole Milky Way, made up on demand.
//!
//! The catalog only has the stars near the Sun. Everywhere else the stars come from a model of
//! the galaxy: how many stars there are at each place, and how bright stars are. Space is cut
//! into cubes, and the stars of a cube are computed from its coordinates alone, so the same
//! stars are found again whenever the camera returns. Only the cubes around the camera exist
//! at any time.
//!
//! The cubes come in sizes that double. A cube of a given size holds the stars that can be seen
//! from about two cube sides away and no further: small cubes hold the many faint stars, large
//! cubes the few bright ones. This way the number of stars in memory stays small wherever the
//! camera is, although the model has tens of billions of them.
//!
//! Coordinates here are galactic, in parsecs, with the Sun at the origin, +Z towards the galactic
//! centre and +Y to the galactic north pole.

use std::collections::HashMap;
use std::sync::Arc;

use glam::{DMat3, DVec3};

/// Faintest apparent magnitude that the star shader draws.
const LIMIT_MAG: f64 = 8.8;
/// A star that would be this bright from the Sun is left to the catalog, so the sky of the Earth
/// stays the real one.
const CATALOG_MAG: f64 = 9.5;
/// Side of the smallest cubes, and how many sizes there are.
const SMALLEST_PC: f64 = 2.0;
const LEVELS: usize = 15;
/// A cube is alive while the camera is within this many sides of it.
const REACH: f64 = 2.0;
/// Absolute magnitudes the model covers.
const BRIGHTEST: f64 = -7.0;
/// Stars drawn in one cube at most. Where there are more, the brightest of the cube are the ones
/// drawn; the rest is in the glow of the volumes. All the cubes alive together stay within what
/// the renderer can draw (`MAX_STARS`).
pub const MAX_PER_CELL: usize = 1000;
/// Cubes alive at most: those within reach of the camera, at every size.
pub const MAX_STARS: usize = MAX_PER_CELL * LEVELS * 125;
/// Cubes computed per frame at most, so that a fast trip does not stall the picture.
const CELLS_PER_UPDATE: usize = 96;
/// A cube with more stars than MAX_PER_CELL (the middle of a galaxy) cannot keep them all: from
/// afar it keeps its brightest MAX_PER_CELL, spread over it. Near the camera that would leave out
/// most of the stars seen brightest there (a cube of the bright magnitudes is hundreds of parsecs
/// wide, and its stars near the camera are thinned as much as its far ones). So a cube that the
/// camera is near is cut into parts of about MAX_PER_CELL stars each, and the parts within
/// NEAR_PARTS of their sides of the camera are made whole: the stars seen brightest kept first.
const NEAR_PARTS: f64 = 3.0;
const MOST_PARTS: usize = 64;
/// How far ahead the cubes are made (seconds of the camera's way), and at most how many of a
/// level's reaches long (Field::update_moving).
const LOOK_AHEAD_S: f64 = 1.5;
const AHEAD_REACHES: f64 = 3.0;
/// A star's index in a part: its part (1 and on; 0 a whole cube) above these bits.
const PART_SHIFT: u32 = 10;

/// What a star is.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Kind {
    MainSequence,
    Giant,
    WhiteDwarf,
}

/// The address of a made-up star: its cube, and its place in the list of the cube. Nothing is
/// stored under it: the same address always gives the same star, computed again.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct StarId {
    pub level: u8,
    pub cell: [i64; 3],
    pub index: u32,
}

impl StarId {
    /// A number of its own for the star, for everything that is made up about it.
    pub fn seed(&self, model: &Model) -> u64 {
        let mut rng = Rng::new(model.seed ^ 0x5741_5253, self.level as usize, self.cell);
        (rng.next() * (1u64 << 53) as f64) as u64 ^ (self.index as u64 + 1).wrapping_mul(0x9E37_79B9_7F4A_7C15)
    }
}

/// One star: position (galactic, parsecs), absolute magnitude, color, what it is and its address.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Star {
    pub position: DVec3,
    pub abs_mag: f32,
    pub color: [f32; 3],
    pub kind: Kind,
    pub id: StarId,
}

/// A star as the field keeps it, in 20 bytes: position (galactic, parsecs, from the corner of its
/// cube), absolute magnitude and color. Good for drawing; `generate` gives the exact place again
/// when it is needed.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Kept {
    pub position: [f32; 3],
    pub abs_mag: f32,
    pub color: [u8; 4],
}

impl Kept {
    fn of(star: &Star, corner: DVec3) -> Kept {
        let p = (star.position - corner).as_vec3();
        Kept { position: [p.x, p.y, p.z], abs_mag: star.abs_mag, color: [star.color[0], star.color[1], star.color[2], 1.0].map(|c| (c.clamp(0.0, 1.0) * 255.0 + 0.5) as u8) }
    }

    /// Where the star is, given the corner of its cube.
    pub fn place(&self, corner: DVec3) -> DVec3 {
        corner + DVec3::new(self.position[0] as f64, self.position[1] as f64, self.position[2] as f64)
    }
}

/// The shape of the galaxy: how many stars there are at each place. The numbers come from the
/// map. Lengths are parsecs.
#[derive(Clone, Debug, PartialEq)]
pub struct Model {
    /// Changes which stars come out, not how many or where on average.
    pub seed: u64,
    /// The galactic centre.
    pub centre: DVec3,
    /// Stars per cubic parsec in the disc at the Sun's distance from the centre.
    pub local_density: f64,
    /// Distances over which the disc thins out by a factor e, outwards and upwards.
    pub disc_length: f64,
    pub disc_height: f64,
    /// Radii between which the disc ends.
    pub disc_edge: [f64; 2],
    pub arm_count: f64,
    /// How much denser the arms are than the space between them, 0 for no arms.
    pub arm_strength: f64,
    /// Tangent of the pitch angle: how tightly the arms are wound.
    pub arm_pitch: f64,
    /// Radius where the arms begin.
    pub arm_start: f64,
    /// Stars per cubic parsec in the middle of the bulge.
    pub bulge_density: f64,
    pub bulge_radius: f64,
    /// How much shorter the bulge is than it is wide, squared.
    pub bulge_flattening: f64,
    /// How bright the stars are: absolute magnitudes with the share of all stars that are
    /// brighter, from bright to faint; the last share is 1. Between two entries the share grows
    /// evenly on a logarithmic scale. The entries of the default are where the spectral classes
    /// meet on the main sequence.
    pub brightness: Vec<(f64, f64)>,
    /// Share of the stars that are white dwarfs, taken from the faint ones.
    pub white_dwarfs: f64,
    /// How often stars have companions, by the mass of the first star, lightest first.
    pub companions: Vec<crate::systems::Multiplicity>,
    /// How many planets stars have.
    pub planets: crate::planets::Abundance,
    /// Other galaxies, whose stars are spread as the light of a grid made from a picture.
    pub galaxies: Vec<GridGalaxy>,
    /// Other galaxies as formulas, fitted to their photographs (fitted.rs).
    pub fitted: Vec<crate::fitted::ParamGalaxy>,
    /// The Milky Way's own stars spread by three pictures instead of the disc, arms and bulge
    /// above (the map's starField "photo": tools/galaxy_xyz.py from R. Hurt's map seen from
    /// above, its edge views made from it and the sky's panorama). None: the formula.
    pub home: Option<HomeVolume>,
}

/// The Milky Way from three pictures (photo_volume.rs), as many stars per cubic parsec at the
/// Sun as `local_density` says, the measured number there: from the Sun out its stars follow the
/// pictures, scaled once so that the mean in the plane within NEAR_SUN_PC of the Sun is it.
#[derive(Clone)]
pub struct HomeVolume {
    pub volume: Arc<crate::photo_volume::PhotoVolume>,
    /// Galactic axes (from the galactic centre) to the pictures' own: their X, the face view's
    /// columns, is toward longitude 270 (its map has longitude 90 on the left), their Z, its
    /// rows down, away from the centre toward the Sun, their Y the galactic north.
    pub to_box: DMat3,
    /// Its stars per cubic parsec times this are systems per cubic parsec.
    pub scale: f64,
    /// The thick disc, which the edge views do not have (drawn, a disc's light falls off faster
    /// than its old stars do): this share of each column's stars, spread in height as
    /// exp(-|y| / thick_height), the rest as the pictures have it.
    pub thick_share: f64,
    pub thick_height: f64,
}

/// The Sun's neighbourhood the Milky Way's pictures are scaled to: within this of it in the plane.
pub const NEAR_SUN_PC: f64 = 500.0;

impl HomeVolume {
    pub fn new(volume: Arc<crate::photo_volume::PhotoVolume>, centre: DVec3, local_density: f64, thick_share: f64, thick_height: f64) -> HomeVolume {
        let to_box = DMat3::from_diagonal(DVec3::new(-1.0, 1.0, -1.0));
        let unscaled = HomeVolume { volume: volume.clone(), to_box, scale: 1.0, thick_share, thick_height };
        let sun = to_box * (DVec3::ZERO - centre);
        // the mean in the plane round the Sun (its own mid-plane rows, within a pixel of it)
        let (mut sum, mut n) = (0.0, 0usize);
        let step = volume.pc.max(10.0);
        let k = (NEAR_SUN_PC / step) as i64;
        for i in -k..=k {
            for j in -k..=k {
                let (dx, dz) = (i as f64 * step, j as f64 * step);
                if dx * dx + dz * dz > NEAR_SUN_PC * NEAR_SUN_PC {
                    continue;
                }
                sum += unscaled.at(sun + DVec3::new(dx, 0.0, dz));
                n += 1;
            }
        }
        let mean = sum / n.max(1) as f64;
        let scale = if mean > 0.0 { local_density / mean } else { 0.0 };
        HomeVolume { volume, to_box, scale, thick_share, thick_height }
    }

    /// Stars per cubic parsec at a place of the pictures' own axes, unscaled: the thin disc as
    /// they have it, the thick one from their columns.
    fn at(&self, q: DVec3) -> f64 {
        let thin = self.volume.density(q) * (1.0 - self.thick_share);
        if self.thick_share <= 0.0 {
            return thin;
        }
        thin + self.thick_share * self.volume.column(q) * (-q.y.abs() / self.thick_height).exp() / (2.0 * self.thick_height)
    }

    pub fn density(&self, q_from_centre: DVec3) -> f64 {
        self.at(self.to_box * q_from_centre) * self.scale
    }

    /// Systems in all: its pictures' stars, scaled.
    pub fn systems(&self) -> f64 {
        self.volume.total() * self.scale
    }
}

impl PartialEq for HomeVolume {
    fn eq(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.volume, &other.volume) && self.scale == other.scale
    }
}

impl std::fmt::Debug for HomeVolume {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "HomeVolume({:.3e} systems, scale {:.3e})", self.systems(), self.scale)
    }
}

/// A galaxy whose stars follow a grid of light (a volume of the map, such as the one
/// tools/photo_galaxy.py makes from a photograph): where the grid is brighter there are more
/// stars. Only the grid and a count of stars are needed: no star is stored.
#[derive(Clone, PartialEq)]
pub struct GridGalaxy {
    pub name: String,
    /// Systems in all, spread over the grid.
    pub systems: f64,
    /// Middle of the box (galactic, parsecs) and the rows that turn galactic axes into its own.
    pub centre: DVec3,
    pub to_box: DMat3,
    /// The box along its own axes, parsecs.
    pub lo: DVec3,
    pub hi: DVec3,
    /// Cells along the three axes of the file; the third runs fastest.
    pub size: [usize; 3],
    /// Systems per cubic parsec in each cell.
    pub density: Arc<Vec<f32>>,
    /// A round bulge besides the grid, which a photograph cannot give (its middle is burnt
    /// out): systems per cubic parsec in its middle, and the radius where that falls by e.
    pub bulge_density: f64,
    pub bulge_radius: f64,
    /// The colour of the light of each cell, from the picture.
    pub colour: Arc<Vec<[u8; 3]>>,
    /// Absolute magnitude of the whole galaxy, for drawing its glow.
    pub abs_mag: f64,
    /// Its stars spread by three pictures, one seen along each of its axes, instead of the grid
    /// (photo_volume.rs: the map's "photo"); `systems` is still how many there are. Read only
    /// while the camera is near enough for its stars to be seen (`Field::update_moving`): until
    /// then it has none.
    pub photo: Option<Arc<crate::photo_volume::LazyVolume>>,
    /// The finer face-on map of the disc (tools/photo_plane.py), as a factor on the grid: within
    /// each cell the stars follow the photograph, lanes and all. Empty if there is none.
    pub plane: Arc<Plane>,
    /// Drawn (its light, galaxy_look.gd) but with no stars of its own: the Milky Way's pictures,
    /// whose stars are the model's (Model::home), and which has no address of its own.
    pub light_only: bool,
    /// The volume (map "volumes", by name) whose dust grid dims its light where it is drawn as
    /// voxels (galaxy_voxels.gdshader): the map's "dust" of the galaxy. None: no dust.
    pub dust: Option<String>,
}

/// A factor on the density of a galaxy, across its disc: rows along the box's first axis,
/// columns along its third, over the box's extent in those axes.
#[derive(Clone, PartialEq, Default)]
pub struct Plane {
    pub n: usize,
    pub factor: Vec<f32>,
    /// For each cell of the grid's plane: the factor on average, and at most.
    pub cell_mean: Vec<f32>,
    pub cell_max: Vec<f32>,
}

impl Plane {
    /// The factor that makes the light of the grid, seen face on, that of the map: the map over
    /// the grid's light summed through the disc, both on the map's pixels. (tools/photo_plane.py
    /// fits the map so that the stars, seen from the Earth, look like the photograph.)
    fn load(file: &std::path::Path, density: &[f32], size: [usize; 3]) -> Result<Plane, String> {
        let img = image::open(file).map_err(|e| format!("{}: {e}", file.display()))?.to_luma16();
        let (w, h) = img.dimensions();
        if w != h {
            return Err(format!("{}: the map must be square", file.display()));
        }
        let n = w as usize;
        let map: Vec<f64> = img.pixels().map(|p| p.0[0] as f64 / 65535.0).collect();
        // the grid seen face on: its light summed through the disc, cell by cell
        let (s0, s2) = (size[0], size[2]);
        let mut column = vec![0.0f64; s0 * s2];
        for (i, &d) in density.iter().enumerate() {
            column[(i / (size[1] * s2)) * s2 + i % s2] += d as f64;
        }
        let grid_at = |x: f64, z: f64| {
            // x, z in 0..1 across the box; the cells blended
            let (fx, fz) = (x * s0 as f64 - 0.5, z * s2 as f64 - 0.5);
            let (x0, z0) = (fx.floor(), fz.floor());
            let (tx, tz) = (fx - x0, fz - z0);
            let at = |a: f64, b: f64| column[(a.clamp(0.0, s0 as f64 - 1.0) as usize) * s2 + b.clamp(0.0, s2 as f64 - 1.0) as usize];
            at(x0, z0) * (1.0 - tx) * (1.0 - tz) + at(x0 + 1.0, z0) * tx * (1.0 - tz) + at(x0, z0 + 1.0) * (1.0 - tx) * tz + at(x0 + 1.0, z0 + 1.0) * tx * tz
        };
        let grid: Vec<f64> = (0..n * n).map(|k| grid_at(((k / n) as f64 + 0.5) / n as f64, ((k % n) as f64 + 0.5) / n as f64)).collect();
        let ratio = grid.iter().sum::<f64>() / map.iter().sum::<f64>().max(1e-30);
        let factor: Vec<f32> = (0..n * n)
            .map(|k| {
                (if grid[k] > 1e-9 { (map[k] * ratio / grid[k]).min(30.0) } else { 0.0 }) as f32
            })
            .collect();
        let (mut cell_mean, mut cell_max) = (vec![0.0f32; s0 * s2], vec![0.0f32; s0 * s2]);
        let mut count = vec![0u32; s0 * s2];
        for (k, &f) in factor.iter().enumerate() {
            let c = ((k / n) * s0 / n) * s2 + (k % n) * s2 / n;
            cell_mean[c] += f;
            cell_max[c] = cell_max[c].max(f);
            count[c] += 1;
        }
        cell_mean.iter_mut().zip(&count).for_each(|(m, &c)| *m /= c.max(1) as f32);
        Ok(Plane { n, factor, cell_mean, cell_max })
    }
}

impl std::fmt::Debug for GridGalaxy {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "GridGalaxy({}, {:e} systems, {:?} cells, bulge {} /pc3 over {} pc)", self.name, self.systems, self.size, self.bulge_density, self.bulge_radius)
    }
}

impl GridGalaxy {
    /// From a light grid of RGBA bytes (alpha: the light, to the power 1/gamma) and its metadata
    /// ("Grid size:" and "Volume bounds:", in units of `size_pc`). `bulge` is the share of the
    /// systems in the round bulge and its radius; the rest follows the grid.
    #[allow(clippy::too_many_arguments)]
    pub fn load(name: &str, light: &std::path::Path, metadata: &std::path::Path, size_pc: f64, centre: DVec3, axes: [DVec3; 3], gamma: f64, systems: f64, bulge: (f64, f64), plane: Option<&std::path::Path>) -> Result<GridGalaxy, String> {
        use std::io::Read;
        let meta = std::fs::read_to_string(metadata).map_err(|e| format!("{}: {e}", metadata.display()))?;
        let numbers = |key: &str| -> Vec<f64> {
            meta.lines()
                .find_map(|l| l.strip_prefix(key))
                .map(|l| l.split(|c: char| !(c.is_ascii_digit() || c == '.' || c == '-' || c == 'e')).filter_map(|t| t.parse().ok()).collect())
                .unwrap_or_default()
        };
        let (grid, b) = (numbers("Grid size:"), numbers("Volume bounds:"));
        if grid.len() != 3 || b.len() != 6 {
            return Err(format!("{}: grid size or bounds not found", metadata.display()));
        }
        let size = [grid[0] as usize, grid[1] as usize, grid[2] as usize];
        let cells = size[0] * size[1] * size[2];
        let file = std::fs::File::open(light).map_err(|e| format!("{}: {e}", light.display()))?;
        let mut bytes = Vec::with_capacity(cells * 4);
        flate2::read::GzDecoder::new(file).read_to_end(&mut bytes).map_err(|e| format!("{}: {e}", light.display()))?;
        if bytes.len() != cells * 4 {
            return Err(format!("{}: {} bytes, expected {} for a grid of {size:?}", light.display(), bytes.len(), cells * 4));
        }
        let (lo, hi) = (DVec3::new(b[0], b[1], b[2]) * size_pc, DVec3::new(b[3], b[4], b[5]) * size_pc);
        let cell_volume = ((hi - lo) / DVec3::new(size[0] as f64, size[1] as f64, size[2] as f64)).element_product();
        let mut density: Vec<f32> = bytes.chunks_exact(4).map(|v| (v[3] as f64 / 255.0).powf(gamma) as f32).collect();
        let colour: Vec<[u8; 3]> = bytes.chunks_exact(4).map(|v| [v[0], v[1], v[2]]).collect();
        let sum: f64 = density.iter().map(|&d| d as f64).sum();
        if sum <= 0.0 {
            return Err(format!("{}: the grid has no light", light.display()));
        }
        let (share, radius) = (bulge.0.clamp(0.0, 1.0), bulge.1.max(1.0));
        let plane = match plane {
            Some(file) => Plane::load(file, &density, size)?,
            None => Plane::default(),
        };
        // the systems in all, with the map's factor (which moves stars about) taken in
        let sum = if plane.n > 0 {
            let mut with = 0.0;
            for (i, &d) in density.iter().enumerate() {
                with += d as f64 * plane.cell_mean[(i / (size[1] * size[2])) * size[2] + i % size[2]] as f64;
            }
            with
        } else {
            sum
        };
        let scale = systems * (1.0 - share) / sum / cell_volume;
        density.iter_mut().for_each(|d| *d = (*d as f64 * scale) as f32);
        // exp(-(r/R)^2) adds up to pi^(3/2) R^3.
        let bulge_density = systems * share / (std::f64::consts::PI.powf(1.5) * radius.powi(3));
        Ok(GridGalaxy { name: name.to_owned(), systems, centre, to_box: DMat3::from_cols(axes[0], axes[1], axes[2]).transpose(), lo, hi, size, density: Arc::new(density), bulge_density, bulge_radius: radius, colour: Arc::new(colour), abs_mag: -21.0, plane: Arc::new(plane), photo: None, light_only: false, dust: None })
    }

    /// Systems per cubic parsec at a place (galactic, parsecs), the cells blended.
    pub fn density(&self, p: DVec3) -> f64 {
        if self.light_only {
            return 0.0;
        }
        let q = self.to_box * (p - self.centre);
        if let Some(lazy) = &self.photo {
            return lazy.get().map_or(0.0, |v| v.density(q) * self.systems / v.stars);
        }
        let bulge = self.bulge_density * (-q.length_squared() / (self.bulge_radius * self.bulge_radius)).exp();
        if q.cmplt(self.lo).any() || q.cmpge(self.hi).any() {
            return bulge;
        }
        let n = DVec3::new(self.size[0] as f64, self.size[1] as f64, self.size[2] as f64);
        // Cell centres sit at half cells: blend the eight around the place.
        let f = (q - self.lo) / (self.hi - self.lo) * n - 0.5;
        let base = f.floor();
        let t = f - base;
        let mut sum = 0.0;
        for corner in 0..8 {
            let o = [corner & 1, (corner >> 1) & 1, (corner >> 2) & 1];
            let mut w = 1.0;
            let mut index = 0;
            for k in 0..3 {
                let i = (base[k] as i64 + o[k] as i64).clamp(0, self.size[k] as i64 - 1) as usize;
                index = index * self.size[k] + i;
                w *= if o[k] == 1 { t[k] } else { 1.0 - t[k] };
            }
            sum += w * self.density[index] as f64;
        }
        sum * self.plane_factor(q) + bulge
    }

    /// The glow of the galaxy as `count` points, each the light of many stars: spread as its stars
    /// are (the grid and the bulge), with the colour of the picture there. Places along the axes
    /// of the box (parsecs, from its middle), colours, and the absolute magnitude of each point.
    pub fn glow(&self, count: usize, seed: u64) -> (Vec<DVec3>, Vec<[f32; 3]>, f64) {
        let mut rng = Rng::seeded(seed ^ 0x6C6F_7767);
        let table = self.cumulative();
        let (mut places, mut colours) = (Vec::with_capacity(count), Vec::with_capacity(count));
        for _ in 0..count {
            let (q, cell) = self.place(&table, &mut rng);
            places.push(q);
            colours.push(match cell {
                None => [1.0, 0.86, 0.7],
                Some(i) => {
                    let c = self.photo.as_ref().and_then(|l| l.get()).map_or_else(|| self.colour[i], |v| v.top_rgb[i]);
                    [c[0] as f32 / 255.0, c[1] as f32 / 255.0, c[2] as f32 / 255.0]
                }
            });
        }
        (places, colours, self.abs_mag + 2.5 * (count.max(1) as f64).log10())
    }

    /// What `place` draws by: how many systems each cell has, added up, and the bulge's share
    /// of all. Made once for many draws (it goes through the whole grid).
    pub fn cumulative(&self) -> (Vec<f64>, f64) {
        if let Some(lazy) = &self.photo {
            return (lazy.load().map(|v| v.cumulative()).unwrap_or_default(), 0.0);
        }
        let n = DVec3::new(self.size[0] as f64, self.size[1] as f64, self.size[2] as f64);
        let cell_volume = ((self.hi - self.lo) / n).element_product();
        let mut cumulative = Vec::with_capacity(self.density.len());
        let mut sum = 0.0;
        for (i, &d) in self.density.iter().enumerate() {
            sum += d as f64 * cell_volume * self.cell_mean(i);
            cumulative.push(sum);
        }
        let in_bulge = self.bulge_density * std::f64::consts::PI.powf(1.5) * self.bulge_radius.powi(3);
        (cumulative, in_bulge / (in_bulge + sum))
    }

    /// A place of one of its stars, drawn as the grid has them (its own axes, parsecs from its
    /// middle), and the cell it is in (None: the bulge). `table` is `cumulative()`.
    pub fn place(&self, table: &(Vec<f64>, f64), rng: &mut Rng) -> (DVec3, Option<usize>) {
        let (cumulative, bulge_share) = table;
        if let Some(v) = self.photo.as_ref().and_then(|l| l.load().ok()) {
            let (q, pixel) = v.place(cumulative, rng);
            return (q, Some(pixel));
        }
        let sum = cumulative.last().copied().unwrap_or(0.0);
        if rng.next() < *bulge_share {
            // exp(-(r/R)^2) is a normal spread of R / sqrt 2 along each axis.
            let mut normal = || (-2.0 * rng.next().max(1e-12).ln()).sqrt() * (std::f64::consts::TAU * rng.next()).cos();
            return (DVec3::new(normal(), normal(), normal()) * self.bulge_radius / 2f64.sqrt(), None);
        }
        let n = DVec3::new(self.size[0] as f64, self.size[1] as f64, self.size[2] as f64);
        let cell = (self.hi - self.lo) / n;
        let u = rng.next() * sum;
        let i = cumulative.partition_point(|&c| c < u).min(cumulative.len() - 1);
        let (i0, rest) = (i / (self.size[1] * self.size[2]), i % (self.size[1] * self.size[2]));
        let (i1, i2) = (rest / self.size[2], rest % self.size[2]);
        // within the cell, as the finer map has it
        let top = self.cell_max(i);
        let mut q = self.lo + DVec3::new(i0 as f64 + rng.next(), i1 as f64 + rng.next(), i2 as f64 + rng.next()) * cell;
        for _ in 0..24 {
            if top <= 0.0 || rng.next() * top <= self.plane_factor(q) {
                break;
            }
            q = self.lo + DVec3::new(i0 as f64 + rng.next(), i1 as f64 + rng.next(), i2 as f64 + rng.next()) * cell;
        }
        (q, Some(i))
    }

    /// The finer map's factor at a place along the box's axes; 1 without a map.
    pub fn plane_factor(&self, q: DVec3) -> f64 {
        let p = &self.plane;
        if p.n == 0 {
            return 1.0;
        }
        let n = p.n as f64;
        let fx = (q.x - self.lo.x) / (self.hi.x - self.lo.x) * n - 0.5;
        let fz = (q.z - self.lo.z) / (self.hi.z - self.lo.z) * n - 0.5;
        let (x0, z0) = (fx.floor(), fz.floor());
        let (tx, tz) = (fx - x0, fz - z0);
        let at = |x: f64, z: f64| p.factor[(x.clamp(0.0, n - 1.0) as usize) * p.n + z.clamp(0.0, n - 1.0) as usize] as f64;
        at(x0, z0) * (1.0 - tx) * (1.0 - tz) + at(x0 + 1.0, z0) * tx * (1.0 - tz) + at(x0, z0 + 1.0) * (1.0 - tx) * tz + at(x0 + 1.0, z0 + 1.0) * tx * tz
    }

    fn cell_mean(&self, i: usize) -> f64 {
        if self.plane.n == 0 { 1.0 } else { self.plane.cell_mean[(i / (self.size[1] * self.size[2])) * self.size[2] + i % self.size[2]] as f64 }
    }

    fn cell_max(&self, i: usize) -> f64 {
        if self.plane.n == 0 { 1.0 } else { self.plane.cell_max[(i / (self.size[1] * self.size[2])) * self.size[2] + i % self.size[2]] as f64 }
    }

    /// Whether a cube (corner, side) could hold any of its stars.
    fn touches(&self, corner: DVec3, side: f64) -> bool {
        let mid = corner + DVec3::splat(side * 0.5);
        let q = self.to_box * (mid - self.centre);
        let r = side * 0.87;
        !(q.cmplt(self.lo - r).any() || q.cmpgt(self.hi + r).any())
    }
}

/// The spectral classes, each with the absolute magnitude where it ends on the main sequence.
pub const CLASSES: [(char, f64); 7] = [('O', -3.9), ('B', 0.99), ('A', 2.57), ('F', 4.48), ('G', 5.78), ('K', 8.8), ('M', 19.75)];
/// White dwarfs are this faint.
const WHITE_DWARF_MAGS: (f64, f64) = (10.0, 15.5);

impl Default for Model {
    fn default() -> Model {
        Model {
            seed: 1,
            centre: DVec3::new(0.0, 0.0, 8000.0),
            local_density: 0.14,
            disc_length: 2100.0,
            disc_height: 350.0,
            disc_edge: [16_000.0, 20_000.0],
            arm_count: 2.0,
            arm_strength: 0.7,
            arm_pitch: 0.23,
            arm_start: 3000.0,
            bulge_density: 21.6,
            bulge_radius: 1200.0,
            bulge_flattening: 2.0,
            // Shares of the classes O, B, A, F, G, K, M among the stars near the Sun: see map/SOURCES.md.
            brightness: vec![(BRIGHTEST, 1e-10), (-3.9, 3e-7), (0.99, 0.0012), (2.57, 0.0073), (4.48, 0.0373), (5.78, 0.1133), (8.8, 0.2333), (19.75, 1.0)],
            white_dwarfs: 0.025,
            companions: crate::systems::MULTIPLICITY.to_vec(),
            planets: crate::planets::ABUNDANCE,
            galaxies: Vec::new(),
            fitted: Vec::new(),
            home: None,
        }
    }
}

impl Model {
    /// Reads the pictures of the galaxies made of three (`GridGalaxy::photo`) that a cube with
    /// stars round the camera may reach (`star_reach`, and the reach of their pictures), and
    /// lets go of those far beyond that (half as far again, so that going to and fro at the edge
    /// does not read them over and over). Before that no cube can touch them, so none was made
    /// without their stars. From afar such a galaxy is its photograph (GalaxyLook).
    pub fn keep_near(&self, camera: DVec3) {
        let stars = star_reach();
        for g in &self.galaxies {
            let Some(lazy) = &g.photo else { continue };
            let near = lazy.reach.0.hypot(lazy.reach.1) + stars;
            let d = camera.distance(g.centre);
            if d < near {
                if lazy.get().is_none() {
                    let t = std::time::Instant::now();
                    match lazy.load() {
                        Ok(_) => eprintln!("{}: its pictures read ({:.0} kpc away) in {:?}", lazy.name, d / 1000.0, t.elapsed()),
                        Err(e) => eprintln!("{}: {e}", lazy.name),
                    }
                }
            } else if d > 1.5 * near && lazy.get().is_some() {
                lazy.release();
                eprintln!("{}: its pictures let go ({:.0} kpc away)", lazy.name, d / 1000.0);
            }
        }
    }

    /// Stars per cubic parsec at a place.
    pub fn density(&self, p: DVec3) -> f64 {
        let q = p - self.centre;
        if let Some(home) = &self.home {
            return home.density(q) + self.galaxies.iter().map(|g| g.density(p)).sum::<f64>() + self.fitted.iter().map(|g| g.density(p)).sum::<f64>();
        }
        let r = (q.x * q.x + q.z * q.z).sqrt();
        let sun = self.centre.length();
        let mut disc = self.local_density * (-(r - sun) / self.disc_length).exp() * (-q.y.abs() / self.disc_height).exp();
        if r > self.arm_start && self.arm_strength != 0.0 {
            let angle = q.x.atan2(q.z);
            disc *= 1.0 + self.arm_strength * (self.arm_count * (angle - (r / self.arm_start).ln() / self.arm_pitch)).cos();
        }
        let [inner, outer] = self.disc_edge;
        disc *= 1.0 - ((r - inner) / (outer - inner).max(1.0)).clamp(0.0, 1.0);
        let d2 = (q.x * q.x + q.z * q.z + q.y * q.y * self.bulge_flattening) / (self.bulge_radius * self.bulge_radius);
        disc.max(0.0) + self.bulge_density * (-d2).exp() + self.galaxies.iter().map(|g| g.density(p)).sum::<f64>() + self.fitted.iter().map(|g| g.density(p)).sum::<f64>()
    }

    /// Share of all stars that are brighter than an absolute magnitude.
    pub fn brighter_than(&self, m: f64) -> f64 {
        let b = &self.brightness;
        if m <= b[0].0 {
            return b[0].1;
        }
        match b.windows(2).find(|w| m < w[1].0) {
            Some(w) => {
                let t = (m - w[0].0) / (w[1].0 - w[0].0);
                10f64.powf(w[0].1.log10() + (w[1].1.log10() - w[0].1.log10()) * t)
            }
            None => 1.0,
        }
    }

    /// Of all the light of the stars, the share in those fainter than each of `n` absolute
    /// magnitudes from `lo` to `hi` (evenly): what is left as a haze once the stars brighter than
    /// a magnitude are drawn one by one. 1 at the bright end, 0 at the faint one.
    pub fn light_fainter_than(&self, lo: f64, hi: f64, n: usize) -> Vec<f32> {
        // the stars by their share brighter (u), on a logarithmic scale: the bright few hold
        // much of the light
        let first = self.brightness[0].1.max(1e-15);
        let steps = 40_000;
        let mut mags = Vec::with_capacity(steps);
        let mut lights = Vec::with_capacity(steps);
        for i in 0..steps {
            let (a, b) = (first * (1.0 / first).powf(i as f64 / steps as f64), first * (1.0 / first).powf((i + 1) as f64 / steps as f64));
            let m = self.mag_of_share((a * b).sqrt());
            mags.push(m);
            lights.push((b - a) * 10f64.powf(-0.4 * m));
        }
        let total: f64 = lights.iter().sum::<f64>().max(1e-300);
        (0..n)
            .map(|k| {
                let at = lo + (hi - lo) * k as f64 / (n - 1).max(1) as f64;
                (mags.iter().zip(&lights).filter(|(m, _)| **m > at).map(|(_, l)| l).sum::<f64>() / total) as f32
            })
            .collect()
    }

    /// The inverse of `brighter_than`.
    pub fn mag_of_share(&self, u: f64) -> f64 {
        let b = &self.brightness;
        match b.windows(2).find(|w| u < w[1].1) {
            Some(w) if u > w[0].1 => w[0].0 + (w[1].0 - w[0].0) * (u.log10() - w[0].1.log10()) / (w[1].1.log10() - w[0].1.log10()),
            Some(w) => w[0].0,
            None => b[b.len() - 1].0,
        }
    }

    /// A star for a draw `u` of the share of stars brighter than it, and two more draws: its
    /// absolute magnitude and what it is.
    pub fn star_of_share(&self, u: f64, draw_giant: f64, draw_dwarf: f64) -> (f64, Kind) {
        let m = self.mag_of_share(u);
        if draw_giant < 0.3 && (-1.5..2.0).contains(&m) {
            return (m, Kind::Giant);
        }
        // White dwarfs take the place of some of the stars as faint as they are.
        let (from, to) = WHITE_DWARF_MAGS;
        if (from..to).contains(&m) && draw_dwarf * (self.brighter_than(to) - self.brighter_than(from)) < self.white_dwarfs {
            return (m, Kind::WhiteDwarf);
        }
        (m, Kind::MainSequence)
    }

    /// How many stars a system has on average: a sample of systems, counted.
    pub fn stars_per_system(&self) -> f64 {
        let mut rng = Rng::seeded(self.seed);
        let n = 20_000;
        let mut stars = 0;
        for i in 0..n {
            let (m, kind) = self.star_of_share(rng.next(), rng.next(), rng.next());
            let star = Star { position: DVec3::ZERO, abs_mag: m as f32, color: [1.0; 3], kind, id: StarId { level: 0, cell: [0; 3], index: i } };
            stars += crate::systems::system_of(&star, rng.next().to_bits(), &self.companions).stars();
        }
        stars as f64 / n as f64
    }

    /// How many systems (single stars, or stars with their companions) the model has in all:
    /// the density added up over the galaxy.
    pub fn total(&self) -> f64 {
        const STEP: f64 = 50.0;
        const HEIGHT_STEP: f64 = 20.0;
        const TURN: usize = 32;
        let reach = self.disc_edge[1].max(6.0 * self.bulge_radius);
        let mut sum = 0.0;
        for i in 0..(reach / STEP) as usize {
            let r = (i as f64 + 0.5) * STEP;
            for j in -300..300 {
                let y = (j as f64 + 0.5) * HEIGHT_STEP;
                for k in 0..TURN {
                    let angle = (k as f64 + 0.5) / TURN as f64 * std::f64::consts::TAU;
                    let p = self.centre + DVec3::new(r * angle.sin(), y, r * angle.cos());
                    sum += self.density(p) * r * STEP * HEIGHT_STEP * std::f64::consts::TAU / TURN as f64;
                }
            }
        }
        sum
    }
}

pub fn color(m: f64, kind: Kind) -> [f32; 3] {
    match kind {
        Kind::Giant => return [1.0, 0.78, 0.55],
        Kind::WhiteDwarf => return [0.8, 0.87, 1.0],
        Kind::MainSequence => {}
    }
    // From blue giants to red dwarfs.
    const STOPS: [(f64, [f32; 3]); 6] =
        [(-3.0, [0.62, 0.74, 1.0]), (0.0, [0.8, 0.87, 1.0]), (2.5, [1.0, 0.98, 0.95]), (5.0, [1.0, 0.93, 0.8]), (7.5, [1.0, 0.82, 0.62]), (10.0, [1.0, 0.68, 0.45])];
    let i = STOPS.iter().position(|s| m < s.0).unwrap_or(STOPS.len());
    if i == 0 {
        return STOPS[0].1;
    }
    if i == STOPS.len() {
        return STOPS[i - 1].1;
    }
    let ((m0, c0), (m1, c1)) = (STOPS[i - 1], STOPS[i]);
    let t = ((m - m0) / (m1 - m0)) as f32;
    [0, 1, 2].map(|k| c0[k] + (c1[k] - c0[k]) * t)
}

/// A small random number generator, seeded from the cube.
pub struct Rng(u64);

impl Rng {
    /// From a number alone, such as the one of a star.
    pub fn seeded(seed: u64) -> Rng {
        let mut rng = Rng((seed ^ 0x9E37_79B9_7F4A_7C15).wrapping_mul(0xBF58_476D_1CE4_E5B9) | 1);
        rng.next();
        rng
    }

    fn new(seed: u64, level: usize, cell: [i64; 3]) -> Rng {
        let mut h = 0x9E37_79B9_7F4A_7C15u64 ^ level as u64 ^ seed.wrapping_mul(0x94D0_49BB_1331_11EB);
        for c in cell {
            h = (h ^ c as u64).wrapping_mul(0xBF58_476D_1CE4_E5B9);
            h ^= h >> 31;
        }
        Rng(h | 1)
    }

    /// Uniform in [0, 1).
    pub fn next(&mut self) -> f64 {
        self.0 ^= self.0 >> 12;
        self.0 ^= self.0 << 25;
        self.0 ^= self.0 >> 27;
        (self.0.wrapping_mul(0x2545_F491_4F6C_DD1D) >> 11) as f64 / (1u64 << 53) as f64
    }
}

/// The corner of a cube with the smallest coordinates.
fn corner_of(level: usize, cell: [i64; 3]) -> DVec3 {
    DVec3::new(cell[0] as f64, cell[1] as f64, cell[2] as f64) * side(level)
}

/// How far from the camera (parsecs) a cube with stars can be: the largest such cube's reach, the
/// way ahead of it (`update_moving`: up to AHEAD_REACHES of its reach) and its own half diagonal.
pub fn star_reach() -> f64 {
    (0..LEVELS)
        .filter(|&level| faintest(level) > faintest(level + 1).max(BRIGHTEST))
        .map(|level| {
            let s = side(level);
            (REACH * s + 0.87 * s) * (1.0 + AHEAD_REACHES) + 0.87 * s
        })
        .fold(0.0, f64::max)
}

fn side(level: usize) -> f64 {
    SMALLEST_PC * (1u64 << level) as f64
}

/// Faintest absolute magnitude that can be seen from `REACH` sides of a cube of this level.
/// The smallest cubes hold every star left, however faint: next to the camera all stars exist.
fn faintest(level: usize) -> f64 {
    if level == 0 { f64::INFINITY } else { LIMIT_MAG + 5.0 - 5.0 * (REACH * side(level)).log10() }
}

/// A cube's magnitudes (the shares of stars brighter than its brightest and its faintest), its
/// corner and side, and how many stars of them it holds on average, from its density sampled over
/// it, with the most of that density and the generator the samples left (for its stars). None if
/// it holds none.
struct Sampled {
    share_bright: f64,
    share_faint: f64,
    corner: DVec3,
    side: f64,
    mean: f64,
    max: f64,
    rng: Rng,
}

fn sample(model: &Model, level: usize, cell: [i64; 3], part: Option<(usize, [usize; 3])>) -> Option<Sampled> {
    let whole = side(level);
    // This level holds the stars too faint for the next one up.
    let (bright, faint) = (faintest(level + 1).max(BRIGHTEST), faintest(level));
    if faint <= bright {
        return None;
    }
    let (share_bright, share_faint) = (model.brighter_than(bright), model.brighter_than(faint));
    let cell_corner = DVec3::new(cell[0] as f64, cell[1] as f64, cell[2] as f64) * whole;
    // Far from the Milky Way and from every other galaxy there is nothing to sample.
    let home = (cell_corner + DVec3::splat(whole * 0.5)).distance(model.centre) - whole * 0.87;
    if home > model.disc_edge[1].max(8.0 * model.bulge_radius) && !model.galaxies.iter().any(|g| g.touches(cell_corner, whole)) && !model.fitted.iter().any(|g| g.touches(cell_corner, whole)) {
        return None;
    }
    // a part of a cube that has too many stars (Field): a cube of its own, seeded by where it is
    let (corner, s, mut rng) = match part {
        None => (cell_corner, whole, Rng::new(model.seed, level, cell)),
        Some((k, at)) => {
            let ps = whole / k as f64;
            let lin = ((at[2] * k + at[1]) * k + at[0]) as u64;
            (cell_corner + DVec3::new(at[0] as f64, at[1] as f64, at[2] as f64) * ps, ps, Rng::new(model.seed ^ (lin + 1).wrapping_mul(0xD1B5_4A32_D192_ED03), level, cell))
        }
    };

    // How dense the cube is, from samples spread over it. Large cubes are much thicker than the
    // disc, so they need more samples.
    let n = if s > 400.0 { 6 } else { 2 };
    let (mut sum, mut max) = (0.0, 0.0f64);
    for i in 0..n {
        for j in 0..n {
            for k in 0..n {
                let jitter = DVec3::new(i as f64 + rng.next(), j as f64 + rng.next(), k as f64 + rng.next()) / n as f64;
                let d = model.density(corner + jitter * s);
                sum += d;
                max = max.max(d);
            }
        }
    }
    let mean = sum / (n * n * n) as f64 * s * s * s * (share_faint - share_bright);
    (mean > 0.0).then_some(Sampled { share_bright, share_faint, corner, side: s, mean, max, rng })
}

/// The stars of one cube. The same cube always gives the same stars.
pub fn generate(model: &Model, level: usize, cell: [i64; 3]) -> Vec<Star> {
    sample(model, level, cell, None).map_or_else(Vec::new, |c| make(model, level, cell, c, 0))
}

/// Into how many parts each way a cube is cut for its stars near the camera (Field): 1 if it holds
/// no more than can be drawn (MAX_PER_CELL), else parts that each hold about that many.
pub fn parts_of(model: &Model, level: usize, cell: [i64; 3]) -> usize {
    match sample(model, level, cell, None) {
        Some(c) if c.mean > MAX_PER_CELL as f64 => ((c.mean / MAX_PER_CELL as f64).cbrt().ceil() as usize).clamp(2, MOST_PARTS),
        _ => 1,
    }
}

/// The stars of one part of a cube cut into k^3 (`parts_of`): all of its stars of the cube's
/// magnitudes (but for MAX_PER_CELL). The same part always gives the same stars.
pub fn generate_part(model: &Model, level: usize, cell: [i64; 3], k: usize, at: [usize; 3]) -> Vec<Star> {
    let lin = (at[2] * k + at[1]) * k + at[0];
    sample(model, level, cell, Some((k, at))).map_or_else(Vec::new, |c| make(model, level, cell, c, lin as u32 + 1))
}

fn make(model: &Model, level: usize, cell: [i64; 3], c: Sampled, part: u32) -> Vec<Star> {
    let Sampled { share_bright, share_faint, corner, side: s, mean, max, mut rng } = c;
    // A whole number of stars, right on average. A cube with too many to draw keeps its
    // brightest ones: the part of its range of magnitudes that holds as many as can be drawn.
    let share_faint = if mean > MAX_PER_CELL as f64 { share_bright + (share_faint - share_bright) * MAX_PER_CELL as f64 / mean } else { share_faint };
    let mean = mean.min(MAX_PER_CELL as f64);
    let count = (mean.floor() as usize) + (rng.next() < mean.fract()) as usize;

    let mut stars = Vec::with_capacity(count);
    let mut tries = 0;
    while stars.len() < count && tries < count * 40 + 100 {
        tries += 1;
        let p = corner + DVec3::new(rng.next(), rng.next(), rng.next()) * s;
        // More stars where the cube is denser.
        if rng.next() * max * 1.5 > model.density(p) {
            continue;
        }
        let (m, kind) = model.star_of_share(share_bright + rng.next() * (share_faint - share_bright), rng.next(), rng.next());
        // The catalog has the stars that are bright from the Sun.
        if m + 5.0 * (p.length().max(1e-3) / 10.0).log10() < CATALOG_MAG {
            continue;
        }
        let id = StarId { level: level as u8, cell, index: part << PART_SHIFT | stars.len() as u32 };
        stars.push(Star { position: p, abs_mag: m as f32, color: color(m, kind), kind, id });
    }
    stars
}

/// The star at an address, computed again. None if the address is of no star.
pub fn find(model: &Model, id: StarId) -> Option<Star> {
    let (level, part, i) = (id.level as usize, id.index >> PART_SHIFT, (id.index & ((1 << PART_SHIFT) - 1)) as usize);
    if part == 0 {
        return generate(model, level, id.cell).get(i).copied();
    }
    let k = parts_of(model, level, id.cell);
    let lin = part as usize - 1;
    if lin >= k * k * k {
        return None;
    }
    generate_part(model, level, id.cell, k, [lin % k, lin / k % k, lin / (k * k)]).get(i).copied()
}

/// The cubes that exist now, around the camera.
#[derive(Default)]
pub struct Field {
    model: Model,
    /// By level, cube and part (0: the whole cube, else 1 + its index among k^3): its stars.
    cells: HashMap<(usize, [i64; 3], u32), Vec<Kept>>,
    /// Into how many parts each way each cube in reach is cut (`parts_of`).
    parts: HashMap<(usize, [i64; 3]), usize>,
    /// The cubes cut into parts near the camera: how far round it (parsecs) their parts are made.
    near: HashMap<(usize, [i64; 3]), f64>,
}

impl Field {
    pub fn new(model: Model) -> Field {
        Field { model, cells: HashMap::new(), parts: HashMap::new(), near: HashMap::new() }
    }

    /// Brings the cubes in line with the camera position. Computes at most `budget` new cubes,
    /// the nearest to the camera first; call again to go on. Returns whether anything changed.
    pub fn update(&mut self, camera: DVec3, budget: usize) -> bool {
        self.update_moving(camera, DVec3::ZERO, budget)
    }

    /// ... the camera going at `velocity` (parsecs a second): each level's cubes are those within
    /// its reach not of the camera alone (a ball) but of the way it goes in the next LOOK_AHEAD_S
    /// (a long capsule, the faster the longer): the stars it is coming to are made before it is
    /// there (they take a while to make, and at thousands of light years a second the ball was
    /// always behind). No longer than AHEAD_REACHES of the level's reach: the small cubes' stars,
    /// near the ship, are passed in a moment; a long way of them would be thousands of cubes.
    pub fn update_moving(&mut self, camera: DVec3, velocity: DVec3, budget: usize) -> bool {
        self.model.keep_near(camera);
        let mut wanted: Vec<(usize, [i64; 3])> = Vec::new();
        let ahead_dir = velocity.normalize_or_zero();
        for level in (0..LEVELS).rev() {
            let s = side(level);
            if faintest(level) <= faintest(level + 1).max(BRIGHTEST) {
                continue;
            }
            let reach = REACH * s + 0.87 * s;
            let length = (velocity.length() * LOOK_AHEAD_S).min(AHEAD_REACHES * reach);
            let end = camera + ahead_dir * length;
            // the cubes round the way from the camera to `end`
            let lo = camera.min(end) - DVec3::splat(reach);
            let hi = camera.max(end) + DVec3::splat(reach);
            let cell_of = |x: f64| (x / s).floor() as i64;
            for cx in cell_of(lo.x)..=cell_of(hi.x) {
                for cy in cell_of(lo.y)..=cell_of(hi.y) {
                    for cz in cell_of(lo.z)..=cell_of(hi.z) {
                        let cell = [cx, cy, cz];
                        let mid = (DVec3::new(cx as f64, cy as f64, cz as f64) + 0.5) * s;
                        let along = (mid - camera).dot(ahead_dir).clamp(0.0, length);
                        if mid.distance(camera + ahead_dir * along) <= reach {
                            wanted.push((level, cell));
                        }
                    }
                }
            }
        }
        // a cube with too many stars, near the camera: its parts near it instead
        let in_reach: std::collections::HashSet<_> = wanted.iter().copied().collect();
        self.parts.retain(|k, _| in_reach.contains(k));
        let mut sampled = 0;
        let mut keys: Vec<(usize, [i64; 3], u32)> = Vec::with_capacity(wanted.len());
        self.near.clear();
        // nearest the camera first (sampled and made): the cubes there are those that would be
        // seen come in late, bright beside it, if the far ones went first
        let box_from = |corner: DVec3, s: f64| camera.clamp(corner, corner + DVec3::splat(s)).distance(camera);
        let mut wanted = wanted;
        wanted.sort_by(|a, b| box_from(corner_of(a.0, a.1), side(a.0)).total_cmp(&box_from(corner_of(b.0, b.1), side(b.0))));
        for (level, cell) in wanted {
            let k = match self.parts.get(&(level, cell)) {
                Some(&k) => k,
                None if sampled < budget * 8 => {
                    sampled += 1;
                    let k = parts_of(&self.model, level, cell);
                    self.parts.insert((level, cell), k);
                    k
                }
                None => 1,
            };
            if k == 1 {
                keys.push((level, cell, 0));
                continue;
            }
            let ps = side(level) / k as f64;
            let within = NEAR_PARTS * ps + 0.87 * ps;
            let local = (camera - corner_of(level, cell)) / ps;
            let span = within / ps;
            let lo = |x: f64| ((x - span).floor().max(0.0) as usize).min(k);
            let hi = |x: f64| (((x + span).ceil()).max(0.0) as usize).min(k);
            let mut near = Vec::new();
            for z in lo(local.z)..hi(local.z) {
                for y in lo(local.y)..hi(local.y) {
                    for x in lo(local.x)..hi(local.x) {
                        let mid = (DVec3::new(x as f64, y as f64, z as f64) + 0.5) * ps + corner_of(level, cell);
                        if mid.distance(camera) <= within {
                            near.push((level, cell, ((z * k + y) * k + x) as u32 + 1));
                        }
                    }
                }
            }
            // the whole cube (its brightest, spread thin) is kept too: beyond the parts it shows,
            // nearer it gives way to them, softly (the stars of each faded by their distance from
            // the camera, `near_radius`), so that no star comes in or goes out at once as the
            // parts round the camera change
            keys.push((level, cell, 0));
            if !near.is_empty() {
                // (faded within the radius the parts cover whole, not to the last part's middle)
                self.near.insert((level, cell), NEAR_PARTS * ps);
                keys.extend(near);
            }
        }
        let part_from = |&(level, cell, part): &(usize, [i64; 3], u32)| {
            if part == 0 {
                // a whole cube cut into parts near the camera: its stars spread thin over all of it,
                // seen beyond its parts; made after them (as far as its parts reach), else its far
                // stars came in first and the near ones after
                let near = self.near.get(&(level, cell)).copied().unwrap_or(0.0);
                return box_from(corner_of(level, cell), side(level)).max(near);
            }
            let k = self.parts.get(&(level, cell)).copied().unwrap_or(1);
            let ps = side(level) / k as f64;
            let lin = part as usize - 1;
            let at = DVec3::new((lin % k) as f64, (lin / k % k) as f64, (lin / (k * k)) as f64);
            box_from(corner_of(level, cell) + at * ps, ps)
        };
        keys.sort_by(|a, b| part_from(a).total_cmp(&part_from(b)));
        let before = self.cells.len();
        let keep: std::collections::HashSet<_> = keys.iter().copied().collect();
        self.cells.retain(|k, _| keep.contains(k));
        let mut changed = self.cells.len() != before;
        let mut made = 0;
        for key in keys {
            if made >= budget {
                break;
            }
            if !self.cells.contains_key(&key) {
                let (level, cell, part) = key;
                let corner = corner_of(level, cell);
                let stars = if part == 0 {
                    generate(&self.model, level, cell)
                } else {
                    let k = self.parts[&(level, cell)];
                    let lin = part as usize - 1;
                    generate_part(&self.model, level, cell, k, [lin % k, lin / k % k, lin / (k * k)])
                };
                self.cells.insert(key, stars.iter().map(|s| Kept::of(s, corner)).collect());
                made += 1;
                changed = true;
            }
        }
        changed || sampled > 0
    }

    /// Into how many parts each way a cube in reach is cut (1: whole; `parts_of`).
    pub fn parts_of_cube(&self, level: usize, cell: [i64; 3]) -> usize {
        self.parts.get(&(level, cell)).copied().unwrap_or(1)
    }

    /// The cubes that exist now: each its address (level, cube, part: 0 the whole cube), its
    /// corner (galactic, parsecs) and its stars.
    pub fn cubes(&self) -> impl Iterator<Item = ((usize, [i64; 3], u32), DVec3, &[Kept])> {
        self.cells.iter().map(|(&key, stars)| (key, corner_of(key.0, key.1), stars.as_slice()))
    }

    /// Into how many parts each way a cube in reach is cut (1: whole; `parts_of`).
    pub fn parts_of(&self, level: usize, cell: [i64; 3]) -> usize {
        self.parts.get(&(level, cell)).copied().unwrap_or(1)
    }

    /// How far round the camera a cube's parts are made (parsecs), if it is cut into parts near
    /// it now: its parts' stars fade out toward there, its whole cube's fade in beyond.
    pub fn near_radius(&self, level: usize, cell: [i64; 3]) -> Option<f64> {
        self.near.get(&(level, cell)).copied()
    }

    /// The stars that exist now, each with its place (galactic, parsecs).
    pub fn stars(&self) -> impl Iterator<Item = (DVec3, &Kept)> {
        self.cells.iter().flat_map(|(&(level, cell, _), stars)| {
            let corner = corner_of(level, cell);
            stars.iter().map(move |s| (s.place(corner), s))
        })
    }

    /// The addresses of the stars within a distance of a place.
    pub fn within(&self, place: DVec3, radius: f64) -> Vec<StarId> {
        let mut found = Vec::new();
        for (&(level, cell, part), stars) in &self.cells {
            // Only the cubes that reach into the sphere.
            let (corner, s) = (corner_of(level, cell), side(level));
            let nearest = place.clamp(corner, corner + DVec3::splat(s));
            if nearest.distance(place) > radius {
                continue;
            }
            found.extend(stars.iter().enumerate().filter(|(_, k)| k.place(corner).distance(place) <= radius).map(|(i, _)| StarId { level: level as u8, cell, index: part << PART_SHIFT | i as u32 }));
        }
        found
    }

    /// The distance in parsecs to the nearest made-up star that lies ahead: within a cone around
    /// a direction, given by the cosine of its half angle. None if there is no star in the cone.
    pub fn ahead(&self, place: DVec3, dir: DVec3, cos_half_angle: f64) -> Option<f64> {
        let mut best = f64::INFINITY;
        for (&(level, cell, _), stars) in &self.cells {
            // A cube that is farther than the best so far has nothing nearer.
            let (corner, s) = (corner_of(level, cell), side(level));
            if place.clamp(corner, corner + DVec3::splat(s)).distance(place) >= best {
                continue;
            }
            for star in stars {
                let to = star.place(corner) - place;
                let d = to.length();
                if d < best && to.dot(dir) >= cos_half_angle * d {
                    best = d;
                }
            }
        }
        best.is_finite().then_some(best)
    }

    /// The made-up stars nearest to a place, nearest first, with their distances in parsecs.
    pub fn nearest(&self, place: DVec3, count: usize) -> Vec<(StarId, f64)> {
        let mut found: Vec<(StarId, f64)> = self
            .cells
            .iter()
            .flat_map(|(&(level, cell, part), stars)| {
                let corner = corner_of(level, cell);
                stars.iter().enumerate().map(move |(i, s)| (StarId { level: level as u8, cell, index: part << PART_SHIFT | i as u32 }, s.place(corner).distance(place)))
            })
            .collect();
        found.sort_by(|a, b| a.1.total_cmp(&b.1));
        found.truncate(count);
        found
    }

    pub fn model(&self) -> &Model {
        &self.model
    }

    pub fn budget() -> usize {
        CELLS_PER_UPDATE
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Going fast, the cubes ahead are made before the camera is there: the way ahead in the
    /// next LOOK_AHEAD_S is in reach, behind it only the ball round it; round the camera as before.
    #[test]
    fn going_fast_the_way_ahead_is_made() {
        let model = Model::default();
        let camera = model.centre + DVec3::new(0.0, 20.0, -4000.0);
        let v = DVec3::new(0.0, 0.0, 300.0); // pc a second, toward the middle
        let count_near = |field: &Field, at: DVec3| field.stars().filter(|(p, _)| p.distance(at) < 40.0).count();
        let mut still = Field::new(model.clone());
        while still.update(camera, 10_000) {}
        let mut moving = Field::new(model.clone());
        while moving.update_moving(camera, v, 10_000) {}
        let soon = camera + v * 1.0;
        let (a, b) = (count_near(&still, soon), count_near(&moving, soon));
        let all = |f: &Field| f.stars().count();
        eprintln!("300 pc a second: within 40 pc of where it is in a second, still {a}, going {b}; stars in all {} and {}", all(&still), all(&moving));
        assert!(b > a * 2, "the way ahead made: {b} against {a}");
        assert!(count_near(&moving, camera) >= count_near(&still, camera) * 9 / 10, "and round the camera as before");
    }

    #[test]
    fn cubes_are_repeatable() {
        let model = Model::default();
        let a = generate(&model, 5, [10, -1, 120]);
        assert!(!a.is_empty());
        assert_eq!(a, generate(&model, 5, [10, -1, 120]));
        assert_ne!(a, generate(&model, 5, [11, -1, 120]));
        assert_ne!(a, generate(&Model { seed: 2, ..Model::default() }, 5, [10, -1, 120]), "another seed, other stars");
        assert!(generate(&Model { local_density: 0.0, bulge_density: 0.0, ..Model::default() }, 5, [10, -1, 120]).is_empty());
        let s = side(5);
        for star in &a {
            assert!(star.position.x >= 10.0 * s && star.position.x < 11.0 * s);
            assert!(star.position.y >= -s && star.position.y < 0.0);
        }
    }

    #[test]
    fn brightness_shares_invert() {
        let model = Model::default();
        for m in [-6.0, -1.0, 4.9, 5.0, 5.1, 9.0, 19.0] {
            assert!((model.mag_of_share(model.brighter_than(m)) - m).abs() < 1e-9, "{m}");
        }
        assert!((model.brighter_than(4.83) - 0.05).abs() < 0.01, "about 5% of stars outshine the Sun");
        assert!(model.brighter_than(-5.0) < 1e-6);
        assert_eq!(model.brighter_than(40.0), 1.0);
    }

    #[test]
    fn the_stars_come_in_the_shares_of_the_table() {
        // A million stars drawn like the cubes draw them, counted by spectral class.
        let model = Model::default();
        let mut rng = Rng::new(7, 0, [0, 0, 0]);
        let n = 1_000_000;
        let (mut by_class, mut giants, mut dwarfs) = ([0u32; 7], 0u32, 0u32);
        for _ in 0..n {
            let (m, kind) = model.star_of_share(rng.next(), rng.next(), rng.next());
            match kind {
                Kind::Giant => giants += 1,
                Kind::WhiteDwarf => dwarfs += 1,
                Kind::MainSequence => by_class[CLASSES.iter().position(|c| m < c.1).unwrap_or(6)] += 1,
            }
        }
        let share = |count: u32| count as f64 / n as f64;
        // M 74% (76.7% less the white dwarfs among them), K 12%, G 7.6%, F 3%.
        for (class, expected, tolerance) in [(6, 0.742, 0.01), (5, 0.120, 0.005), (4, 0.076, 0.004), (3, 0.030, 0.003)] {
            assert!((share(by_class[class]) - expected).abs() < tolerance, "class {}: {}", CLASSES[class].0, share(by_class[class]));
        }
        assert!((share(dwarfs) - model.white_dwarfs).abs() < 0.003, "white dwarfs: {}", share(dwarfs));
        assert!(share(by_class[1]) < 0.002 && share(by_class[2]) < 0.008 && by_class[0] < 5);
        assert!(giants > 0 && share(giants) < 0.003, "giants: {}", share(giants));
    }

    #[test]
    fn the_star_ahead_is_the_nearest_in_its_direction() {
        let mut field = Field::default();
        let camera = DVec3::new(30.0, 5.0, 2000.0);
        while field.update(camera, 10_000) {}
        let (id, distance) = field.nearest(camera, 1)[0];
        let star = find(field.model(), id).unwrap();
        let towards = (star.position - camera).normalize();
        // Looking at the nearest star, it is the one ahead; looking away from it, another, farther one is.
        let ahead = field.ahead(camera, towards, 0.9).unwrap();
        assert!((ahead - distance).abs() < 1e-3 * distance.max(1e-3), "{ahead} {distance}");
        assert!(field.ahead(camera, -towards, 0.9).unwrap() > distance);
        // A camera that keeps flying at it with a speed set by that distance gets ever closer
        // and never past it: each step covers a part of what is left. (Closer than this, the
        // viewer goes by the exact places of the stars of the system, not by this search.)
        let mut place = camera;
        for _ in 0..800 {
            let left = field.ahead(place, towards, 0.9).unwrap();
            place += towards * left * 0.01;
        }
        let left = star.position.distance(place);
        assert!(left < 1e-3 * distance && (star.position - place).dot(towards) > 0.0, "{left} of {distance}");
    }

    #[test]
    fn an_address_gives_its_star_back() {
        let model = Model::default();
        let stars = generate(&model, 4, [3, 0, 40]);
        assert!(stars.len() > 3, "{}", stars.len());
        for (i, s) in stars.iter().enumerate() {
            assert_eq!(s.id, StarId { level: 4, cell: [3, 0, 40], index: i as u32 });
            assert_eq!(find(&model, s.id), Some(*s));
        }
        assert_eq!(find(&model, StarId { level: 4, cell: [3, 0, 40], index: stars.len() as u32 }), None);
        assert_ne!(stars[0].id.seed(&model), stars[1].id.seed(&model));
        // The smallest cubes have the faint stars too, down to the end of the table.
        let near = generate(&model, 0, [1, 0, 3500]);
        assert!(near.iter().any(|s| s.abs_mag > 12.0), "{} stars", near.len());
    }

    #[test]
    fn the_model_has_four_hundred_billion_stars() {
        let model = Model::default();
        // 270 billion systems, which with their companions are about 360 billion stars.
        let total = model.total();
        assert!((total / 2.7e11 - 1.0).abs() < 0.03, "{total:e} systems");
        let per_system = model.stars_per_system();
        assert!((1.25..1.45).contains(&per_system), "{per_system} stars a system");
        // At the Sun, between two arms of the model, the density is close to the real one
        // (about 0.1 stars per cubic parsec).
        let here = model.density(DVec3::ZERO);
        assert!((0.06..0.15).contains(&here), "{here} stars per cubic parsec");
        // Without the bulge and the arms the disc alone adds up as its formula says.
        let plain = Model { bulge_density: 0.0, arm_strength: 0.0, disc_edge: [1e9, 2e9], ..model };
        let expected = std::f64::consts::TAU * plain.local_density * 2.0 * plain.disc_height * (8000.0 / plain.disc_length).exp() * plain.disc_length.powi(2);
        let reach = Model { disc_edge: [1e9, 60_000.0], ..plain }.total();
        assert!((reach / expected - 1.0).abs() < 0.02, "{reach:e} against {expected:e}");
    }

    #[test]
    fn nothing_new_in_the_sky_of_the_earth() {
        let mut field = Field::default();
        while field.update(DVec3::ZERO, 10_000) {}
        let stars: Vec<(DVec3, &Kept)> = field.stars().collect();
        assert!(stars.len() > 1000, "{}", stars.len());
        for (place, s) in &stars {
            let m = s.abs_mag as f64 + 5.0 * (place.length() / 10.0).log10();
            assert!(m >= CATALOG_MAG - 1e-4, "a made-up star of magnitude {m} in the sky");
        }
    }

    #[test]
    fn stars_appear_far_from_the_sun_and_stay_few() {
        // Across the galaxy, at the place the far-side trip ends.
        let camera = DVec3::new(0.0, 2500.0, 16_000.0);
        let mut field = Field::default();
        while field.update(camera, 10_000) {}
        let visible = field.stars().filter(|(place, s)| (s.abs_mag as f64 + 5.0 * (place.distance(camera) / 10.0).log10()) < LIMIT_MAG).count();
        let total = field.stars().count();
        assert!(visible > 500, "{visible} stars to see");
        assert!(total < 2_000_000, "{total} stars in memory");
        // In the middle of the disc there are more.
        let camera = DVec3::new(0.0, 0.0, 13_000.0);
        let mut field = Field::default();
        while field.update(camera, 10_000) {}
        let dense = field.stars().filter(|(place, s)| (s.abs_mag as f64 + 5.0 * (place.distance(camera) / 10.0).log10()) < LIMIT_MAG).count();
        assert!(dense > visible, "{dense} in the disc, {visible} above it");
    }
}

#[cfg(test)]
mod light_tests {
    use super::*;

    #[test]
    fn the_haze_thins_as_fainter_stars_are_drawn() {
        let m = Model::default();
        let share = m.light_fainter_than(-12.0, 20.0, 33);
        assert!(share[0] > 0.999, "all the light fainter than -12: {}", share[0]);
        assert!(share[32] < 0.001, "none fainter than +20: {}", share[32]);
        assert!(share.windows(2).all(|w| w[1] <= w[0] + 1e-6), "never more for a fainter limit: {share:?}");
        // drawing the stars down to the Sun's brightness (+5) takes most of the light out of the haze
        assert!(share[17] < 0.5 && share[8] > share[17] + 0.3, "fainter than -4: {}, than +5: {}", share[8], share[17]);
    }
}
