//! What a made-up planet is like up close: its style (a handful of parameters: how rough the
//! ground is, how much of it is sea, its colors, its air), and the height of its ground at any
//! point, computed from the planet's own number.
//!
//! A style starts from a named one of the map, made after a planet of the solar system (an
//! Earth, a Mars, a Jupiter...), and is varied a little for each planet, so that the same few
//! styles serve every star. planet.wgsl draws the ground from these parameters; `height` here
//! computes the same ground on the CPU, so that the camera can land on what is drawn. The two
//! must be changed together.

use bytemuck::{Pod, Zeroable};
use glam::{DMat3, DVec3};

use crate::planets::{Kind, Planet};
use crate::procedural::Rng;

/// A style of planet: the parameters of the map, before they are varied for one planet.
#[derive(Clone, Debug, PartialEq)]
pub struct Preset {
    pub name: String,
    /// Highest ground over lowest, in km on a planet the size of the Earth.
    pub relief_km: f64,
    /// Share of the height range that is under the sea: 0 for a dry planet.
    pub sea: f64,
    /// How many times the main shapes of the ground repeat around the planet.
    pub roughness: f64,
    /// How strong the mountain ranges are, 0 to 1.
    pub mountains: f64,
    /// Where the ice caps begin: the sine of the latitude; over 1 for none.
    pub ice_cap: f64,
    /// Share of the sky covered by clouds.
    pub clouds: f64,
    /// Glow of molten rock in the low ground.
    pub glow: f64,
    /// Colors of the lowland, the upland and the peaks; for a planet of gas: of its dark bands,
    /// its light bands and its storm.
    pub low: [f32; 3],
    pub mid: [f32; 3],
    pub high: [f32; 3],
    pub sea_color: [f32; 3],
    pub air: [f32; 3],
    /// How deep the air is, in km on a planet the size of the Earth; 0 for none.
    pub air_km: f64,
    /// For a planet of gas: how many bands, and how sharply they differ.
    pub bands: f64,
    pub contrast: f64,
    /// With measured heights: how high the computed small shapes are that fill in between the
    /// measured points, in km (on the roughest ground; flat land gets a part of it).
    pub detail_km: f64,
    /// Plants: 0 for a bare planet, 1 for one where they grow wherever the climate lets them.
    pub life: f64,
    /// Colors of forest, of grassland, of sand and of bare rock.
    pub forest: [f32; 3],
    pub grass: [f32; 3],
    pub sand: [f32; 3],
    pub rock: [f32; 3],
    /// The climate: the mean temperature at sea level on the equator and at the poles, in
    /// degrees Celsius; how much colder it gets with each km of height; and how wet the planet
    /// is as a whole (1: like the Earth).
    pub temp_equator_c: f64,
    pub temp_pole_c: f64,
    /// How the cold is spread between equator and poles: 1 evenly, more keeps the warmth far
    /// from the equator and the cold near the poles.
    pub temp_power: f64,
    pub lapse_c_per_km: f64,
    pub rain: f64,
    /// How high the clouds float, in km on a planet the size of the Earth.
    pub cloud_km: f64,
    /// How much of the light of the star the air scatters on the way straight down, for the
    /// color it scatters most (the Earth's blue: 0.27).
    pub air_thickness: f64,
}

/// Measured heights of a planet: whole metres on a grid of longitude and latitude, as
/// tools/make_terrain.py writes them. Rows go from the north pole to the south, each row
/// eastwards from longitude -180.
///
/// The graphics card holds the whole grid for drawing. Here it is kept packed, in tiles that are
/// unpacked when the camera asks for the height under it: that takes a fifth of the memory, and
/// the camera only ever asks near where it is.
pub struct Terrain {
    pub name: String,
    pub width: usize,
    pub height: usize,
    /// The highest height of the grid, in metres.
    pub highest: f64,
    tiles: Vec<Vec<u8>>,
    unpacked: std::sync::Mutex<Vec<(usize, std::sync::Arc<Vec<i16>>)>>,
}

const TILE: usize = 512;

impl std::fmt::Debug for Terrain {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Terrain({}, {} x {})", self.name, self.width, self.height)
    }
}

impl Terrain {
    /// Reads `<path>.i16.gz` with its `<path>_metadata.txt`. Returns the terrain and the plain
    /// grid, which goes to the graphics card and is then let go.
    pub fn load(path: &std::path::Path) -> Result<(Terrain, Vec<i16>), String> {
        use std::io::Read;
        let name = path.file_name().map(|n| n.to_string_lossy().into_owned()).unwrap_or_default();
        let with = |suffix: &str| std::path::PathBuf::from(format!("{}{suffix}", path.display()));
        let meta = std::fs::read_to_string(with("_metadata.txt")).map_err(|e| format!("{}: {e}", with("_metadata.txt").display()))?;
        let size: Vec<usize> = meta.lines().find_map(|l| l.strip_prefix("Grid size:")).map(|l| l.split_whitespace().filter_map(|n| n.parse().ok()).collect()).unwrap_or_default();
        let [width, height] = size[..] else { return Err(format!("{}: no 'Grid size: W H'", with("_metadata.txt").display())) };
        let file = std::fs::File::open(with(".i16.gz")).map_err(|e| format!("{}: {e}", with(".i16.gz").display()))?;
        let mut bytes = Vec::with_capacity(width * height * 2);
        flate2::read::GzDecoder::new(file).read_to_end(&mut bytes).map_err(|e| format!("{}: {e}", with(".i16.gz").display()))?;
        if bytes.len() != width * height * 2 {
            return Err(format!("{}: {} bytes, expected {width} x {height} heights", with(".i16.gz").display(), bytes.len()));
        }
        let grid: Vec<i16> = bytes.chunks_exact(2).map(|b| i16::from_le_bytes([b[0], b[1]])).collect();
        drop(bytes);
        Ok((Terrain::from_grid(name, width, height, &grid), grid))
    }

    pub fn from_grid(name: String, width: usize, height: usize, grid: &[i16]) -> Terrain {
        use std::io::Write;
        let (across, down) = (width.div_ceil(TILE), height.div_ceil(TILE));
        let pack = |index: usize| {
            let (tx, ty) = (index % across, index / across);
            let mut raw = Vec::with_capacity(TILE * TILE * 2);
            for y in ty * TILE..((ty + 1) * TILE).min(height) {
                for x in tx * TILE..((tx + 1) * TILE).min(width) {
                    raw.extend_from_slice(&grid[y * width + x].to_le_bytes());
                }
            }
            let mut packer = flate2::write::DeflateEncoder::new(Vec::new(), flate2::Compression::fast());
            packer.write_all(&raw).expect("packing in memory");
            packer.finish().expect("packing in memory")
        };
        // The tiles are packed at the same time, a share to each thread.
        let count = across * down;
        #[cfg(not(target_family = "wasm"))]
        let threads = std::thread::available_parallelism().map_or(4, |n| n.get()).min(count.max(1));
        let mut tiles: Vec<Vec<u8>> = vec![Vec::new(); count];
        // on the web (no threads) one after another
        #[cfg(target_family = "wasm")]
        for (i, tile) in tiles.iter_mut().enumerate() {
            *tile = pack(i);
        }
        #[cfg(not(target_family = "wasm"))]
        std::thread::scope(|s| {
            for (k, share) in tiles.chunks_mut(count.div_ceil(threads).max(1)).enumerate() {
                let first = k * count.div_ceil(threads).max(1);
                let pack = &pack;
                s.spawn(move || {
                    for (i, tile) in share.iter_mut().enumerate() {
                        *tile = pack(first + i);
                    }
                });
            }
        });
        Terrain { name, width, height, highest: grid.iter().copied().max().unwrap_or(0) as f64, tiles, unpacked: std::sync::Mutex::new(Vec::new()) }
    }

    /// The height at a cell of the grid, in metres. Columns wrap around the planet; rows stop at
    /// the poles.
    pub fn at(&self, x: i64, y: i64) -> f64 {
        use std::io::Read;
        let x = x.rem_euclid(self.width as i64) as usize;
        let y = y.clamp(0, self.height as i64 - 1) as usize;
        let across = self.width.div_ceil(TILE);
        let index = (y / TILE) * across + x / TILE;
        let wide = TILE.min(self.width - (x / TILE) * TILE);
        let mut unpacked = self.unpacked.lock().expect("terrain tiles");
        let tile = match unpacked.iter().position(|(i, _)| *i == index) {
            Some(at) => {
                // The one just used goes to the front; the last is the first to go.
                let hit = unpacked.remove(at);
                unpacked.insert(0, hit);
                unpacked[0].1.clone()
            }
            None => {
                let mut raw = Vec::new();
                flate2::read::DeflateDecoder::new(&self.tiles[index][..]).read_to_end(&mut raw).expect("a tile packed here");
                let tile = std::sync::Arc::new(raw.chunks_exact(2).map(|b| i16::from_le_bytes([b[0], b[1]])).collect::<Vec<i16>>());
                unpacked.insert(0, (index, tile.clone()));
                unpacked.truncate(6);
                tile
            }
        };
        tile[(y % TILE) * wide + x % TILE] as f64
    }

    /// The height at a point of the planet (unit vector, the planet's axes: pole +Y, longitude 0
    /// along +Z, east towards +X), in metres. Between the points of the grid the height is
    /// blended smoothly, the same way planet.wgsl does it.
    pub fn metres(&self, n: DVec3) -> f64 {
        let lon = n.x.atan2(n.z);
        let lat = n.y.clamp(-1.0, 1.0).asin();
        let x = (lon / std::f64::consts::TAU + 0.5) * self.width as f64 - 0.5;
        let y = (0.5 - lat / std::f64::consts::PI) * self.height as f64 - 0.5;
        let (fx, fy) = (x.floor(), y.floor());
        let (tx, ty) = (x - fx, y - fy);
        let (ux, uy) = (tx * tx * (3.0 - 2.0 * tx), ty * ty * (3.0 - 2.0 * ty));
        let (ix, iy) = (fx as i64, fy as i64);
        let top = self.at(ix, iy) + (self.at(ix + 1, iy) - self.at(ix, iy)) * ux;
        let bottom = self.at(ix, iy + 1) + (self.at(ix + 1, iy + 1) - self.at(ix, iy + 1)) * ux;
        top + (bottom - top) * uy
    }
}

/// The styles of the default map, after the planets of the solar system.
pub fn default_presets() -> Vec<Preset> {
    let p = |name: &str, relief_km, sea, roughness, mountains, ice_cap, clouds, glow, low, mid, high, sea_color, air, air_km, bands, contrast| Preset {
        name: name.to_owned(),
        relief_km,
        sea,
        roughness,
        mountains,
        ice_cap,
        clouds,
        glow,
        low,
        mid,
        high,
        sea_color,
        air,
        air_km,
        bands,
        contrast,
        detail_km: 0.5,
        life: 0.0,
        forest: [0.17, 0.26, 0.09],
        grass: [0.50, 0.34, 0.26],
        sand: [0.70, 0.57, 0.41],
        rock: [0.42, 0.41, 0.50],
        temp_equator_c: 32.5,
        temp_pole_c: -14.2,
        temp_power: 3.27,
        lapse_c_per_km: 6.5,
        rain: 1.0,
        cloud_km: 5.0,
        air_thickness: if air_km > 0.0 { 0.27 } else { 0.0 },
    };
    let mut styles = vec![
        p("earth", 20.0, 0.52, 2.2, 0.9, 0.86, 0.5, 0.0, [0.16, 0.36, 0.12], [0.45, 0.38, 0.24], [0.95, 0.95, 0.97], [0.03, 0.16, 0.38], [0.35, 0.55, 0.95], 100.0, 0.0, 0.0),
        p("mars", 30.0, 0.0, 2.6, 0.7, 0.93, 0.03, 0.0, [0.62, 0.30, 0.16], [0.72, 0.42, 0.25], [0.86, 0.66, 0.50], [0.0, 0.0, 0.0], [0.85, 0.62, 0.45], 40.0, 0.0, 0.0),
        p("moon", 16.0, 0.0, 3.4, 0.4, 2.0, 0.0, 0.0, [0.30, 0.30, 0.31], [0.50, 0.49, 0.47], [0.72, 0.71, 0.69], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], 0.0, 0.0, 0.0),
        p("europa", 4.0, 0.0, 3.0, 0.25, 2.0, 0.0, 0.0, [0.74, 0.80, 0.88], [0.90, 0.93, 0.97], [0.99, 0.99, 1.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], 0.0, 0.0, 0.0),
        p("lava", 12.0, 0.0, 2.8, 0.6, 2.0, 0.0, 1.0, [0.10, 0.06, 0.05], [0.20, 0.13, 0.10], [0.32, 0.25, 0.22], [0.0, 0.0, 0.0], [0.6, 0.3, 0.2], 30.0, 0.0, 0.0),
        p("jupiter", 0.0, 0.0, 1.0, 0.0, 2.0, 0.0, 0.0, [0.55, 0.40, 0.28], [0.90, 0.82, 0.68], [0.72, 0.30, 0.18], [0.0, 0.0, 0.0], [0.85, 0.75, 0.6], 600.0, 17.0, 0.6),
        p("neptune", 0.0, 0.0, 1.0, 0.0, 2.0, 0.0, 0.0, [0.25, 0.42, 0.85], [0.40, 0.60, 0.95], [0.85, 0.92, 1.0], [0.0, 0.0, 0.0], [0.45, 0.65, 0.95], 400.0, 7.0, 0.9),
    ];
    for style in &mut styles {
        match style.name.as_str() {
            // The air of the Earth scatters blue the most: these are the shares of red, green
            // and blue that the air of the Earth scatters.
            "earth" => {
                style.life = 1.0;
                style.air = [0.175, 0.41, 1.0];
            }
            // Thin air full of dust.
            "mars" => style.air_thickness = 0.06,
            "lava" => style.air_thickness = 0.15,
            _ => {}
        }
    }
    styles
}

/// What one planet is like: a preset varied by the planet's own number.
#[derive(Clone, Debug, PartialEq)]
pub struct Style {
    pub preset: String,
    pub kind: Kind,
    pub seed: u32,
    /// Highest ground over lowest, in radii of the planet.
    pub relief: f64,
    pub sea: f64,
    pub roughness: f64,
    pub mountains: f64,
    pub ice_cap: f64,
    pub clouds: f64,
    pub glow: f64,
    pub low: [f32; 3],
    pub mid: [f32; 3],
    pub high: [f32; 3],
    pub sea_color: [f32; 3],
    pub air: [f32; 3],
    /// Depth of the air in radii of the planet.
    pub air_depth: f64,
    pub bands: f64,
    pub contrast: f64,
    /// The turn from the axes of space (galactic) to the planet's own, whose pole is +Y.
    pub to_planet: DMat3,
    /// A planet of gas has a station floating over its clouds: where it is from the centre of
    /// the planet (the planet's axes, radii) and how large (the radius of its wheel, in radii).
    pub station: Option<(DVec3, f64)>,
    /// Measured heights, if the planet has them: the ground is then theirs, with computed detail
    /// between the measured points; and how high that detail is, in metres.
    pub terrain: Option<std::sync::Arc<Terrain>>,
    pub detail_m: f64,
    /// Which picture of heights the shader reads: 0 for none, else one more than the place of
    /// the terrain in the renderer's list.
    pub terrain_slot: usize,
    /// The radius of the planet in metres.
    pub radius_m: f64,
    /// Plants, climate, clouds and air: as in the preset, with lengths in radii of the planet.
    pub life: f64,
    pub forest: [f32; 3],
    pub grass: [f32; 3],
    pub sand: [f32; 3],
    pub rock: [f32; 3],
    pub temp_equator_c: f64,
    pub temp_pole_c: f64,
    pub temp_power: f64,
    pub lapse_c_per_km: f64,
    pub rain: f64,
    pub cloud_height: f64,
    pub air_thickness: f64,
}

impl PartialEq for Terrain {
    fn eq(&self, other: &Terrain) -> bool {
        self.name == other.name && self.width == other.width && self.height == other.height
    }
}

const EARTH_RADIUS_KM: f64 = 6371.0;
/// The wheel of a station is this many km across from its middle.
const STATION_KM: f64 = 1.2;

impl Style {
    /// The style of a planet. `seed` is the planet's own number.
    pub fn of(planet: &Planet, seed: u64, presets: &[Preset]) -> Style {
        let mut rng = Rng::seeded(seed ^ 0x5354_594c_45);
        // A small rocky planet has lost its air, like the Moon; a larger one keeps some, like Mars.
        let name = match planet.kind {
            Kind::Temperate => "earth",
            Kind::Rocky if planet.mass < 0.15 => "moon",
            Kind::Rocky => "mars",
            Kind::Icy => "europa",
            Kind::Lava => "lava",
            Kind::GasGiant => "jupiter",
            Kind::IceGiant | Kind::SubNeptune => "neptune",
        };
        let fallback = default_presets();
        let preset = presets.iter().find(|p| p.name == name).or_else(|| fallback.iter().find(|p| p.name == name)).expect("a built-in style").clone();
        let radius_km = planet.radius * EARTH_RADIUS_KM;
        // Each planet differs a little from its style: tints, seas, heights.
        // (The draws are made up front, in a fixed order, so that a planet always gets the same.)
        let draws: Vec<f64> = (0..16).map(|_| rng.next()).collect();
        let mut next = draws.iter().copied();
        let mut vary = |v: f64, by: f64| v * (1.0 + by * (next.next().unwrap_or(0.5) * 2.0 - 1.0));
        let tint: [f32; 3] = [vary(1.0, 0.12) as f32, vary(1.0, 0.12) as f32, vary(1.0, 0.12) as f32];
        let tinted = |c: [f32; 3]| [0, 1, 2].map(|k| (c[k] * tint[k]).clamp(0.0, 1.0));
        // Small planets have higher mountains for their size: less gravity to pull them down.
        let relief = vary(preset.relief_km, 0.3) / radius_km * (1.0 / planet.radius.max(0.3)).sqrt();
        let sea = if preset.sea > 0.0 { (vary(preset.sea, 0.25)).clamp(0.2, 0.8) } else { 0.0 };
        // A pole tilted by up to 35 degrees, any way round.
        let (tilt, turn) = (draws[12] * 0.6, draws[13] * std::f64::consts::TAU);
        let to_planet = DMat3::from_rotation_x(tilt) * DMat3::from_rotation_y(turn);
        // Each planet has a climate of its own: warmer or colder, with more or less difference
        // between equator and poles, wetter or drier.
        let climate = ((draws[9] - 0.5) * 16.0, (draws[10] - 0.5) * 20.0, 0.6 + 0.8 * draws[11]);
        let gas = matches!(planet.kind, Kind::GasGiant | Kind::IceGiant | Kind::SubNeptune);
        let station = gas.then(|| {
            // Over the equator region, 60 km over the cloud tops.
            let (lon, lat) = (draws[14] * std::f64::consts::TAU, (draws[15] - 0.5) * 0.5);
            let dir = DVec3::new(lat.cos() * lon.cos(), lat.sin(), lat.cos() * lon.sin());
            (dir * (1.0 + 60.0 / radius_km), STATION_KM / radius_km)
        });
        Style {
            preset: preset.name.clone(),
            kind: planet.kind,
            seed: (seed ^ (seed >> 32)) as u32,
            relief,
            sea,
            roughness: vary(preset.roughness, 0.2),
            mountains: vary(preset.mountains, 0.3).clamp(0.0, 1.0),
            ice_cap: if preset.ice_cap > 1.0 { preset.ice_cap } else { vary(preset.ice_cap, 0.08).clamp(0.5, 0.99) },
            clouds: vary(preset.clouds, 0.4).clamp(0.0, 0.85),
            glow: preset.glow,
            low: tinted(preset.low),
            mid: tinted(preset.mid),
            high: tinted(preset.high),
            sea_color: tinted(preset.sea_color),
            air: preset.air,
            air_depth: preset.air_km / radius_km * (1.0 / planet.radius.max(0.3)).sqrt(),
            bands: vary(preset.bands, 0.3),
            contrast: preset.contrast,
            to_planet,
            station,
            terrain: None,
            detail_m: 0.0,
            terrain_slot: 0,
            radius_m: radius_km * 1000.0,
            life: preset.life,
            forest: tinted(preset.forest),
            grass: tinted(preset.grass),
            sand: tinted(preset.sand),
            rock: tinted(preset.rock),
            temp_equator_c: preset.temp_equator_c + climate.0,
            temp_pole_c: preset.temp_pole_c + climate.0 - climate.1,
            temp_power: preset.temp_power,
            lapse_c_per_km: preset.lapse_c_per_km,
            rain: preset.rain * climate.2,
            cloud_height: preset.cloud_km / radius_km * (1.0 / planet.radius.max(0.3)).sqrt(),
            air_thickness: preset.air_thickness,
        }
    }

    /// The style of a body that is known: a named style as it is, not varied, on a planet of
    /// the given radius; with measured heights if there are any.
    pub fn of_body(preset: &Preset, kind: Kind, radius_km: f64, seed: u32, terrain: Option<(std::sync::Arc<Terrain>, usize)>) -> Style {
        let size = radius_km / EARTH_RADIUS_KM;
        Style {
            preset: preset.name.clone(),
            kind,
            seed,
            relief: preset.relief_km / radius_km,
            sea: preset.sea,
            roughness: preset.roughness,
            mountains: preset.mountains,
            ice_cap: preset.ice_cap,
            clouds: preset.clouds,
            glow: preset.glow,
            low: preset.low,
            mid: preset.mid,
            high: preset.high,
            sea_color: preset.sea_color,
            air: preset.air,
            air_depth: preset.air_km / radius_km * (1.0 / size.max(0.3)).sqrt(),
            bands: preset.bands,
            contrast: preset.contrast,
            to_planet: DMat3::IDENTITY,
            station: None,
            detail_m: preset.detail_km * 1000.0,
            terrain_slot: terrain.as_ref().map_or(0, |t| t.1 + 1),
            terrain: terrain.map(|t| t.0),
            radius_m: radius_km * 1000.0,
            life: preset.life,
            forest: preset.forest,
            grass: preset.grass,
            sand: preset.sand,
            rock: preset.rock,
            temp_equator_c: preset.temp_equator_c,
            temp_pole_c: preset.temp_pole_c,
            temp_power: preset.temp_power,
            lapse_c_per_km: preset.lapse_c_per_km,
            rain: preset.rain,
            cloud_height: preset.cloud_km / radius_km * (1.0 / size.max(0.3)).sqrt(),
            air_thickness: preset.air_thickness,
        }
    }

    /// The dock of the station: where to stand on its deck and which way is ahead and up, from
    /// the centre of the planet (the planet's axes, radii). The station stands upright over the
    /// planet: its own axes are east, up and north there, as in planet.wgsl.
    pub fn dock(&self) -> Option<(DVec3, DVec3, DVec3)> {
        let (at, size) = self.station?;
        let up = at.normalize();
        let east = DVec3::Y.cross(up).normalize();
        // On the deck under the hub, near its edge, a little over the floor; looking out along it.
        Some((at + (east * 0.3 + up * (-0.34 + 0.03 + 0.012)) * size, east, up))
    }

    pub fn is_gas(&self) -> bool {
        matches!(self.kind, Kind::GasGiant | Kind::IceGiant | Kind::SubNeptune)
    }

    /// How high the ground is at a point of the planet (unit vector, the planet's axes): 0 the
    /// lowest, about 1 the highest. The same as `height` in planet.wgsl with all its detail.
    pub fn height(&self, n: DVec3) -> f64 {
        if let Some(terrain) = &self.terrain {
            let metres = terrain.metres(n);
            let scale = 1.0 / (self.relief * self.radius_m);
            return self.sea + metres * scale + measured_detail(n, metres, self.seed) * self.detail_m * scale;
        }
        let p = n * self.roughness;
        let continents = rough(p, 5, self.seed);
        let ridged = 1.0 - (2.0 * rough(p * 4.0 + DVec3::splat(17.0), 4, self.seed.wrapping_add(100)) - 1.0).abs();
        let inland = smoothstep(self.sea, self.sea + 0.25, continents);
        continents + self.mountains * ridged * ridged * inland * 0.35 + (rougher(p * 24.0, 10, self.seed.wrapping_add(200), 0.62) - 0.5) * 0.14 * (0.25 + inland)
    }

    /// The height of the surface over the mean radius at a point, in radii: the ground, or the
    /// sea over it; for a planet of gas, its cloud tops.
    pub fn surface(&self, n: DVec3) -> f64 {
        if self.is_gas() { 0.0 } else { self.relief * (self.height(n) - self.sea).max(0.0) }
    }

    /// The numbers of the planet for the shader. `camera` is the camera from the centre of the
    /// planet (galactic axes, radii), `light` the direction of its star, `brightness` how
    /// strongly the star shines on it.
    /// `hours` is the time, which moves the clouds; `star_radius` how large the star looks from
    /// the planet, as an angle in radians.
    pub fn uniform(&self, camera: DVec3, light: DVec3, brightness: f64, hours: f64, star_radius: f64) -> PlanetUniform {
        let v3 = |v: DVec3, w: f64| [v.x as f32, v.y as f32, v.z as f32, w as f32];
        let rgb = |c: [f32; 3], w: f64| [c[0], c[1], c[2], w as f32];
        let rows = self.to_planet.transpose();
        let (station, size) = self.station.unwrap_or((DVec3::ZERO, 0.0));
        PlanetUniform {
            cam: v3(camera, camera.length() - 1.0),
            light: v3(light, brightness),
            to_planet: [v3(rows.x_axis, 0.0), v3(rows.y_axis, 0.0), v3(rows.z_axis, 0.0)],
            shape: [self.relief as f32, self.sea as f32, self.roughness as f32, self.mountains as f32],
            look: [self.kind.code() as f32, self.ice_cap as f32, self.clouds as f32, self.glow as f32],
            low: rgb(self.low, self.bands),
            mid: rgb(self.mid, self.contrast),
            high: rgb(self.high, 0.0),
            sea: rgb(self.sea_color, 0.0),
            air: rgb(self.air, self.air_depth),
            station: v3(station, size),
            terrain: [self.terrain.is_some() as u32 as f32, (1.0 / (self.relief * self.radius_m)) as f32, self.detail_m as f32, self.terrain.as_ref().map_or(0.0, |t| (t.highest.max(0.0) + self.detail_m) / self.radius_m) as f32],
            climate: [self.temp_equator_c as f32, self.temp_pole_c as f32, self.lapse_c_per_km as f32, self.rain as f32],
            forest: rgb(self.forest, self.life),
            grass: rgb(self.grass, 0.0),
            sand: rgb(self.sand, self.air_thickness),
            rock: rgb(self.rock, 0.0),
            misc: [self.cloud_height as f32, (hours % 100_000.0) as f32, star_radius as f32, self.temp_power as f32],
            seed: [self.seed, self.seed.rotate_left(13) ^ 0x9e37_79b9, 0, 0],
        }
    }
}

/// A planet as planet.wgsl reads it.
#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable)]
pub struct PlanetUniform {
    cam: [f32; 4],
    light: [f32; 4],
    to_planet: [[f32; 4]; 3],
    shape: [f32; 4],
    look: [f32; 4],
    low: [f32; 4],
    mid: [f32; 4],
    high: [f32; 4],
    sea: [f32; 4],
    air: [f32; 4],
    station: [f32; 4],
    terrain: [f32; 4],
    climate: [f32; 4],
    forest: [f32; 4],
    grass: [f32; 4],
    sand: [f32; 4],
    rock: [f32; 4],
    misc: [f32; 4],
    seed: [u32; 4],
}

/// The computed small shapes between the points of measured heights: -0.5 to 0.5, to be
/// multiplied by how high they are. Mountains are rough, plains nearly flat, the sea floor is
/// left as measured. The same as in planet.wgsl with all its detail.
fn measured_detail(n: DVec3, metres: f64, seed: u32) -> f64 {
    let rough_here = 0.12 + 0.88 * smoothstep(100.0, 2500.0, metres);
    (rougher(n * 1600.0, 10, seed.wrapping_add(200), 0.62) - 0.5) * rough_here * smoothstep(-30.0, 30.0, metres)
}

fn smoothstep(a: f64, b: f64, v: f64) -> f64 {
    let t = ((v - a) / (b - a)).clamp(0.0, 1.0);
    t * t * (3.0 - 2.0 * t)
}

/// Whole numbers in, a number from 0 to 1 out: the same in the shader, bit for bit.
fn lattice(i: [i32; 3], seed: u32) -> f64 {
    let mut h = (i[0] as u32).wrapping_mul(0x8da6_b343) ^ (i[1] as u32).wrapping_mul(0xd816_3841) ^ (i[2] as u32).wrapping_mul(0xcb1a_b31f) ^ seed.wrapping_mul(0x9e37_79b9);
    h ^= h >> 15;
    h = h.wrapping_mul(0x2c1b_3c6d);
    h ^= h >> 12;
    h = h.wrapping_mul(0x297a_2d39);
    h ^= h >> 15;
    (h >> 8) as f64 / 16_777_216.0
}

fn noise(p: DVec3, seed: u32) -> f64 {
    let f = p.floor();
    let i = [f.x as i32, f.y as i32, f.z as i32];
    let t = p - f;
    let u = t * t * (DVec3::splat(3.0) - 2.0 * t);
    let at = |dx: i32, dy: i32, dz: i32| lattice([i[0] + dx, i[1] + dy, i[2] + dz], seed);
    let mix = |a: f64, b: f64, t: f64| a + (b - a) * t;
    mix(
        mix(mix(at(0, 0, 0), at(1, 0, 0), u.x), mix(at(0, 1, 0), at(1, 1, 0), u.x), u.y),
        mix(mix(at(0, 0, 1), at(1, 0, 1), u.x), mix(at(0, 1, 1), at(1, 1, 1), u.x), u.y),
        u.z,
    )
}

/// Noise with detail on `octaves` scales, each twice as fine and `gain` times as strong; 0 to 1.
fn rougher(p: DVec3, octaves: u32, seed: u32, gain: f64) -> f64 {
    let (mut sum, mut size, mut total, mut q) = (0.0, 0.5, 0.0, p);
    for k in 0..octaves {
        sum += noise(q, seed.wrapping_add(k)) * size;
        total += size;
        q *= 2.0;
        size *= gain;
    }
    sum / total
}

fn rough(p: DVec3, octaves: u32, seed: u32) -> f64 {
    rougher(p, octaves, seed, 0.5)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::planets::Host;
    use crate::systems::orbit;

    fn planet(kind: Kind, radius: f64, mass: f64) -> Planet {
        let mut rng = Rng::seeded(3);
        Planet { kind, host: Host::Star(0), mass, radius, orbit: orbit(1.0, 0.0, 1.0, &mut rng), temperature: 280.0, moons: 0, rings: false }
    }

    #[test]
    fn styles_follow_the_kind_and_differ_between_planets() {
        let presets = default_presets();
        let earth = Style::of(&planet(Kind::Temperate, 1.0, 1.0), 11, &presets);
        assert_eq!(earth.preset, "earth");
        assert!(earth.sea > 0.2 && earth.air_depth > 0.0 && earth.station.is_none());
        // About 20 km from the deepest sea floor to the highest peak, on a planet like the Earth.
        assert!((earth.relief * EARTH_RADIUS_KM - 20.0).abs() < 7.0, "{} km", earth.relief * EARTH_RADIUS_KM);
        assert_eq!(Style::of(&planet(Kind::Rocky, 0.4, 0.05), 11, &presets).preset, "moon");
        assert_eq!(Style::of(&planet(Kind::Rocky, 0.9, 0.7), 11, &presets).preset, "mars");
        let jupiter = Style::of(&planet(Kind::GasGiant, 11.0, 318.0), 11, &presets);
        assert!(jupiter.is_gas() && jupiter.station.is_some() && jupiter.bands > 5.0);
        // The same planet, the same style; another planet, another.
        assert_eq!(earth, Style::of(&planet(Kind::Temperate, 1.0, 1.0), 11, &presets));
        assert_ne!(earth, Style::of(&planet(Kind::Temperate, 1.0, 1.0), 12, &presets));
    }

    #[test]
    fn the_ground_is_continuous_and_has_land_and_sea() {
        let earth = Style::of(&planet(Kind::Temperate, 1.0, 1.0), 5, &default_presets());
        let (mut land, mut total, mut highest) = (0, 0, 0.0f64);
        let mut rng = Rng::seeded(9);
        for _ in 0..4000 {
            let n = DVec3::new(rng.next() - 0.5, rng.next() - 0.5, rng.next() - 0.5).normalize();
            let s = earth.surface(n);
            land += (s > 0.0) as u32;
            total += 1;
            highest = highest.max(s);
            // A step of a hundred metres changes the height by far less than the relief.
            let near = (n + DVec3::new(1.6e-5, 0.0, 0.0)).normalize();
            assert!((earth.surface(near) - s).abs() < 0.1 * earth.relief, "a cliff at {n:?}");
        }
        let share = land as f64 / total as f64;
        assert!((0.15..0.85).contains(&share), "{share} of the planet is land");
        assert!(highest > 0.1 * earth.relief && highest < 1.3 * earth.relief);
        // A planet of gas has no ground: its surface is its cloud tops.
        assert_eq!(Style::of(&planet(Kind::GasGiant, 11.0, 318.0), 5, &default_presets()).surface(DVec3::X), 0.0);
    }
}
