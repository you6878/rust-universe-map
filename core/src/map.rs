//! The map of the universe: a JSON file that lists every layer the viewer draws and the files it
//! is made of. See `map/universe.json`.

use std::path::{Path, PathBuf};

use glam::DVec3;
use serde_json::Value;

use crate::procedural::Model;

/// The map that ships with the viewer, for when no file is given.
const DEFAULT: &str = include_str!("../../map/universe.json");

/// How many volumes the shader marches through at once.
pub const MAX_VOLUMES: usize = 6;

/// Gas, starlight and dust as grids of voxels.
pub struct VolumeSpec {
    pub name: String,
    /// The grid of light and the grid of dust. A volume may have only one of them.
    pub light: Option<PathBuf>,
    pub dust: Option<PathBuf>,
    pub metadata: PathBuf,
    /// Length in parsecs of one unit of the box in the metadata.
    pub size_pc: f64,
    /// Middle of the box, galactic axes, parsecs.
    pub centre_pc: DVec3,
    /// Where the three axes of the box point, in galactic axes. Without them in the map, the box
    /// lies along the galactic axes: to the centre, north, to longitude 90.
    pub axes: [DVec3; 3],
    /// What a kiloparsec of full density does: the brightness it adds, the optical depth of its
    /// dust, and the brightness its dust gives back.
    pub light_per_kpc: f32,
    pub dust_per_kpc: f32,
    pub dust_glow_per_kpc: f32,
    /// The density is the stored value raised to this power. More than 1 gives the eight bits of
    /// a voxel more range for thin material.
    pub gamma: f32,
    /// Distances of the camera from the centre between which the volume fades away, in parsecs.
    /// For what only holds as seen from near its centre.
    pub fade_out_pc: Option<(f64, f64)>,
    /// Fine clumps for the dust, finer than its voxels: the noise of the map repeats every
    /// `detail_pc` parsecs, and makes the dust uneven by `detail_strength` (0 for none).
    pub detail_pc: f32,
    pub detail_strength: f32,
    /// Whether the white balance of the map's exposure is for the light of this volume. It is
    /// made for the galaxy the reference place is in; another galaxy keeps its own colors.
    pub white_balance: bool,
}

impl VolumeSpec {
    /// How much of the volume is there for a camera at a place (galactic parsecs): 1 to 0.
    pub fn presence(&self, camera_pc: DVec3) -> f32 {
        let Some((from, to)) = self.fade_out_pc else { return 1.0 };
        let t = ((camera_pc.distance(self.centre_pc) - from) / (to - from).max(1e-9)).clamp(0.0, 1.0);
        (1.0 - t * t * (3.0 - 2.0 * t)) as f32
    }
}

/// The exposure of the volumes: how the camera adapts to how bright they are around it.
pub struct Exposure {
    /// Where the exposure is 1, in galactic parsecs.
    pub at_pc: DVec3,
    pub min: f32,
    pub max: f32,
    /// What red, green and blue are multiplied by at that place, and the distances from it
    /// between which this white balance goes back to none.
    pub white_balance: [f32; 3],
    pub white_balance_fade_pc: (f64, f64),
    /// A narrow view gathers light like a telescope: below this field of view (degrees, top to
    /// bottom) the exposure grows as the view narrows, up to `zoom_max` times.
    pub zoom_from_fov_deg: f64,
    pub zoom_max: f32,
}

impl Exposure {
    /// How much more light a view of a field (radians, top to bottom) gathers.
    pub fn zoom(&self, fov_y: f64) -> f32 {
        ((self.zoom_from_fov_deg.to_radians() / fov_y.max(1e-6)) as f32).clamp(1.0, self.zoom_max.max(1.0))
    }

    /// The white balance for a camera at a place, in galactic parsecs.
    pub fn balance(&self, camera_pc: DVec3) -> [f32; 3] {
        let (from, to) = self.white_balance_fade_pc;
        let t = ((camera_pc.distance(self.at_pc) - from) / (to - from).max(1e-9)).clamp(0.0, 1.0);
        let there = (1.0 - t * t * (3.0 - 2.0 * t)) as f32;
        self.white_balance.map(|w| 1.0 + (w - 1.0) * there)
    }
}

pub struct Map {
    /// The data folder, which `$data` in the paths of the catalogs stands for.
    pub data: PathBuf,
    pub star_catalogs: Vec<PathBuf>,
    pub star_field: Option<Model>,
    /// Every galaxy of starField.fitted, those with a grid of their own too: what they look like
    /// from afar (fitted.rs) is drawn from these.
    pub galaxy_looks: Vec<crate::fitted::ParamGalaxy>,
    pub volumes: Vec<VolumeSpec>,
    /// A cube of noise that repeats, one byte a cell, for the fine clumps of the dust.
    pub detail_noise: Option<PathBuf>,
    pub exposure: Option<Exposure>,
    pub planets: PathBuf,
    pub moons: PathBuf,
    pub vsop87: PathBuf,
    pub earth_day: PathBuf,
    pub earth_night: PathBuf,
    pub earth_clouds: PathBuf,
    /// Bodies of the Solar System that are drawn by the planet object: the name of the body, the
    /// style it has, and its measured heights if there are any (`<path>.i16.gz`).
    pub body_planets: Vec<(String, String, Option<PathBuf>)>,
    /// The styles of planets: named sets of parameters, after planets of the solar system.
    pub planet_styles: Vec<crate::ground::Preset>,
    /// Kinds of rare objects: name, and how many in all or how many per star.
    pub rare: Vec<(String, Option<f64>, Option<f64>)>,
    pub picture_catalogs: Vec<PathBuf>,
    /// Pictures with a longer side than this many pixels are made smaller when read.
    pub picture_max_side: Option<u32>,
    /// The folder of the ship that is in the scene (`world/ships/<name>`), if there is one.
    pub ship: Option<PathBuf>,
}

/// Turns a path of the map into a real one.
fn resolve(spec: &str, data: &Path, base: &Path) -> PathBuf {
    let path = if let Some(rest) = spec.strip_prefix("$data/") {
        data.join(rest)
    } else if let Some(rest) = spec.strip_prefix("~/") {
        PathBuf::from(std::env::var("HOME").unwrap_or_default()).join(rest)
    } else {
        base.join(spec)
    };
    if !spec.contains('*') {
        return path;
    }
    // The star is a quality suffix. The best one installed wins.
    let text = path.to_string_lossy().into_owned();
    ["-high", "-med", "-low", ""].iter().map(|q| PathBuf::from(text.replace('*', q))).find(|p| p.exists()).unwrap_or(path)
}

struct Reader<'a> {
    data: PathBuf,
    base: &'a Path,
}

impl Reader<'_> {
    fn path(&self, v: &Value, key: &str, context: &str) -> Result<PathBuf, String> {
        let spec = v.get(key).and_then(Value::as_str).ok_or(format!("the map has no \"{key}\" in {context}"))?;
        let path = resolve(spec, &self.data, self.base);
        if !path.exists() {
            return Err(format!("{context}.{key}: {} does not exist", path.display()));
        }
        Ok(path)
    }

    fn paths(&self, v: Option<&Value>, context: &str) -> Result<Vec<PathBuf>, String> {
        let Some(list) = v.and_then(|v| v.get("catalogs")).and_then(Value::as_array) else { return Ok(Vec::new()) };
        list.iter()
            .map(|s| {
                let path = resolve(s.as_str().ok_or(format!("{context}.catalogs must be paths"))?, &self.data, self.base);
                if path.exists() { Ok(path) } else { Err(format!("{context}.catalogs: {} does not exist", path.display())) }
            })
            .collect()
    }
}

fn number(v: &Value, key: &str, default: f64) -> f64 {
    v.get(key).and_then(Value::as_f64).unwrap_or(default)
}

fn vector(v: &Value, key: &str, default: DVec3) -> DVec3 {
    match v.get(key).and_then(Value::as_array) {
        Some(a) if a.len() == 3 => DVec3::new(a[0].as_f64().unwrap_or(0.0), a[1].as_f64().unwrap_or(0.0), a[2].as_f64().unwrap_or(0.0)),
        _ => default,
    }
}

fn star_field(v: &Value) -> Model {
    let d = Model::default();
    let pair = |key: &str, default: [f64; 2]| match v.get(key).and_then(Value::as_array) {
        Some(a) if a.len() == 2 => [a[0].as_f64().unwrap_or(default[0]), a[1].as_f64().unwrap_or(default[1])],
        _ => default,
    };
    let empty = Value::Null;
    let arms = v.get("arms").unwrap_or(&empty);
    let bulge = v.get("bulge").unwrap_or(&empty);
    let planets = v.get("planets").unwrap_or(&empty);
    Model {
        seed: v.get("seed").and_then(Value::as_u64).unwrap_or(d.seed),
        centre: vector(v, "centrePc", d.centre),
        local_density: number(v, "localDensity", d.local_density),
        disc_length: number(v, "discScaleLengthPc", d.disc_length),
        disc_height: number(v, "discScaleHeightPc", d.disc_height),
        disc_edge: pair("discEdgePc", d.disc_edge),
        arm_count: number(arms, "count", d.arm_count),
        arm_strength: number(arms, "strength", d.arm_strength),
        arm_pitch: number(arms, "pitch", d.arm_pitch),
        arm_start: number(arms, "startPc", d.arm_start),
        bulge_density: number(bulge, "density", d.bulge_density),
        bulge_radius: number(bulge, "radiusPc", d.bulge_radius),
        bulge_flattening: number(bulge, "flattening", d.bulge_flattening),
        // The table of the map: [absolute magnitude, share of the stars brighter] pairs.
        brightness: match v.get("brightness").and_then(Value::as_array) {
            Some(rows) => {
                let table: Vec<(f64, f64)> = rows.iter().filter_map(|r| Some((r.get(0)?.as_f64()?, r.get(1)?.as_f64()?))).collect();
                let rising = table.windows(2).all(|w| w[0].0 < w[1].0 && w[0].1 < w[1].1);
                if table.len() >= 2 && rising && table[0].1 > 0.0 { table } else { d.brightness }
            }
            None => d.brightness,
        },
        white_dwarfs: number(v, "whiteDwarfShare", d.white_dwarfs),
        planets: crate::planets::Abundance {
            mean_planets: number(planets, "meanPerStar", d.planets.mean_planets),
            giants_at_one_solar_mass: number(planets, "giantsAtOneSolarMass", d.planets.giants_at_one_solar_mass),
            hot_share: number(planets, "hotShareOfGiants", d.planets.hot_share),
        },
        // Rows of [least mass of the first star, share of doubles, of triples, of more].
        companions: match v.get("companions").and_then(Value::as_array) {
            Some(rows) => {
                let table: Vec<_> = rows
                    .iter()
                    .filter_map(|r| Some(crate::systems::Multiplicity { min_mass: r.get(0)?.as_f64()?, double: r.get(1)?.as_f64()?, triple: r.get(2)?.as_f64()?, more: r.get(3)?.as_f64()? }))
                    .collect();
                if table.is_empty() || table.len() != rows.len() { d.companions } else { table }
            }
            None => d.companions,
        },
        galaxies: Vec::new(),
        fitted: Vec::new(),
        home: None,
    }
}

/// starField.galaxies: [{"name", "volume": the name of a volume with a light grid, "systems"}], or
/// with "photo" (a frame file of tools/galaxy_xyz.py) and "look" (the name of its galaxy in
/// starField.fitted, which says where it is and how it is turned) instead of "volume": its stars
/// spread by three pictures (photo_volume.rs).
/// Their stars are spread as the light of the volume, but for "bulgeShare" of them, which are in a
/// round bulge of "bulgeRadiusPc". "plane" is a finer face-on map of the disc (tools/photo_plane.py) that spreads the stars within the grid's cells as the photograph does. "absMag" is how bright the whole galaxy is, for its glow.
fn grid_galaxies(field: &Value, volumes: &[VolumeSpec], looks: &[crate::fitted::ParamGalaxy], r: &Reader) -> Result<Vec<crate::procedural::GridGalaxy>, String> {
    let mut list = Vec::new();
    for (i, g) in field.get("galaxies").and_then(Value::as_array).into_iter().flatten().enumerate() {
        let context = format!("starField.galaxies[{i}]");
        if let Some(photo) = g.get("photo").and_then(Value::as_str) {
            // its pictures are read when the camera comes near (procedural.rs Model::keep_near)
            let volume = crate::photo_volume::LazyVolume::new(&resolve(photo, &r.data, r.base)).map_err(|e| format!("{context}: {e}"))?;
            let wanted = g.get("look").and_then(Value::as_str).ok_or(format!("{context} has a \"photo\" but no \"look\""))?;
            let look = looks.iter().find(|l| l.name == wanted).ok_or(format!("{context}: no galaxy of starField.fitted is called \"{wanted}\""))?;
            let (along, up) = volume.reach;
            list.push(crate::procedural::GridGalaxy {
                name: g.get("name").and_then(Value::as_str).unwrap_or(wanted).to_owned(),
                systems: number(g, "systems", 0.0),
                centre: look.centre,
                to_box: look.to_box,
                lo: DVec3::new(-along, -up, -along),
                hi: DVec3::new(along, up, along),
                size: [1, 1, 1],
                density: std::sync::Arc::new(vec![0.0]),
                bulge_density: 0.0,
                bulge_radius: 1.0,
                colour: std::sync::Arc::new(vec![[0; 3]]),
                abs_mag: number(g, "absMag", -21.0),
                plane: std::sync::Arc::new(Default::default()),
                photo: Some(std::sync::Arc::new(volume)),
                light_only: g.get("lightOnly").and_then(Value::as_bool).unwrap_or(false),
                dust: match g.get("dust").and_then(Value::as_str) {
                    Some(d) if volumes.iter().any(|v| v.name == d && v.dust.is_some()) => Some(d.to_owned()),
                    Some(d) => return Err(format!("{context}: no volume called \"{d}\" has a dust grid")),
                    None => None,
                },
            });
            continue;
        }
        let wanted = g.get("volume").and_then(Value::as_str).ok_or(format!("{context} has no \"volume\""))?;
        let vol = volumes.iter().find(|v| v.name == wanted).ok_or(format!("{context}: no volume is called \"{wanted}\""))?;
        let light = vol.light.as_ref().ok_or(format!("{context}: the volume \"{wanted}\" has no light grid"))?;
        let systems = number(g, "systems", 0.0);
        let name = g.get("name").and_then(Value::as_str).unwrap_or(wanted);
        let plane = g.get("plane").and_then(Value::as_str).map(|p| resolve(p, &r.data, r.base));
        let galaxy = crate::procedural::GridGalaxy::load(name, light, &vol.metadata, vol.size_pc, vol.centre_pc, vol.axes, vol.gamma as f64, systems, (number(g, "bulgeShare", 0.0), number(g, "bulgeRadiusPc", 1000.0)), plane.as_deref()).map_err(|e| format!("{context}: {e}"))?;
        list.push(crate::procedural::GridGalaxy { abs_mag: number(g, "absMag", -21.0), light_only: false, ..galaxy });
    }
    Ok(list)
}

impl Map {
    /// Reads a map. Without a file, the map that ships with the viewer is used. `data` replaces
    /// the data folder named in the map.
    pub fn load(file: Option<&Path>, data: Option<&Path>) -> Result<Map, String> {
        let (text, base, name) = match file {
            Some(path) => {
                let text = std::fs::read_to_string(path).map_err(|e| format!("{}: {e}", path.display()))?;
                (text, path.parent().unwrap_or(Path::new(".")).to_path_buf(), path.display().to_string())
            }
            // The built-in map is the file of the source tree, and its paths are relative to it.
            None => (DEFAULT.to_owned(), Path::new(env!("CARGO_MANIFEST_DIR")).join("../map"), "the built-in map".to_owned()),
        };
        Map::parse(&text, &base, data).map_err(|e| format!("{name}: {e}"))
    }

    fn parse(text: &str, base: &Path, data: Option<&Path>) -> Result<Map, String> {
        let v: Value = serde_json::from_str(text).map_err(|e| e.to_string())?;
        let data = match data {
            Some(d) => d.to_path_buf(),
            None => resolve(v.get("data").and_then(Value::as_str).ok_or("the map has no \"data\" folder")?, Path::new(""), base),
        };
        let r = Reader { data, base };
        // universe-map: a map without "bodies" has only what is computed (no Solar System, which
        // needs the data of Gaia Sky); its paths are then empty.
        let none = Value::Null;
        let has_bodies = v.get("bodies").is_some();
        let bodies = v.get("bodies").unwrap_or(&none);
        let earth = if has_bodies { bodies.get("earth").ok_or("the map has no \"bodies.earth\"")? } else { &none };
        let mut volumes = Vec::new();
        for (i, vol) in v.get("volumes").and_then(Value::as_array).into_iter().flatten().enumerate() {
            let context = format!("volumes[{i}]");
            volumes.push(VolumeSpec {
                name: vol.get("name").and_then(Value::as_str).unwrap_or("volume").to_owned(),
                light: vol.get("light").map(|_| r.path(vol, "light", &context)).transpose()?,
                dust: vol.get("dust").map(|_| r.path(vol, "dust", &context)).transpose()?,
                metadata: r.path(vol, "metadata", &context)?,
                size_pc: number(vol, "sizePc", 1.0),
                centre_pc: vector(vol, "centrePc", DVec3::ZERO),
                axes: match vol.get("axes").and_then(Value::as_array) {
                    Some(rows) if rows.len() == 3 => {
                        let row = |k: usize| rows[k].as_array().filter(|r| r.len() == 3).map(|r| DVec3::new(r[0].as_f64().unwrap_or(0.0), r[1].as_f64().unwrap_or(0.0), r[2].as_f64().unwrap_or(0.0)).normalize_or_zero());
                        match (row(0), row(1), row(2)) {
                            (Some(a), Some(b), Some(c)) if a != DVec3::ZERO && b != DVec3::ZERO && c != DVec3::ZERO => [a, b, c],
                            _ => return Err(format!("{context}.axes must be three vectors of three numbers")),
                        }
                    }
                    Some(_) => return Err(format!("{context}.axes must be three vectors of three numbers")),
                    None => [DVec3::Z, DVec3::Y, DVec3::X],
                },
                light_per_kpc: number(vol, "lightPerKpc", 0.0) as f32,
                dust_per_kpc: number(vol, "dustPerKpc", 0.0) as f32,
                dust_glow_per_kpc: number(vol, "dustGlowPerKpc", 0.0) as f32,
                gamma: number(vol, "gamma", 1.0) as f32,
                white_balance: vol.get("whiteBalance").and_then(Value::as_bool).unwrap_or(true),
                detail_pc: number(vol, "detailPc", 400.0) as f32,
                detail_strength: number(vol, "detailStrength", 0.0) as f32,
                fade_out_pc: vol.get("fadeOutPc").and_then(Value::as_array).and_then(|a| Some((a.first()?.as_f64()?, a.get(1)?.as_f64()?))),
            });
        }
        if volumes.len() > MAX_VOLUMES {
            return Err(format!("the map has {} volumes, and at most {MAX_VOLUMES} can be drawn", volumes.len()));
        }
        if let Some(v) = volumes.iter().find(|v| v.light.is_none() && v.dust.is_none()) {
            return Err(format!("the volume \"{}\" has neither \"light\" nor \"dust\"", v.name));
        }
        // starField.fitted: one file of galaxies (tools/galaxy_views.py, galaxy_imagine.py) or a list
        let mut looks = Vec::new();
        let files: Vec<&str> = match v.get("starField").and_then(|f| f.get("fitted")) {
            Some(Value::String(f)) => vec![f.as_str()],
            Some(Value::Array(list)) => list.iter().filter_map(Value::as_str).collect(),
            _ => Vec::new(),
        };
        for file in files {
            looks.extend(crate::fitted::load(&resolve(file, &r.data, r.base), &r.data)?);
        }
        Ok(Map {
            star_catalogs: r.paths(v.get("stars"), "stars")?,
            star_field: match v.get("starField") {
                Some(f) => {
                    let galaxies = grid_galaxies(f, &volumes, &looks, &r)?;
                    // the galaxies as formulas, but those that have a grid of their own
                    let fitted = looks.iter().filter(|p| !galaxies.iter().any(|g| p.names.iter().any(|n| n.to_lowercase().contains(&g.name.to_lowercase())))).cloned().collect();
                    let model = star_field(f);
                    // the Milky Way from three pictures, where the map gives them
                    let home = match f.get("photo").and_then(Value::as_str) {
                        Some(file) => {
                            let path = resolve(file, &r.data, r.base);
                            let volume = crate::photo_volume::PhotoVolume::load(&path).map_err(|e| format!("starField.photo: {e}"))?;
                            let thick = f.get("thickDisc").unwrap_or(&Value::Null);
                            Some(crate::procedural::HomeVolume::new(std::sync::Arc::new(volume), model.centre, model.local_density, number(thick, "share", 0.0), number(thick, "heightPc", 900.0)))
                        }
                        None => None,
                    };
                    Some(Model { galaxies, fitted, home, ..model })
                }
                None => None,
            },
            galaxy_looks: looks,
            volumes,
            detail_noise: v.get("detail").map(|d| r.path(d, "noise", "detail")).transpose()?,
            exposure: v.get("exposure").map(|e| {
                let list = |key: &str| e.get(key).and_then(Value::as_array).map(|a| a.iter().filter_map(Value::as_f64).collect::<Vec<_>>());
                let balance = list("whiteBalance").filter(|b| b.len() == 3).map_or([1.0; 3], |b| [b[0] as f32, b[1] as f32, b[2] as f32]);
                let fade = list("whiteBalanceFadePc").filter(|f| f.len() == 2).map_or((f64::MAX / 4.0, f64::MAX / 2.0), |f| (f[0], f[1]));
                Exposure {
                    at_pc: vector(e, "atPc", DVec3::ZERO),
                    min: number(e, "min", 1.0) as f32,
                    max: number(e, "max", 1.0) as f32,
                    white_balance: balance,
                    white_balance_fade_pc: fade,
                    zoom_from_fov_deg: number(e, "zoomFromFovDeg", 0.0),
                    zoom_max: number(e, "zoomMax", 1.0) as f32,
                }
            }),
            planets: if has_bodies { r.path(bodies, "planets", "bodies")? } else { PathBuf::new() },
            moons: if has_bodies { r.path(bodies, "moons", "bodies")? } else { PathBuf::new() },
            vsop87: if has_bodies { r.path(bodies, "vsop87", "bodies")? } else { PathBuf::new() },
            earth_day: if has_bodies { r.path(earth, "day", "bodies.earth")? } else { PathBuf::new() },
            earth_night: if has_bodies { r.path(earth, "night", "bodies.earth")? } else { PathBuf::new() },
            earth_clouds: if has_bodies { r.path(earth, "clouds", "bodies.earth")? } else { PathBuf::new() },
            body_planets: {
                let mut list = Vec::new();
                for (name, given) in bodies.get("planetObjects").and_then(Value::as_object).into_iter().flatten() {
                    if name == "about" {
                        continue;
                    }
                    let style = given.get("style").and_then(Value::as_str).ok_or(format!("bodies.planetObjects.{name} has no \"style\""))?.to_owned();
                    let terrain = match given.get("terrain").and_then(Value::as_str) {
                        Some(spec) => {
                            let path = resolve(spec, &r.data, r.base);
                            let file = PathBuf::from(format!("{}.i16.gz", path.display()));
                            if !file.exists() {
                                return Err(format!("bodies.planetObjects.{name}.terrain: {} does not exist", file.display()));
                            }
                            Some(path)
                        }
                        None => None,
                    };
                    list.push((name.clone(), style, terrain));
                }
                list
            },
            planet_styles: {
                // The built-in styles, each changed by what the map says under its name.
                let mut styles = crate::ground::default_presets();
                for (name, given) in v.get("planetStyles").and_then(|s| s.get("styles")).and_then(Value::as_object).into_iter().flatten() {
                    let color = |key: &str, default: [f32; 3]| match given.get(key).and_then(Value::as_array) {
                        Some(c) if c.len() == 3 => [0, 1, 2].map(|k| c[k].as_f64().unwrap_or(default[k] as f64) as f32),
                        _ => default,
                    };
                    let base = styles.iter().find(|s| &s.name == name).cloned().unwrap_or_else(|| crate::ground::Preset { name: name.clone(), ..styles[0].clone() });
                    let style = crate::ground::Preset {
                        relief_km: number(given, "reliefKm", base.relief_km),
                        sea: number(given, "sea", base.sea),
                        roughness: number(given, "roughness", base.roughness),
                        mountains: number(given, "mountains", base.mountains),
                        ice_cap: number(given, "iceCap", base.ice_cap),
                        clouds: number(given, "clouds", base.clouds),
                        glow: number(given, "glow", base.glow),
                        low: color("low", base.low),
                        mid: color("mid", base.mid),
                        high: color("high", base.high),
                        sea_color: color("seaColor", base.sea_color),
                        air: color("air", base.air),
                        air_km: number(given, "airKm", base.air_km),
                        bands: number(given, "bands", base.bands),
                        contrast: number(given, "contrast", base.contrast),
                        detail_km: number(given, "detailKm", base.detail_km),
                        life: number(given, "life", base.life),
                        forest: color("forest", base.forest),
                        grass: color("grass", base.grass),
                        sand: color("sand", base.sand),
                        rock: color("rock", base.rock),
                        temp_equator_c: number(given, "tempEquatorC", base.temp_equator_c),
                        temp_pole_c: number(given, "tempPoleC", base.temp_pole_c),
                        temp_power: number(given, "tempPower", base.temp_power),
                        lapse_c_per_km: number(given, "lapseCPerKm", base.lapse_c_per_km),
                        rain: number(given, "rain", base.rain),
                        cloud_km: number(given, "cloudKm", base.cloud_km),
                        air_thickness: number(given, "airThickness", base.air_thickness),
                        name: name.clone(),
                    };
                    match styles.iter_mut().find(|s| &s.name == name) {
                        Some(slot) => *slot = style,
                        None => styles.push(style),
                    }
                }
                styles
            },
            rare: v
                .get("rare")
                .and_then(|r| r.get("classes"))
                .and_then(Value::as_array)
                .into_iter()
                .flatten()
                .filter_map(|c| Some((c.get("name")?.as_str()?.to_owned(), c.get("total").and_then(Value::as_f64), c.get("perStar").and_then(Value::as_f64))))
                .collect(),
            picture_catalogs: r.paths(v.get("pictures"), "pictures")?,
            ship: match v.get("ship").and_then(|s| s.get("folder")).and_then(Value::as_str) {
                Some(spec) => {
                    let folder = resolve(spec, &r.data, r.base);
                    if !folder.join("ship.json").exists() {
                        return Err(format!("ship.folder: {} has no ship.json", folder.display()));
                    }
                    Some(folder)
                }
                None => None,
            },
            picture_max_side: v.get("pictures").and_then(|p| p.get("maxSidePx")).and_then(Value::as_u64).map(|n| n as u32),
            data: r.data,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn resolves_paths() {
        let (data, base) = (Path::new("/d"), Path::new("/maps"));
        assert_eq!(resolve("$data/a/b.bin", data, base), Path::new("/d/a/b.bin"));
        assert_eq!(resolve("my/dust.raw.gz", data, base), Path::new("/maps/my/dust.raw.gz"));
        assert_eq!(resolve("/abs/x", data, base), Path::new("/abs/x"));
        // Nothing installed: the path is left as it is, to be reported as missing.
        assert_eq!(resolve("$data/tex/earth*", data, base), Path::new("/d/tex/earth*"));
    }

    #[test]
    fn reports_what_is_missing() {
        let e = Map::parse(r#"{"data": "/nowhere", "bodies": {"planets": "$data/p.json"}}"#, Path::new("."), None).err().unwrap();
        assert!(e.contains("bodies.earth"), "{e}");
        let e = Map::parse(r#"{"data": "/nowhere", "bodies": {"earth": {}, "planets": "$data/p.json"}}"#, Path::new("."), None).err().unwrap();
        assert!(e.contains("/nowhere/p.json does not exist"), "{e}");
        assert!(Map::parse("{", Path::new("."), None).is_err());
    }

    /// The galaxies made of three pictures are read only when the camera comes near enough for a
    /// cube with stars to reach them (none on loading the map), and let go when it is far; once read,
    /// their stars are made round the camera.
    #[test]
    fn xyz_galaxies_are_read_when_near() {
        let Some(home) = std::env::var_os("HOME") else { return };
        if !Path::new(&home).join(".gaiasky/data/default-data/planets.json").exists() {
            return;
        }
        let t = std::time::Instant::now();
        let map = Map::load(None, None).unwrap();
        eprintln!("the map in {:?}", t.elapsed());
        let model = map.star_field.clone().unwrap();
        let read = |m: &crate::procedural::Model| m.galaxies.iter().filter_map(|g| g.photo.as_ref()).filter(|l| l.get().is_some()).map(|l| l.name.clone()).collect::<Vec<_>>();
        assert!(model.galaxies.iter().filter(|g| g.photo.is_some()).count() >= 13);
        assert!(read(&model).is_empty(), "none read on loading: {:?}", read(&model));
        let stars = crate::procedural::star_reach();
        let wanted_at = |camera: DVec3| -> Vec<String> {
            model.galaxies.iter().filter_map(|g| g.photo.as_ref().map(|l| (g, l))).filter(|(g, l)| camera.distance(g.centre) < l.reach.0.hypot(l.reach.1) + stars).map(|(_, l)| l.name.clone()).collect()
        };
        // from the Sun
        let t = std::time::Instant::now();
        model.keep_near(DVec3::ZERO);
        eprintln!("cubes with stars reach {:.1} kpc; from the Sun read {:?} in {:?}", stars / 1000.0, read(&model), t.elapsed());
        assert_eq!(read(&model), wanted_at(DVec3::ZERO));
        // in Andromeda: read, and its stars there
        let andromeda = model.galaxies.iter().find(|g| g.name == "Andromeda").unwrap().clone();
        assert_eq!(andromeda.density(andromeda.centre), 0.0, "no stars before it is read");
        let at = andromeda.centre + andromeda.to_box.transpose() * DVec3::new(8000.0, 0.0, 0.0);
        model.keep_near(at);
        assert!(read(&model).contains(&"Andromeda Galaxy".to_owned()));
        assert!(!read(&model).contains(&"Large Magellanic Cloud".to_owned()), "the clouds let go: {:?}", read(&model));
        assert!(andromeda.density(at) > 0.0);
        let mut field = crate::procedural::Field::new(model.clone());
        while field.update(at, 100_000) {}
        let near = field.stars().filter(|(p, _)| p.distance(at) < 50.0).count();
        eprintln!("8 kpc out in Andromeda: {near} stars within 50 pc");
        assert!(near > 100);
        // and far again: let go
        model.keep_near(DVec3::ZERO);
        assert!(!read(&model).contains(&"Andromeda Galaxy".to_owned()));
    }

    /// The Milky Way from three pictures (starField.photo): as many stars round the Sun as
    /// measured, about as many in all as the formula has, and where the formula has them in the
    /// large (the disc thinning out from the centre, a bright middle, a thin plane).
    #[test]
    fn the_milky_way_from_three_pictures() {
        let Ok(map) = Map::load(None, None) else { return };
        let model = map.star_field.expect("a star field");
        let Some(home) = model.home.clone() else { panic!("starField.photo is not read") };
        let formula = crate::procedural::Model { home: None, galaxies: Vec::new(), fitted: Vec::new(), ..model.clone() };
        // round the Sun, in the plane: the measured density
        let mut near = 0.0f64;
        let mut n = 0.0f64;
        for i in -5..=5 {
            for k in -5..=5 {
                let p = DVec3::new(i as f64 * 80.0, 0.0, k as f64 * 80.0);
                if p.length() <= 400.0 {
                    near += home.density(p - model.centre);
                    n += 1.0;
                }
            }
        }
        let near: f64 = near / n;
        let systems = home.systems();
        println!("round the Sun {near:.3} systems/pc3 (measured {}), {systems:.3e} systems in all (the formula's ~2.7e11)", model.local_density);
        assert!((near / model.local_density - 1.0).abs() < 0.25, "{near} systems/pc3 round the Sun");
        assert!(systems > 1.0e11 && systems < 6.0e11, "{systems:e} systems");
        // in the large as the formula: along the line from the centre through the Sun, in the plane
        // (each a ring's mean, the pictures having arms where the formula has a smooth disc)
        let ring = |m: &crate::procedural::Model, r: f64, y: f64| -> f64 {
            (0..72).map(|a| {
                let t = a as f64 * std::f64::consts::TAU / 72.0;
                m.density(model.centre + DVec3::new(r * t.sin(), y, -r * t.cos()))
            }).sum::<f64>() / 72.0
        };
        for (r, y) in [(1000.0, 0.0), (4000.0, 0.0), (8000.0, 0.0), (12000.0, 0.0), (8000.0, 300.0), (8000.0, 1000.0)] {
            let (a, b) = (ring(&model, r, y), ring(&formula, r, y));
            println!("  {:5.0} pc out, {:4.0} pc up: pictures {a:.4}, formula {b:.4} systems/pc3 ({:.2})", r, y, a / b);
        }
        let (centre, sun) = (ring(&model, 1000.0, 0.0), ring(&model, 8000.0, 0.0));
        assert!(centre > 5.0 * sun, "the middle {centre} against the Sun's ring {sun}");
        // a thin plane, and the thick disc above it (measured, Juric et al. 2008: about 7% of the
        // plane's stars 1 kpc up at the Sun)
        let up = ring(&model, 8000.0, 1000.0) / sun;
        assert!(up > 0.02 && up < 0.2, "1 kpc up {up} of the plane");
    }

    /// With the Milky Way from its pictures the sky of the Earth is still the catalog's: no star
    /// made up round the Sun is brighter from there than the catalog's limit.
    #[test]
    fn the_milky_way_from_pictures_leaves_the_sky_of_the_earth() {
        let Ok(map) = Map::load(None, None) else { return };
        let model = map.star_field.expect("a star field");
        assert!(model.home.is_some());
        let mut field = crate::procedural::Field::new(model);
        while field.update(DVec3::ZERO, 10_000) {}
        let stars: Vec<_> = field.stars().collect();
        assert!(stars.len() > 1000, "{}", stars.len());
        let brightest = stars.iter().map(|(p, s)| s.abs_mag as f64 + 5.0 * (p.length() / 10.0).log10()).fold(f64::INFINITY, f64::min);
        println!("{} stars round the Sun, the brightest from the Earth {brightest:.2}", stars.len());
        assert!(brightest >= 9.5 - 1e-4, "a made-up star of magnitude {brightest} in the sky of the Earth");
    }

    #[test]
    fn the_built_in_map_reads() {
        let Some(home) = std::env::var_os("HOME") else { return };
        // The data is not part of the repository. Without it there is nothing to check.
        if !Path::new(&home).join(".gaiasky/data/default-data/planets.json").exists() {
            return;
        }
        let map = Map::load(None, None).unwrap();
        assert_eq!(map.star_catalogs.len(), 1);
        assert_eq!(map.volumes.len(), 6);
        assert!(map.volumes[0].light.is_some() && map.volumes[0].dust.is_none(), "the light of the galaxy has no dust");
        assert!(map.volumes[1].light.is_none() && map.volumes[1].dust_per_kpc > 0.0, "the dust of the galaxy has no light");
        assert_eq!(map.volumes[0].size_pc, 1000.0);
        assert!(map.volumes[2].light.is_none() && map.volumes[2].dust.is_some(), "the measured dust has no light");
        assert!(map.volumes[3].light.is_some() && map.volumes[3].fade_out_pc == Some((200.0, 1000.0)) && map.volumes[3].presence(DVec3::ZERO) == 1.0 && map.volumes[3].presence(DVec3::new(0.0, 6000.0, 0.0)) == 0.0);
        let exposure = map.exposure.as_ref().unwrap();
        assert!(exposure.at_pc == DVec3::ZERO && exposure.max > exposure.min);
        assert_eq!(exposure.balance(DVec3::ZERO), exposure.white_balance);
        assert_eq!(exposure.balance(DVec3::new(0.0, 0.0, 8000.0)), [1.0; 3]);
        assert_eq!(map.picture_catalogs.len(), 2);
        assert!(map.earth_day.is_dir());
        let field = map.star_field.unwrap();
        assert_eq!(field.centre, DVec3::new(0.0, 0.0, 8000.0));
        assert_eq!(field.arm_count, 2.0);
        // Andromeda: made of stars spread as the light of its volume, a trillion or so.
        let andromeda = &field.galaxies[0];
        assert_eq!(andromeda.name, "Andromeda");
        field.keep_near(andromeda.centre); // its pictures are read when near
        let middle = andromeda.density(andromeda.centre);
        assert!(middle > 1.0, "{middle} systems per cubic parsec in its middle");
        let at = |kpc: f64| andromeda.density(andromeda.centre + andromeda.to_box.transpose() * DVec3::new(kpc * 1000.0, 0.0, 0.0));
        // as the photograph has it (tools/photo_plane.py): the middle densest, the bright ring at
        // 10 kpc nearly as dense, thinning out beyond it
        assert!(at(10.0) > 0.05 && at(10.0) < middle && at(18.0) < at(10.0), "{} at 10 kpc, {} at 18, {middle} in the middle", at(10.0), at(18.0));
        assert!(at(30.0) < 1e-12, "nothing beyond its box");
        assert_eq!(field.density(andromeda.centre), middle, "the Milky Way adds nothing there");
        let mut stars = crate::procedural::Field::new(field.clone());
        let camera = andromeda.centre + andromeda.to_box.transpose() * DVec3::new(8000.0, 300.0, 0.0);
        while stars.update(camera, 10_000) {}
        let seen = stars.stars().filter(|(p, s)| s.abs_mag as f64 + 5.0 * (p.distance(camera) / 10.0).log10() < 8.8).count();
        assert!(seen > 1000, "{seen} of its stars to see from inside it");
        // of its trillion stars only those around the camera exist, in a few tens of megabytes
        let kept = stars.stars().count();
        eprintln!("inside Andromeda: {kept} stars kept ({:.0} MB), {seen} bright enough to see", kept as f64 * std::mem::size_of::<crate::procedural::Kept>() as f64 / 1e6);
        assert!(kept < crate::procedural::MAX_STARS, "{kept} stars kept");
    }

    /// Andromeda as the star formula spreads it, seen from the Earth, against its photograph:
    /// where its long axis points, how flat it looks, how long it is. Moments of the light: of the
    /// formula's glow points, and of the photograph's pixels (on the scale tools/photo_galaxy.py
    /// gave it: its disc reaches the grid's radius).
    #[test]
    fn andromeda_from_the_earth_looks_like_its_photograph() {
        let Ok(map) = Map::load(None, None) else { return };
        let photo_file = map.data.join("catalog-nbg/tex/extragal/M31.jpg");
        let Ok(photo) = image::open(&photo_file) else { return };
        let field = map.star_field.unwrap();
        let g = &field.galaxies[0];
        let to_eq = crate::render::galactic_to_equatorial();
        let centre = to_eq * g.centre;
        let sight = centre.normalize();
        let north = (DVec3::Y - sight * sight.y).normalize();
        let east = north.cross(sight);
        // (along east, along north, weight) -> centre, angle of the long axis from north through
        // east (degrees, 0..180), long and short spreads
        fn moments(points: &[(f64, f64, f64)]) -> (f64, f64, f64) {
            let w: f64 = points.iter().map(|p| p.2).sum();
            let (mx, my) = (points.iter().map(|p| p.0 * p.2).sum::<f64>() / w, points.iter().map(|p| p.1 * p.2).sum::<f64>() / w);
            let (mut xx, mut yy, mut xy) = (0.0, 0.0, 0.0);
            for &(x, y, k) in points {
                xx += (x - mx) * (x - mx) * k;
                yy += (y - my) * (y - my) * k;
                xy += (x - mx) * (y - my) * k;
            }
            let (xx, yy, xy) = (xx / w, yy / w, xy / w);
            let half = ((xx - yy) * (xx - yy) / 4.0 + xy * xy).sqrt();
            let (long, short) = ((xx + yy) / 2.0 + half, (xx + yy) / 2.0 - half);
            let angle = (0.5 * (2.0 * xy).atan2(yy - xx)).to_degrees().rem_euclid(180.0);
            (angle, long.sqrt(), short.sqrt())
        }
        // the formula: its glow points, seen from the Earth (the Sun, at this distance)
        let (places, _, _) = g.glow(300_000, 7);
        let formula: Vec<(f64, f64, f64)> = places
            .iter()
            .map(|q| {
                let d = (to_eq * (g.centre + g.to_box.transpose() * *q)).normalize();
                (d.dot(east), d.dot(north), 1.0)
            })
            .collect();
        // the photograph: north up, east to the left, light in linear units
        let photo = photo.to_rgb32f();
        let (w, h) = photo.dimensions();
        let lum: Vec<f64> = photo.pixels().map(|p| p.0.iter().map(|&c| (c as f64).powf(2.2)).sum::<f64>() / 3.0).collect();
        let top = lum.iter().cloned().fold(0.0, f64::max);
        // its middle and its scale, as tools/photo_plane.py found them (in the map's metadata)
        let meta = std::fs::read_to_string(Path::new(env!("CARGO_MANIFEST_DIR")).join("../map/grids/andromeda_plane_metadata.txt")).unwrap();
        let numbers = |key: &str| -> Vec<f64> { meta.lines().find_map(|l| l.strip_prefix(key)).unwrap().split_whitespace().filter_map(|t| t.parse().ok()).collect() };
        let (cx, cy) = (numbers("Photograph middle px:")[0], numbers("Photograph middle px:")[1]);
        let rad_per_px = numbers("Photograph kpc per px:")[0] * 1000.0 / centre.length();
        let pixels: Vec<(f64, f64, f64)> = (0..lum.len())
            .filter(|&i| lum[i] > 0.02 * top)
            .map(|i| (-((i as u32 % w) as f64 - cx), -((i as u32 / w) as f64 - cy), lum[i]))
            .collect();
        let (pa, pl, ps) = moments(&pixels);
        let (fa, fl, fs) = moments(&formula);
        eprintln!("photograph: middle ({cx}, {cy})");
        eprintln!("photograph {w}x{h}: long axis {pa:.1} deg, {:.3} by {:.3} deg; formula: {fa:.1} deg, {:.3} by {:.3} deg", (pl * rad_per_px).to_degrees(), (ps * rad_per_px).to_degrees(), fl.to_degrees(), fs.to_degrees());
        // pixel by pixel: the formula's light on the photograph's pixels (in blocks of 4), against
        // the photograph's, both as the square root of the light (as the eye sees it, roughly)
        let block = 4u32;
        let (bw, bh) = ((w / block) as usize, (h / block) as usize);
        let similarity = |galaxy: &crate::procedural::GridGalaxy| {
            let (places, _, _) = galaxy.glow(3_000_000, 11);
            let mut img = vec![0.0f64; bw * bh];
            for q in &places {
                let d = (to_eq * (galaxy.centre + galaxy.to_box.transpose() * *q)).normalize();
                let (col, row) = (cx - d.dot(east) / rad_per_px, cy - d.dot(north) / rad_per_px);
                if col >= 0.0 && row >= 0.0 && (col as u32) < w && (row as u32) < h {
                    img[(row as usize / block as usize) * bw + col as usize / block as usize] += 1.0;
                }
            }
            let mut pic = vec![0.0f64; bw * bh];
            for (i, &l) in lum.iter().enumerate() {
                let (x, y) = (i % w as usize, i / w as usize);
                if x / (block as usize) < bw && y / (block as usize) < bh {
                    pic[(y / block as usize) * bw + x / block as usize] += l;
                }
            }
            if let Ok(dir) = std::env::var("DUMP") {
                let save = |v: &[f64], name: &str| {
                    let m = v.iter().cloned().fold(0.0, f64::max).max(1e-12);
                    let bytes: Vec<u8> = v.iter().map(|x| ((x / m).sqrt() * 255.0) as u8).collect();
                    image::GrayImage::from_raw(bw as u32, bh as u32, bytes).unwrap().save(format!("{dir}/{name}_{}.png", galaxy.plane.n)).unwrap();
                };
                save(&img, "formula");
                save(&pic, "photo");
            }
            let (a, b): (Vec<f64>, Vec<f64>) = img.iter().zip(&pic).filter(|(i, p)| **i > 0.0 || **p > 0.02 * top * 16.0).map(|(i, p)| (i.sqrt(), p.sqrt())).unzip();
            let (ma, mb) = (a.iter().sum::<f64>() / a.len() as f64, b.iter().sum::<f64>() / b.len() as f64);
            let cov: f64 = a.iter().zip(&b).map(|(x, y)| (x - ma) * (y - mb)).sum();
            let va: f64 = a.iter().map(|x| (x - ma).powi(2)).sum();
            let vb: f64 = b.iter().map(|y| (y - mb).powi(2)).sum();
            cov / (va * vb).sqrt()
        };
        let turn = (pa - fa).abs().min(180.0 - (pa - fa).abs());
        if g.photo.is_some() {
            // spread by three pictures (photo_volume.rs) with their light put back: its stars are
            // far more to the middle than a photograph's squeezed light; that it is the
            // photograph seen from the Earth, light as drawn, is photo_volume.rs's test. Here: where
            // its long axis points and how flat it looks, as the photograph.
            eprintln!("spread by three pictures: like the photograph pixel by pixel {:.3} (its light put back)", similarity(g));
            assert!(turn < 3.0, "the long axis is turned by {turn:.1} degrees");
            assert!(((fs / fl) / (ps / pl) - 1.0).abs() < 0.15, "as flat: {:.3} against {:.3}", fs / fl, ps / pl);
            return;
        }
        let plain = crate::procedural::GridGalaxy { plane: std::sync::Arc::new(Default::default()), ..g.clone() };
        let (with_map, without) = (similarity(g), similarity(&plain));
        eprintln!("pixel by pixel like the photograph: {with_map:.3} with the face-on map, {without:.3} with the grid alone");
        assert!(with_map > 0.98 && with_map > without, "{with_map:.3} against {without:.3}");
        assert!(turn < 3.0, "the long axis is turned by {turn:.1} degrees");
        assert!(((fs / fl) / (ps / pl) - 1.0).abs() < 0.15, "as flat: {:.3} against {:.3}", fs / fl, ps / pl);
        assert!((fl / (pl * rad_per_px) - 1.0).abs() < 0.15, "as long: {:.4} against {:.4} rad", fl, pl * rad_per_px);
        let _ = h;
    }

    /// Seen from 200 000 light years away on the Earth's side, Andromeda fills a good part of the
    /// sky: its glow points spread over tens of degrees.
    #[test]
    fn andromeda_is_large_from_near_it() {
        let Ok(map) = Map::load(None, None) else { return };
        let field = map.star_field.unwrap();
        let g = &field.galaxies[0];
        let to_eq = crate::render::galactic_to_equatorial();
        let centre = to_eq * g.centre;
        let eye = centre - centre.normalize() * 200_000.0 / 3.2616;
        let (places, _, _) = g.glow(100_000, 3);
        let mut off: Vec<f64> = places.iter().map(|q| (to_eq * (g.centre + g.to_box.transpose() * *q) - eye).angle_between(centre - eye).to_degrees()).collect();
        off.sort_by(f64::total_cmp);
        let reach = off[(off.len() as f64 * 0.98) as usize];
        eprintln!("from 200 000 ly: 98% of its light within {reach:.1} degrees of its middle");
        assert!((10.0..30.0).contains(&reach), "{reach}");
    }

}
