//! Places by address, as a postal address goes from the province to the street: a galaxy, a star
//! of it, and a planet of that star. "Andromeda Galaxy/1/2" is the second planet of the first star
//! of Andromeda; "Milky Way/Sirius" the star Sirius; "Milky Way/Sun/3" the Earth.
//!
//! The stars of a galaxy are numbered from 1 to the number of star systems it has: the model's
//! own count (the Milky Way's from its formula in procedural.rs, about 2.7e11; the others' from
//! map/galaxies/fitted.json). Star 1 is the galaxy's central star, a hot young star of the cluster
//! round its black hole, a quarter of a parsec from the middle (as lib.rs's CORE_STARS has for
//! Andromeda). Every other number is a star drawn from the galaxy's shape, where its light is
//! (fitted.rs `glow`; the Milky Way's formula, `milky_way_place`), what it is like from the share
//! of stars of each brightness (procedural.rs), its companions (systems.rs) and its planets
//! (planets.rs). All of it follows from the galaxy and the number alone: nothing is stored, and
//! the same address is the same star for everybody.
//!
//! The made-up stars of other galaxies take the brightness shares of the stars near the Sun (the
//! only ones measured star by star); their planets follow what is known of planets in the Milky
//! Way. In the Milky Way a star may also be named: the Sun's planets are the real ones (lib.rs,
//! the bodies), a named star of the catalog gets made-up planets from its name.

use glam::{DMat3, DQuat, DVec3};

use crate::fitted::ParamGalaxy;
use crate::procedural::{self, Kind, Model, Rng, Star, StarId};
use crate::{planets, stars, systems};

pub const MILKY_WAY: &str = "Milky Way";
pub const AU_PER_PC: f64 = 206_264.806;
const LY_PER_PC: f64 = 3.261_563_777;

/// Other names of the galaxies, Korean first (the label shown), for the address bar.
const NAMES: &[(&str, &[&str])] = &[
    (MILKY_WAY, &["우리 은하", "은하수", "Milky Way Galaxy", "Galaxy"]),
    ("Andromeda Galaxy", &["안드로메다 은하", "안드로메다", "Andromeda", "M31", "Messier 31", "NGC 224"]),
    ("Large Magellanic Cloud", &["대마젤란 은하", "대마젤란운"]),
    ("Small Magellanic Cloud", &["소마젤란 은하", "소마젤란운"]),
    ("Triangulum spiral", &["삼각형자리 은하", "삼각형자리", "Triangulum Galaxy"]),
    ("Surfboard Galaxy", &["서핑보드 은하"]),
    ("Vacuum Cleaner Galaxy", &["진공청소기 은하"]),
    ("Whirlpool Galaxy", &["소용돌이 은하"]),
    ("M58", &["M58 은하"]),
    ("M61", &["M61 은하"]),
    ("Bode's Galaxy", &["보데 은하"]),
    ("M82", &["시가 은하"]),
    ("NGC 2403", &["NGC 2403 은하"]),
    ("M65", &["M65 은하"]),
];

/// Where the central star is from the galaxy's middle, in its own axes (pc), and what it is: a
/// young A star, 2.4 times the Sun across, 9,800 K (lib.rs CORE_STARS, Andromeda's P3 cluster).
pub const CORE_OFFSET_PC: [f64; 3] = [0.25, 0.08, -0.15];
const CORE_ABS_MAG: f64 = 0.8;
/// A planet "near the middle" of a galaxy is of a star this near it, pc.
pub const CORE_REACH_PC: f64 = 500.0;

/// From galactic axes (those of procedural.rs and fitted.rs: the Sun at the origin, +Z toward the
/// galactic centre, +Y to the north galactic pole) to equatorial (lib.rs: +Z to the equinox, +X
/// to RA 90°, +Y north).
pub fn galactic_to_equatorial() -> DMat3 {
    let (r, q, p) = (32.93192f64, 27.12825f64, 192.85948f64);
    DMat3::from_rotation_y((90.0 + p).to_radians()) * DMat3::from_rotation_z((90.0 - q).to_radians()) * DMat3::from_rotation_y((-r).to_radians())
}

/// Names compared without case or spaces.
pub fn key(s: &str) -> String {
    s.chars().filter(|c| !c.is_whitespace() && *c != '\'').flat_map(char::to_lowercase).collect()
}

enum Shape {
    MilkyWay,
    Fitted(usize),
    /// A galaxy made from a grid of light (Andromeda), and what its stars are drawn by, made the
    /// first time one is asked for (it goes through the whole grid).
    Grid(usize, std::sync::OnceLock<(Vec<f64>, f64)>),
}

/// A galaxy one can go to by address.
pub struct GalaxyEntry {
    /// The name in addresses.
    pub name: String,
    /// The Korean name shown, and every name it is found by.
    pub label: String,
    pub names: Vec<String>,
    /// Its middle (galactic, parsecs from the Sun) and the rows that turn galactic axes into its
    /// own (x along the disc, y its pole).
    pub centre_pc: DVec3,
    pub to_box: DMat3,
    /// How many star systems it has: the numbers of its stars go up to this.
    pub systems: u64,
    shape: Shape,
}

impl GalaxyEntry {
    pub fn light_years(&self) -> f64 {
        self.centre_pc.length() * LY_PER_PC
    }
}

/// The galaxies of the model, by address.
pub struct Atlas {
    pub model: Model,
    pub galaxies: Vec<GalaxyEntry>,
}

/// A star of an address: by its number, or (in the Milky Way) by its name.
#[derive(Clone, Debug, PartialEq)]
pub enum StarRef {
    Number(u64),
    Named(String),
}

#[derive(Clone, Debug, PartialEq)]
pub struct Address {
    pub galaxy: usize,
    pub star: StarRef,
    pub planet: Option<u32>,
}

/// A star found by its address, with what goes round it.
pub struct MadeStar {
    /// Equatorial, parsecs from the Sun (lib.rs's axes).
    pub position_pc: DVec3,
    pub star: Star,
    pub about: stars::Description,
    pub system: systems::System,
    pub planets: Vec<planets::Planet>,
    /// The plane its planets go round in (their orbits are made in the system's own).
    pub plane: DQuat,
}

impl MadeStar {
    /// Where planet `k` (from 1, the innermost first) is on `days` (since J2000): equatorial
    /// parsecs from the Sun.
    pub fn planet_pc(&self, k: u32, days: f64) -> Option<DVec3> {
        let i = (k as usize).checked_sub(1)?;
        if i >= self.planets.len() {
            return None;
        }
        let places = planets::places_au(&self.system, &self.planets, days / 365.25);
        Some(self.position_pc + self.plane * places[i].0 / AU_PER_PC)
    }
}

impl Atlas {
    /// The Milky Way of `model` and its fitted galaxies.
    pub fn new(model: Model) -> Atlas {
        let mut galaxies = Vec::new();
        // its name in addresses (the one of NAMES it is known by), its label and all its names
        let named = |name: &str, own: &[String]| {
            let mut names: Vec<String> = own.to_vec();
            let (mut canonical, mut label) = (name.to_string(), name.to_string());
            if let Some((n, more)) = NAMES.iter().find(|(n, more)| key(n) == key(name) || more.iter().any(|m| key(m) == key(name))) {
                canonical = n.to_string();
                label = more[0].to_string();
                names.extend(more.iter().map(|s| s.to_string()));
            }
            for n in [name.to_string(), canonical.clone()] {
                if !names.iter().any(|m| *m == n) {
                    names.insert(0, n);
                }
            }
            (canonical, label, names)
        };
        let (name, label, names) = named(MILKY_WAY, &[]);
        galaxies.push(GalaxyEntry { name, label, names, centre_pc: model.centre, to_box: DMat3::IDENTITY, systems: milky_way_systems(&model) as u64, shape: Shape::MilkyWay });
        for (i, g) in model.galaxies.iter().enumerate() {
            if g.light_only {
                continue;
            }
            let (name, label, names) = named(&g.name, &[]);
            galaxies.push(GalaxyEntry { name, label, names, centre_pc: g.centre, to_box: g.to_box, systems: g.systems.max(1.0) as u64, shape: Shape::Grid(i, std::sync::OnceLock::new()) });
        }
        for (i, g) in model.fitted.iter().enumerate() {
            let (name, label, names) = named(&g.name, &g.names);
            galaxies.push(GalaxyEntry { name, label, names, centre_pc: g.centre, to_box: g.to_box, systems: g.systems.max(1.0) as u64, shape: Shape::Fitted(i) });
        }
        Atlas { model, galaxies }
    }

    /// The galaxy of a name (any of its names, whole).
    pub fn galaxy(&self, name: &str) -> Option<usize> {
        let k = key(name);
        self.galaxies.iter().position(|g| g.names.iter().any(|n| key(n) == k))
    }

    /// The galaxies whose names hold what is typed, those that begin with it first, then the
    /// nearest; at most `limit`.
    pub fn suggest(&self, typed: &str, limit: usize) -> Vec<usize> {
        let k = key(typed);
        let mut hits: Vec<(bool, f64, usize)> = self
            .galaxies
            .iter()
            .enumerate()
            .filter(|(_, g)| k.is_empty() || g.names.iter().any(|n| key(n).contains(&k)))
            .map(|(i, g)| (!g.names.iter().any(|n| key(n).starts_with(&k)), g.centre_pc.length(), i))
            .collect();
        hits.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.total_cmp(&b.1)));
        hits.into_iter().take(limit).map(|h| h.2).collect()
    }

    /// Reads an address: "galaxy/star/planet", the star and the planet optional. A star is a
    /// number from 1 to the galaxy's count (commas and spaces in it are fine), or in the Milky Way
    /// a name (checked by lib.rs, which has the catalog).
    pub fn parse(&self, address: &str) -> Result<Address, String> {
        let parts: Vec<&str> = address.split('/').map(str::trim).collect();
        if parts.len() < 2 || parts.len() > 3 {
            return Err("an address is galaxy/star or galaxy/star/planet".into());
        }
        let galaxy = self.galaxy(parts[0]).ok_or(format!("no galaxy called {}", parts[0]))?;
        let g = &self.galaxies[galaxy];
        let digits: String = parts[1].chars().filter(|c| !matches!(c, ',' | '_' | ' ')).collect();
        let star = if !digits.is_empty() && digits.chars().all(|c| c.is_ascii_digit()) {
            let n: u64 = digits.parse().map_err(|_| format!("star {} is past the last of {} ({})", parts[1], g.name, g.systems))?;
            if n == 0 || n > g.systems {
                return Err(format!("{} has stars 1 to {}", g.name, g.systems));
            }
            StarRef::Number(n)
        } else if galaxy == 0 && !parts[1].is_empty() {
            StarRef::Named(parts[1].to_string())
        } else {
            return Err(format!("a star of {} is a number, 1 to {}", g.name, g.systems));
        };
        let planet = match parts.get(2) {
            None => None,
            Some(p) => match p.parse::<u32>() {
                Ok(k) if k >= 1 => Some(k),
                _ => return Err(format!("a planet is a number from 1, not {p}")),
            },
        };
        Ok(Address { galaxy, star, planet })
    }

    /// The address written the one way, for the server and the route.
    pub fn canonical(&self, a: &Address) -> String {
        let star = match &a.star {
            StarRef::Number(n) => n.to_string(),
            StarRef::Named(s) => s.clone(),
        };
        match a.planet {
            Some(k) => format!("{}/{star}/{k}", self.galaxies[a.galaxy].name),
            None => format!("{}/{star}", self.galaxies[a.galaxy].name),
        }
    }

    /// A star by its number in a galaxy.
    pub fn star(&self, galaxy: usize, n: u64) -> MadeStar {
        let g = &self.galaxies[galaxy];
        let seed = seed_of(&g.name, n);
        let mut rng = Rng::seeded(seed);
        let to_eq = galactic_to_equatorial();
        let (place, abs_mag, kind) = if n == 1 {
            (g.centre_pc + g.to_box.transpose() * DVec3::from_array(CORE_OFFSET_PC), CORE_ABS_MAG, Kind::MainSequence)
        } else {
            let place = match &g.shape {
                Shape::MilkyWay => milky_way_place(&self.model, &mut rng),
                Shape::Grid(i, table) => {
                    let grid = &self.model.galaxies[*i];
                    let (own, _) = grid.place(table.get_or_init(|| grid.cumulative()), &mut Rng::seeded(seed ^ 0x6C6F_7767));
                    grid.centre + grid.to_box.transpose() * own
                }
                Shape::Fitted(i) => {
                    let f: &ParamGalaxy = &self.model.fitted[*i];
                    let own = f.glow(1, seed).0.first().copied().unwrap_or(DVec3::ZERO);
                    f.centre + f.to_box.transpose() * own
                }
            };
            let (m, kind) = self.model.star_of_share(rng.next(), rng.next(), rng.next());
            (place, m, kind)
        };
        let star = Star { position: place, abs_mag: abs_mag as f32, color: procedural::color(abs_mag, kind), kind, id: StarId { level: 255, cell: [galaxy as i64, (n >> 32) as i64, (n & 0xFFFF_FFFF) as i64], index: 0 } };
        self.made(star, to_eq * place, seed, n == 1)
    }

    /// A star of the Milky Way's catalog, by its name: its place (equatorial pc) and its brightness
    /// are the catalog's; its planets are made up from its name.
    pub fn named_star(&self, name: &str, position_pc: DVec3, abs_mag: f64) -> MadeStar {
        let seed = seed_of(&key(name), 0);
        let kind = if abs_mag < -1.5 { Kind::Giant } else { Kind::MainSequence };
        let galactic = galactic_to_equatorial().transpose() * position_pc;
        let star = Star { position: galactic, abs_mag: abs_mag as f32, color: procedural::color(abs_mag, kind), kind, id: StarId { level: 254, cell: [0; 3], index: 0 } };
        // single: a real star's companions are stars of the catalog of their own
        self.made(star, position_pc, seed, true)
    }

    fn made(&self, star: Star, position_pc: DVec3, seed: u64, single: bool) -> MadeStar {
        let about = stars::describe(&star, seed);
        // alone: the central star (in the crowd round a black hole wide companions do not last),
        // and a named star (its real companions are named stars of their own)
        const ALONE: [systems::Multiplicity; 1] = [systems::Multiplicity { min_mass: 0.0, double: 0.0, triple: 0.0, more: 0.0 }];
        let table: &[systems::Multiplicity] = if single { &ALONE } else { &self.model.companions };
        let system = systems::system_of(&star, seed, table);
        let planets = planets::planets_of(&system, seed, &self.model.planets);
        let mut rng = Rng::seeded(seed ^ 0x504c_414e_45);
        let axis = DVec3::new(rng.next() - 0.5, rng.next() - 0.5, rng.next() - 0.5).try_normalize().unwrap_or(DVec3::Y);
        let plane = DQuat::from_axis_angle(axis, rng.next() * std::f64::consts::TAU);
        MadeStar { position_pc, star, about, system, planets, plane }
    }

    /// Of a star's planets, the one to go to (from 1): one in the zone where water stays liquid if
    /// it has one, else the largest; None if it has none.
    pub fn best_planet(star: &MadeStar) -> Option<u32> {
        let p = &star.planets;
        let temperate = p.iter().position(|q| q.kind == planets::Kind::Temperate);
        let largest = (0..p.len()).max_by(|&a, &b| p[a].radius.total_cmp(&p[b].radius));
        temperate.or(largest).map(|i| i as u32 + 1)
    }

    /// A planet near the middle of a galaxy, to go to: of its central star if it has planets, else
    /// of the first star by number within CORE_REACH_PC of the middle that has. (star, planet)
    pub fn core_planet(&self, galaxy: usize) -> Option<(u64, u32)> {
        let g = &self.galaxies[galaxy];
        let to_gal = galactic_to_equatorial().transpose();
        let last = g.systems.min(200_000);
        (1..=last).find_map(|n| {
            let s = self.star(galaxy, n);
            if n > 1 && (to_gal * s.position_pc - g.centre_pc).length() > CORE_REACH_PC {
                return None;
            }
            Some((n, Self::best_planet(&s)?))
        })
    }

    /// A star of the Milky Way on its far side: beyond the centre from the Sun, in the disc, the
    /// first by number. For the tour (hud.gd).
    pub fn far_side(&self) -> u64 {
        let c = self.model.centre;
        (2..100_000)
            .find(|&n| {
                let p = galactic_to_equatorial().transpose() * self.star(0, n).position_pc;
                let q = p - c;
                // past the centre, along the line from the Sun, 6 to 10 kpc beyond it, in the plane
                q.dot(c.normalize()) > 6_000.0 && q.dot(c.normalize()) < 10_000.0 && q.y.abs() < 300.0 && (q - c.normalize() * q.dot(c.normalize())).length() < 3_000.0
            })
            .unwrap_or(2)
    }
}

/// How many star systems the Milky Way's formula has (procedural.rs Model::density), added up as
/// formulas: the disc 4π ρ☉ e^(R☉/h) h² z (as its test says) and the bulge π^1.5 R³ / √flattening.
pub fn milky_way_systems(m: &Model) -> f64 {
    disc_systems(m) + bulge_systems(m)
}

fn disc_systems(m: &Model) -> f64 {
    std::f64::consts::TAU * m.local_density * 2.0 * m.disc_height * (m.centre.length() / m.disc_length).exp() * m.disc_length.powi(2)
}

fn bulge_systems(m: &Model) -> f64 {
    m.bulge_density * std::f64::consts::PI.powf(1.5) * m.bulge_radius.powi(3) / m.bulge_flattening.sqrt()
}

/// A place in the Milky Way drawn as its formula has its stars (galactic pc): the bulge or the
/// disc as their shares; in the disc the radius from r e^(-r/h), the arms and the edge kept as the
/// density says, the height from e^(-|y|/z).
pub fn milky_way_place(m: &Model, rng: &mut Rng) -> DVec3 {
    let bulge_share = bulge_systems(m) / milky_way_systems(m);
    let normal = |rng: &mut Rng| (-2.0 * rng.next().max(1e-12).ln()).sqrt() * (std::f64::consts::TAU * rng.next()).cos();
    for _ in 0..1000 {
        if rng.next() < bulge_share {
            // e^(-(x² + z² + f y²) / R²): normal with σ = R/√2 across, R/√(2f) up
            let s = m.bulge_radius / 2f64.sqrt();
            return m.centre + DVec3::new(normal(rng) * s, normal(rng) * s / m.bulge_flattening.sqrt(), normal(rng) * s);
        }
        let r = -m.disc_length * (rng.next().max(1e-12) * rng.next().max(1e-12)).ln();
        let [inner, outer] = m.disc_edge;
        if rng.next() > 1.0 - ((r - inner) / (outer - inner).max(1.0)).clamp(0.0, 1.0) {
            continue;
        }
        let angle = rng.next() * std::f64::consts::TAU;
        if r > m.arm_start && m.arm_strength != 0.0 {
            let arms = 1.0 + m.arm_strength * (m.arm_count * (angle - (r / m.arm_start).ln() / m.arm_pitch)).cos();
            if rng.next() * (1.0 + m.arm_strength) > arms {
                continue;
            }
        }
        let y = -m.disc_height * rng.next().max(1e-12).ln() * if rng.next() < 0.5 { -1.0 } else { 1.0 };
        // as in Model::density: the angle is atan2(x, z)
        return m.centre + DVec3::new(r * angle.sin(), y, r * angle.cos());
    }
    m.centre
}

/// The number of everything made up about a star of a galaxy.
fn seed_of(galaxy: &str, n: u64) -> u64 {
    let mut h = 0xCBF2_9CE4_8422_2325u64;
    for b in key(galaxy).bytes() {
        h = (h ^ b as u64).wrapping_mul(0x0100_0000_01B3);
    }
    let mut x = h ^ n.wrapping_mul(0x9E37_79B9_7F4A_7C15);
    x = (x ^ (x >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    x = (x ^ (x >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    x ^ (x >> 31)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn atlas() -> Atlas {
        Atlas::new(Model::default())
    }

    #[test]
    fn the_milky_way_counts_as_its_model() {
        let m = Model::default();
        let n = milky_way_systems(&m);
        // procedural.rs: "270 billion systems" (its sum over the galaxy, 3% either way); the sum
        // as formulas leaves out the disc's edge and arms, a little more
        assert!((n / 2.7e11 - 1.0).abs() < 0.12, "{n:e}");
        assert_eq!(atlas().galaxies[0].systems, n as u64);
    }

    #[test]
    fn the_galactic_centre_is_where_the_sky_has_it() {
        // Sgr A*: RA 266.405°, Dec -28.936°
        let (ra, dec) = (266.405f64.to_radians(), (-28.936f64).to_radians());
        let seen = DVec3::new(dec.cos() * ra.sin(), dec.sin(), dec.cos() * ra.cos());
        let got = (galactic_to_equatorial() * Model::default().centre).normalize();
        assert!(got.angle_between(seen).to_degrees() < 0.5, "{got} against {seen}");
    }

    #[test]
    fn reads_and_writes_addresses() {
        let a = atlas();
        let x = a.parse("우리 은하 / 1,000 / 2").unwrap();
        assert_eq!(x, Address { galaxy: 0, star: StarRef::Number(1000), planet: Some(2) });
        assert_eq!(a.canonical(&x), "Milky Way/1000/2");
        assert_eq!(a.parse("milkyway/Sirius").unwrap().star, StarRef::Named("Sirius".into()));
        assert!(a.parse("Milky Way/0").is_err(), "stars from 1");
        assert!(a.parse(&format!("Milky Way/{}", a.galaxies[0].systems + 1)).is_err(), "no more than it has");
        assert!(a.parse("Milky Way/99999999999999999999999").is_err());
        assert!(a.parse("Milky Way/5/0").is_err(), "planets from 1");
        assert!(a.parse("Nowhere/5").is_err());
        assert!(a.parse("Milky Way").is_err(), "a galaxy alone is not a place");
    }

    #[test]
    fn the_same_number_is_the_same_star() {
        let a = atlas();
        for n in [2u64, 77, 123_456_789, a.galaxies[0].systems] {
            let (s, t) = (a.star(0, n), a.star(0, n));
            assert_eq!(s.position_pc, t.position_pc);
            assert_eq!(s.about, t.about);
            assert_eq!(s.planets.len(), t.planets.len());
            assert_eq!(s.planet_pc(1, 9000.0), t.planet_pc(1, 9000.0));
        }
        assert_ne!(a.star(0, 2).position_pc, a.star(0, 3).position_pc);
    }

    #[test]
    fn stars_are_where_the_milky_way_has_them() {
        let a = atlas();
        let m = &a.model;
        let to_gal = galactic_to_equatorial().transpose();
        let (mut bulge, mut thin, mut inside) = (0, 0, 0);
        let n = 4000;
        for i in 2..n + 2 {
            let q = to_gal * a.star(0, i).position_pc - m.centre;
            let r = (q.x * q.x + q.z * q.z).sqrt();
            bulge += (q.length() < 2.0 * m.bulge_radius) as usize;
            thin += (q.y.abs() < 3.0 * m.disc_height) as usize;
            inside += (r < m.disc_edge[1]) as usize;
        }
        assert_eq!(inside, n as usize, "none past the disc's edge");
        assert!(thin as f64 > 0.9 * n as f64, "most in the thin disc: {thin}");
        let share = bulge_systems(m) / milky_way_systems(m);
        assert!(bulge as f64 / n as f64 > 0.5 * share, "the bulge has its share: {bulge} of {n}, {share:.3}");
        // most stars are faint red dwarfs, as near the Sun
        let dwarfs = (2..2002u64).filter(|&i| a.star(0, i).about.class.starts_with('M')).count();
        assert!((1100..1700).contains(&dwarfs), "{dwarfs} M dwarfs of 2000");
    }

    #[test]
    fn the_central_star_is_in_the_middle() {
        let a = atlas();
        let s = a.star(0, 1);
        let centre = galactic_to_equatorial() * a.model.centre;
        let off = (s.position_pc - centre).length();
        assert!((off - DVec3::from_array(CORE_OFFSET_PC).length()).abs() < 1e-9, "{off} pc from the middle");
        assert!(s.about.class.starts_with('A') || s.about.class.starts_with('B'), "a hot young star: {}", s.about.class);
        assert_eq!(s.system.stars(), 1);
    }

    #[test]
    fn planets_go_round_their_star() {
        let a = atlas();
        let with = (2..500u64).map(|n| a.star(0, n)).find(|s| s.planets.len() >= 3).expect("a star with planets");
        let p1 = with.planet_pc(1, 0.0).unwrap();
        let p3 = with.planet_pc(3, 0.0).unwrap();
        assert!((p1 - with.position_pc).length() < (p3 - with.position_pc).length() * 1.5, "the first is the innermost, about");
        // it moves along its orbit, and stays as far
        let later = with.planet_pc(1, 30.0).unwrap();
        assert_ne!(p1, later);
        let a1 = with.planets[0].orbit.semi_major_au;
        let d = (p1 - with.position_pc).length() * AU_PER_PC;
        assert!(d < a1 * 2.5 + 50.0, "{d} AU from its star, its orbit {a1} AU");
        assert!(with.planet_pc(0, 0.0).is_none() && with.planet_pc(with.planets.len() as u32 + 1, 0.0).is_none());
        // how many stars have planets, about as planets.rs means
        let some = (2..1002u64).filter(|&n| !a.star(0, n).planets.is_empty()).count();
        assert!(some > 500, "{some} of 1000 stars with planets");
    }

    #[test]
    fn suggests_galaxies_as_one_types() {
        let mut model = Model::default();
        model.fitted = Vec::new();
        let a = Atlas::new(model);
        assert_eq!(a.suggest("우리", 5), vec![0]);
        assert_eq!(a.suggest("milky", 5), vec![0]);
        assert_eq!(a.suggest("", 5), vec![0]);
        assert!(a.suggest("zzz", 5).is_empty());
        assert_eq!(a.galaxies[0].label, "우리 은하");
    }

    #[test]
    fn a_planet_near_the_middle() {
        let a = atlas();
        let (n, k) = a.core_planet(0).expect("a planet near the Milky Way's middle");
        let s = a.star(0, n);
        let to_gal = galactic_to_equatorial().transpose();
        assert!((to_gal * s.position_pc - a.model.centre).length() <= CORE_REACH_PC, "star {n}");
        assert!(k >= 1 && k as usize <= s.planets.len());
        assert_eq!(a.core_planet(0), Some((n, k)), "the same every time");
        // the planet chosen: a temperate one if there is one, else none larger
        let p = &s.planets;
        let chosen = &p[k as usize - 1];
        assert!(chosen.kind == planets::Kind::Temperate || p.iter().all(|q| q.kind != planets::Kind::Temperate && q.radius <= chosen.radius));
    }

    #[test]
    fn the_far_side_is_beyond_the_centre() {
        let a = atlas();
        let n = a.far_side();
        let p = a.star(0, n).position_pc;
        let d = p.length();
        let centre = (galactic_to_equatorial() * a.model.centre).length();
        assert!(d > centre + 5_000.0, "star {n} at {d:.0} pc, the centre at {centre:.0} pc");
    }
}
