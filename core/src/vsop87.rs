//! The VSOP87A planetary theory: heliocentric rectangular positions of the planets.
//!
//! This mirrors `gaiasky.util.coord.vsop87.VSOP87Binary`. The coefficients come from the
//! `vsop87a.bin` file of the Gaia Sky data, a flat array of little-endian doubles.

use std::path::Path;

use glam::DVec3;

use crate::astro::AU_KM;

/// For each of X, Y and Z, six (offset, count) pairs, one per power of time. Offsets count terms
/// of three doubles.
type Offsets = [usize; 36];

const MERCURY: Offsets = [18367, 1449, 19816, 792, 20608, 299, 20907, 54, 20961, 15, 20976, 10, 20986, 1438, 22424, 782, 23206, 299, 23505, 59, 23564, 15, 23579, 10, 23589, 598, 24187, 351, 24538, 143, 24681, 28, 24709, 10, 24719, 7];
const VENUS: Offsets = [40163, 548, 40711, 338, 41049, 99, 41148, 5, 41153, 4, 41157, 3, 41160, 565, 41725, 325, 42050, 99, 42149, 5, 42154, 4, 42158, 3, 42161, 190, 42351, 108, 42459, 45, 42504, 10, 42514, 3, 42517, 3];
const EARTH: Offsets = [0, 843, 843, 491, 1334, 204, 1538, 18, 1556, 15, 1571, 6, 1577, 854, 2431, 496, 2927, 202, 3129, 17, 3146, 15, 3161, 6, 3167, 178, 3345, 120, 3465, 53, 3518, 12, 3530, 6, 3536, 2];
const MARS: Offsets = [11294, 1584, 12878, 956, 13834, 387, 14221, 135, 14356, 41, 14397, 21, 14418, 1612, 16030, 969, 16999, 384, 17383, 136, 17519, 44, 17563, 21, 17584, 355, 17939, 232, 18171, 122, 18293, 51, 18344, 16, 18360, 7];
const JUPITER: Offsets = [6860, 1055, 7915, 488, 8403, 255, 8658, 140, 8798, 58, 8856, 11, 8867, 1037, 9904, 499, 10403, 259, 10662, 136, 10798, 60, 10858, 11, 10869, 216, 11085, 104, 11189, 65, 11254, 27, 11281, 10, 11291, 3];
const SATURN: Offsets = [27362, 1652, 29014, 892, 29906, 481, 30387, 215, 30602, 87, 30689, 31, 30720, 1658, 32378, 917, 33295, 465, 33760, 201, 33961, 88, 34049, 32, 34081, 420, 34501, 217, 34718, 87, 34805, 44, 34849, 19, 34868, 6];
const URANUS: Offsets = [34874, 1464, 36338, 649, 36987, 249, 37236, 84, 37320, 12, 0, 0, 37332, 1447, 38779, 659, 39438, 255, 39693, 80, 39773, 12, 0, 0, 39785, 235, 40020, 98, 40118, 33, 40151, 12, 0, 0, 0, 0];
const NEPTUNE: Offsets = [24726, 772, 25498, 330, 25828, 102, 25930, 33, 25963, 7, 0, 0, 25970, 746, 26716, 325, 27041, 97, 27138, 34, 27172, 7, 0, 0, 27179, 133, 27312, 37, 27349, 11, 27360, 2, 0, 0, 0, 0];
const EMB: Offsets = [3538, 793, 4331, 478, 4809, 185, 4994, 18, 5012, 10, 5022, 6, 5028, 804, 5832, 482, 6314, 184, 6498, 17, 6515, 10, 6525, 6, 6531, 154, 6685, 113, 6798, 46, 6844, 10, 6854, 4, 6858, 2];

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Planet {
    Mercury,
    Venus,
    Earth,
    Mars,
    Jupiter,
    Saturn,
    Uranus,
    Neptune,
}

impl Planet {
    pub fn from_name(name: &str) -> Option<Planet> {
        Some(match name {
            "Mercury" => Planet::Mercury,
            "Venus" => Planet::Venus,
            "Earth" => Planet::Earth,
            "Mars" => Planet::Mars,
            "Jupiter" => Planet::Jupiter,
            "Saturn" => Planet::Saturn,
            "Uranus" => Planet::Uranus,
            "Neptune" => Planet::Neptune,
            _ => return None,
        })
    }

    fn offsets(self) -> &'static Offsets {
        match self {
            Planet::Mercury => &MERCURY,
            Planet::Venus => &VENUS,
            Planet::Earth => &EARTH,
            Planet::Mars => &MARS,
            Planet::Jupiter => &JUPITER,
            Planet::Saturn => &SATURN,
            Planet::Uranus => &URANUS,
            Planet::Neptune => &NEPTUNE,
        }
    }
}

/// Obliquity of the ecliptic at J2000, as in `Coordinates.OBLIQUITY_DEG_J2000`.
const OBLIQUITY_DEG: f64 = 23.4392808;
/// Mass of the Moon over the mass of the Earth.
const MOON_EARTH_MASS: f64 = 0.01230073677;

pub struct Vsop87 {
    data: Vec<f64>,
}

impl Vsop87 {
    pub fn load(path: &Path) -> Result<Vsop87, String> {
        let bytes = std::fs::read(path).map_err(|e| format!("{}: {e}", path.display()))?;
        let data: Vec<f64> = bytes.chunks_exact(8).map(|c| f64::from_le_bytes(c.try_into().unwrap())).collect();
        // Every table must lie inside the file.
        for o in [&MERCURY, &VENUS, &EARTH, &MARS, &JUPITER, &SATURN, &URANUS, &NEPTUNE, &EMB] {
            for pair in o.chunks_exact(2) {
                if (pair[0] + pair[1]) * 3 > data.len() {
                    return Err(format!("{}: too short for the VSOP87A tables", path.display()));
                }
            }
        }
        Ok(Vsop87 { data })
    }

    /// Ecliptic rectangular coordinates in AU. `days` are days since J2000.
    fn ecliptic_au(&self, offsets: &Offsets, days: f64) -> [f64; 3] {
        let t = days / 365_250.0; // Julian millennia
        let mut p = [0.0; 3];
        for (i, value) in p.iter_mut().enumerate() {
            for j in 0..6 {
                let start = offsets[i * 12 + j * 2] * 3;
                let count = offsets[i * 12 + j * 2 + 1];
                let sum: f64 = self.data[start..start + count * 3].chunks_exact(3).map(|c| c[0] * (c[1] + c[2] * t).cos()).sum();
                *value += sum * t.powi(j as i32);
            }
        }
        p
    }

    /// Heliocentric position of a planet in km, in the internal equatorial axes.
    pub fn planet_km(&self, planet: Planet, days: f64) -> DVec3 {
        to_internal(self.ecliptic_au(planet.offsets(), days))
    }

    /// Position of the Moon relative to the Earth in km, from the Earth and the barycentre of
    /// the Earth and the Moon.
    pub fn moon_km(&self, days: f64) -> DVec3 {
        let earth = self.ecliptic_au(&EARTH, days);
        let emb = self.ecliptic_au(&EMB, days);
        let k = 1.0 + 1.0 / MOON_EARTH_MASS;
        to_internal([(emb[0] - earth[0]) * k, (emb[1] - earth[1]) * k, (emb[2] - earth[2]) * k])
    }
}

/// From ecliptic X, Y, Z in AU to km in the internal axes: (to RA 90°, north pole, to the equinox).
pub fn to_internal(ecl: [f64; 3]) -> DVec3 {
    let (sin, cos) = OBLIQUITY_DEG.to_radians().sin_cos();
    let [x, y, z] = ecl;
    DVec3::new(y * cos - z * sin, y * sin + z * cos, x) * AU_KM
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::astro;

    fn load() -> Option<Vsop87> {
        let home = std::env::var("HOME").ok()?;
        let path = Path::new(&home).join(".gaiasky/data/default-data/vsop87/vsop87a.bin");
        // The data is not part of the repository. Without it there is nothing to check.
        path.exists().then(|| Vsop87::load(&path).unwrap())
    }

    #[test]
    fn earth_agrees_with_the_almanac_formula() {
        let Some(v) = load() else { return };
        for date in ["2000-01-01T12:00:00", "2026-03-20T14:46:00", "2026-11-14T12:00:00", "1987-07-04"] {
            let days = astro::parse_time(date).unwrap();
            let a = v.planet_km(Planet::Earth, days);
            let b = astro::earth_position_km(days);
            let angle = a.angle_between(b).to_degrees();
            // VSOP87A is in the frame of J2000 and the formula in the frame of the date. The two
            // drift apart by the precession, 50.3" a year.
            let precession = (days / 365.25).abs() * 50.3 / 3600.0;
            assert!((angle - precession).abs() < 0.02, "{date}: {angle}°, precession {precession}°");
            assert!((a.length() - b.length()).abs() / AU_KM < 2e-4, "{date}");
        }
    }

    #[test]
    fn distances_from_the_sun() {
        let Some(v) = load() else { return };
        let days = astro::parse_time("2026-11-14").unwrap();
        // Perihelion and aphelion of each planet in AU.
        let ranges = [
            (Planet::Mercury, 0.307, 0.467),
            (Planet::Venus, 0.718, 0.729),
            (Planet::Earth, 0.983, 1.017),
            (Planet::Mars, 1.381, 1.666),
            (Planet::Jupiter, 4.95, 5.46),
            (Planet::Saturn, 9.0, 10.1),
            (Planet::Uranus, 18.2, 20.1),
            (Planet::Neptune, 29.7, 30.4),
        ];
        for (p, min, max) in ranges {
            let au = v.planet_km(p, days).length() / AU_KM;
            assert!(au > min && au < max, "{p:?}: {au} AU");
        }
    }

    #[test]
    fn moon_distance_and_full_moon() {
        let Some(v) = load() else { return };
        for date in ["2026-01-10", "2026-06-10", "2026-11-14"] {
            let d = v.moon_km(astro::parse_time(date).unwrap()).length();
            assert!((356_000.0..407_000.0).contains(&d), "{date}: {d} km");
        }
        // Total lunar eclipse of 3 March 2026, greatest at 11:34 UTC: the Moon is opposite the Sun.
        let days = astro::parse_time("2026-03-03T11:34:00").unwrap();
        let earth = v.planet_km(Planet::Earth, days);
        let angle = v.moon_km(days).angle_between(earth).to_degrees();
        assert!(angle < 0.5, "{angle}°");
    }
}
