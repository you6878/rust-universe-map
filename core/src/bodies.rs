//! The Sun, the planets and the Moon, as described by the JSON files of the Gaia Sky data.

use std::path::{Path, PathBuf};

use glam::{DMat3, DVec3};
use serde_json::Value;

use crate::astro;
use crate::vsop87::{self, Planet, Vsop87};

/// Where a body is.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Orbit {
    /// At the origin.
    Sun,
    Planet(Planet),
    /// Around the Earth.
    Moon,
    Pluto,
}

/// Rotation of a body, as in the `rotation` block of the JSON. Angles in degrees.
#[derive(Debug, Clone, Copy, Default)]
pub struct Rotation {
    /// Sidereal period in hours.
    pub period_h: f64,
    pub axial_tilt: f64,
    pub inclination: f64,
    pub ascending_node: f64,
    pub meridian_angle: f64,
}

pub struct Ring {
    /// In planet radii.
    pub inner: f64,
    pub outer: f64,
    pub texture: PathBuf,
}

pub struct Body {
    pub name: String,
    pub radius_km: f64,
    pub orbit: Orbit,
    pub rotation: Rotation,
    /// Magnitude at 1 AU from the Sun and from the observer, for drawing the body as a point.
    pub abs_mag: f64,
    pub color: [f32; 3],
    /// Equirectangular surface map. The Earth and the Sun have none: the Earth uses its cube maps.
    pub texture: Option<PathBuf>,
    pub ring: Option<Ring>,
}

fn objects(path: &Path) -> Result<Vec<Value>, String> {
    let text = std::fs::read_to_string(path).map_err(|e| format!("{}: {e}", path.display()))?;
    let v: Value = serde_json::from_str(&text).map_err(|e| format!("{}: {e}", path.display()))?;
    v.get("objects").and_then(Value::as_array).cloned().ok_or(format!("{}: no \"objects\" array", path.display()))
}

/// Looks a key up in an object. The JSON files are not consistent about the case of the keys.
fn get<'a>(v: &'a Value, key: &str) -> Option<&'a Value> {
    v.as_object()?.iter().find(|(k, _)| k.eq_ignore_ascii_case(key)).map(|(_, v)| v)
}

fn num(v: &Value, key: &str) -> Option<f64> {
    get(v, key)?.as_f64()
}

/// Turns a path of the data files into a real one. `$data` is the data folder and `*` stands for
/// the quality suffix; the best quality that is installed wins.
fn resolve(data: &Path, spec: &str) -> Option<PathBuf> {
    let rel = spec.strip_prefix("$data/")?;
    ["-high", "-med", "-low", ""].iter().map(|q| data.join(rel.replace('*', q))).find(|p| p.exists())
}

fn body(v: &Value, data: &Path, orbit: Orbit) -> Option<Body> {
    let name = get(v, "name")?.as_str()?.to_owned();
    let rot = get(v, "rotation");
    let r = |key| rot.and_then(|r| num(r, key)).unwrap_or(0.0);
    let model = get(v, "model");
    let material = model.and_then(|m| get(m, "material"));
    let tex = |key| material.and_then(|m| get(m, key)).and_then(Value::as_str).and_then(|s| resolve(data, s));
    let params = model.and_then(|m| get(m, "params"));
    // Ring radii are in units of the diameter.
    let ring = tex("ring").and_then(|texture| {
        Some(Ring { inner: num(params?, "innerradius")? * 2.0, outer: num(params?, "outerradius")? * 2.0, texture })
    });
    let color = get(v, "color").and_then(Value::as_array).map(|c| [0, 1, 2].map(|i| c.get(i).and_then(Value::as_f64).unwrap_or(1.0) as f32));
    Some(Body {
        radius_km: num(v, "size")?,
        orbit,
        rotation: Rotation {
            period_h: r("period"),
            axial_tilt: r("axialTilt"),
            inclination: r("inclination"),
            ascending_node: r("ascendingNode"),
            meridian_angle: r("meridianAngle"),
        },
        abs_mag: num(v, "absMag").unwrap_or(0.0),
        color: color.unwrap_or([1.0; 3]),
        texture: tex("diffuse").or_else(|| tex("base")),
        ring,
        name,
    })
}

/// Reads the Sun, the planets and the Moon from the two files of the map. `data` is the folder
/// that `$data` stands for in the texture paths.
pub fn load(planets: &Path, moons: &Path, data: &Path) -> Result<Vec<Body>, String> {
    let mut bodies = vec![Body {
        name: "Sun".into(),
        radius_km: 695_700.0,
        orbit: Orbit::Sun,
        rotation: Rotation::default(),
        abs_mag: 0.0,
        color: [1.0, 0.92, 0.82],
        texture: None,
        ring: None,
    }];
    for v in objects(planets)? {
        let planet = get(&v, "name").and_then(Value::as_str).and_then(Planet::from_name);
        if let Some(b) = planet.and_then(|p| body(&v, data, Orbit::Planet(p))) {
            bodies.push(b);
        }
    }
    for v in objects(moons)? {
        if get(&v, "name").and_then(Value::as_str) == Some("Moon") {
            // The Moon has no magnitude in the file. 0.2 is its magnitude at 1 AU.
            bodies.extend(body(&v, data, Orbit::Moon).map(|b| Body { abs_mag: 0.2, ..b }));
        }
    }
    // Pluto goes last, so that the planets keep the numbers 1 to 8 and the Moon the 9.
    for v in objects(planets)? {
        if get(&v, "name").and_then(Value::as_str) == Some("Pluto") {
            bodies.extend(body(&v, data, Orbit::Pluto));
        }
    }
    if bodies.len() < 10 {
        return Err(format!("{}: expected the eight planets and the Moon, found {} bodies", planets.display(), bodies.len() - 1));
    }
    Ok(bodies)
}

/// Heliocentric position of Pluto in km, as `AstroUtils.plutoEclipticCoordinates`: a fit by
/// Paul Schlyter, good to about a minute of arc between 1800 and 2100.
fn pluto_km(d: f64) -> DVec3 {
    let s = (50.03 + 0.033459652 * d).to_radians();
    let p = (238.95 + 0.003968789 * d).to_radians();
    let (sin, cos) = (|k: f64| (k * p).sin(), |k: f64| (k * p).cos());
    let lon = 238.9508 + 0.00400703 * d - 19.799 * sin(1.0) + 19.848 * cos(1.0) + 0.897 * sin(2.0) - 4.956 * cos(2.0) + 0.610 * sin(3.0) + 1.211 * cos(3.0)
        - 0.341 * sin(4.0)
        - 0.190 * cos(4.0)
        + 0.128 * sin(5.0)
        - 0.034 * cos(5.0)
        - 0.038 * sin(6.0)
        + 0.031 * cos(6.0)
        + 0.020 * (s - p).sin()
        - 0.010 * (s - p).cos();
    let lat = -3.9082 - 5.453 * sin(1.0) - 14.975 * cos(1.0) + 3.527 * sin(2.0) + 1.673 * cos(2.0) - 1.051 * sin(3.0) + 0.328 * cos(3.0) + 0.179 * sin(4.0)
        - 0.292 * cos(4.0)
        + 0.019 * sin(5.0)
        + 0.100 * cos(5.0)
        - 0.031 * sin(6.0)
        - 0.026 * cos(6.0)
        + 0.011 * (s - p).cos();
    let r = 40.72 + 6.68 * sin(1.0) + 6.90 * cos(1.0) - 1.18 * sin(2.0) - 0.03 * cos(2.0) + 0.15 * sin(3.0) - 0.14 * cos(3.0);
    let (lon, lat) = (lon.to_radians(), lat.to_radians());
    vsop87::to_internal([r * lat.cos() * lon.cos(), r * lat.cos() * lon.sin(), r * lat.sin()])
}

impl Body {
    /// Heliocentric position in km.
    pub fn position_km(&self, vsop: &Vsop87, days: f64) -> DVec3 {
        match self.orbit {
            Orbit::Sun => DVec3::ZERO,
            Orbit::Planet(p) => vsop.planet_km(p, days),
            Orbit::Moon => vsop.planet_km(Planet::Earth, days) + vsop.moon_km(days),
            Orbit::Pluto => pluto_km(days),
        }
    }

    /// Orientation of the body: from its own axes (+Y is the north pole) to the internal axes.
    ///
    /// As in `ModelUpdater`: node, ecliptic to equatorial, inclination plus tilt, spin.
    pub fn orientation(&self, days: f64) -> DMat3 {
        if self.orbit == Orbit::Planet(Planet::Earth) {
            // The Earth's textures are cube maps with the prime meridian on +Z.
            return DMat3::from_rotation_y(astro::earth_rotation(days));
        }
        let r = &self.rotation;
        let spin = if r.period_h != 0.0 { r.meridian_angle + 360.0 / r.period_h * days * 24.0 } else { 0.0 };
        DMat3::from_rotation_y(r.ascending_node.to_radians())
            * DMat3::from_rotation_z(astro::OBLIQUITY_DEG.to_radians())
            * DMat3::from_rotation_z((r.inclination + r.axial_tilt).to_radians())
            * DMat3::from_rotation_y(spin.rem_euclid(360.0).to_radians())
    }

    /// Apparent magnitude seen from `observer`, both positions heliocentric in km.
    pub fn apparent_mag(&self, position: DVec3, observer: DVec3) -> f64 {
        if self.orbit == Orbit::Sun {
            // The Sun is in the star list already.
            return f64::INFINITY;
        }
        let r = (position.length() / astro::AU_KM).max(1e-6);
        let d = ((position - observer).length() / astro::AU_KM).max(1e-9);
        // How much of the lit half faces the observer, for a ball that scatters evenly.
        let alpha = (-position).angle_between(observer - position);
        let phase = ((alpha.sin() + (std::f64::consts::PI - alpha) * alpha.cos()) / std::f64::consts::PI).max(1e-4);
        self.abs_mag + 5.0 * (r * d).log10() - 2.5 * phase.log10()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn setup() -> Option<(Vec<Body>, Vsop87)> {
        let data = Path::new(&std::env::var("HOME").ok()?).join(".gaiasky/data");
        // The data is not part of the repository. Without it there is nothing to check.
        data.join("default-data/planets.json").exists().then(|| (load(&data.join("default-data/planets.json"), &data.join("default-data/moons.json"), &data).unwrap(), Vsop87::load(&data.join("default-data/vsop87/vsop87a.bin")).unwrap()))
    }

    #[test]
    fn reads_the_bodies() {
        let Some((bodies, _)) = setup() else { return };
        let names: Vec<&str> = bodies.iter().map(|b| b.name.as_str()).collect();
        assert_eq!(names, ["Sun", "Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune", "Moon", "Pluto"]);
        let saturn = &bodies[6];
        assert_eq!(saturn.radius_km, 58232.0);
        assert!(saturn.texture.is_some());
        let ring = saturn.ring.as_ref().unwrap();
        assert!((ring.inner - 1.279).abs() < 1e-3 && (ring.outer - 2.408).abs() < 1e-3);
        assert!(bodies[3].texture.is_none(), "the Earth uses cube maps");
        assert!(bodies[9].texture.is_some(), "the Moon's map is under \"base\"");
        assert_eq!(bodies[9].rotation.axial_tilt, -6.68, "lower-case keys are read too");
    }

    #[test]
    fn poles_point_the_right_way() {
        let Some((bodies, _)) = setup() else { return };
        let pole = |i: usize| bodies[i].orientation(1234.5) * DVec3::Y;
        // The Earth's pole is the celestial pole.
        assert!(pole(3).y > 0.999_999);
        // The poles of the others are tilted from the pole of the ecliptic by the inclination plus
        // the axial tilt of the file: 24.24° for Saturn, and Uranus lies on its side.
        let eps = astro::OBLIQUITY_DEG.to_radians();
        let ecliptic_pole = DVec3::new(-eps.sin(), eps.cos(), 0.0);
        let tilt = |i: usize| pole(i).angle_between(ecliptic_pole).to_degrees();
        assert!((tilt(6) - 24.24476).abs() < 1e-6, "{}", tilt(6));
        assert!((tilt(7) - 96.997444).abs() < 1e-6, "{}", tilt(7));
    }

    #[test]
    fn pluto_is_where_it_should_be() {
        // In November 2026 Pluto is in Capricornus, at right ascension 20h29m and declination
        // -23.3°, 35.6 AU from the Sun.
        let p = pluto_km(astro::parse_time("2026-11-14").unwrap());
        assert!((p.length() / astro::AU_KM - 35.6).abs() < 0.5, "{} AU", p.length() / astro::AU_KM);
        let ra = p.x.atan2(p.z).to_degrees().rem_euclid(360.0);
        let dec = (p.y / p.length()).asin().to_degrees();
        assert!((ra - 307.0).abs() < 3.0 && (dec + 23.0).abs() < 2.0, "ra {ra} dec {dec}");
    }

    #[test]
    fn magnitudes_are_plausible() {
        let Some((bodies, vsop)) = setup() else { return };
        let days = astro::parse_time("2026-11-14").unwrap();
        let earth = bodies[3].position_km(&vsop, days);
        let mag = |i: usize| bodies[i].apparent_mag(bodies[i].position_km(&vsop, days), earth);
        assert!((-3.0..-1.5).contains(&mag(5)), "Jupiter {}", mag(5));
        // Venus is a thin crescent near the Earth on this date.
        assert!((-5.0..-3.0).contains(&mag(2)), "Venus {}", mag(2));
        // The Moon is a thin crescent on this date, and full on 3 March 2026.
        assert!((-10.5..-9.0).contains(&mag(9)), "crescent Moon {}", mag(9));
        let full = astro::parse_time("2026-03-03T11:34:00").unwrap();
        let m = bodies[9].apparent_mag(bodies[9].position_km(&vsop, full), bodies[3].position_km(&vsop, full));
        assert!((-13.2..-12.2).contains(&m), "full Moon {m}");
    }
}
