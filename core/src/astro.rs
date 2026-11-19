//! Time, the Earth's position and its rotation.
//!
//! Axes are the ones Gaia Sky uses internally: equatorial, with +Y to the north celestial pole,
//! +Z to the vernal equinox and +X to right ascension 90°.

use std::time::{SystemTime, UNIX_EPOCH};

use glam::{DMat3, DQuat, DVec3};

pub const AU_KM: f64 = 149_597_870.7;
pub const PC_KM: f64 = 3.08567758149137e13;

/// Unix time of the J2000 epoch, 2000-01-01T12:00:00. The 64 s between TT and UTC are ignored.
const J2000_UNIX: f64 = 946_728_000.0;
pub const OBLIQUITY_DEG: f64 = 23.4392808;

pub fn days_since_j2000_now() -> f64 {
    let unix = SystemTime::now().duration_since(UNIX_EPOCH).map(|d| d.as_secs_f64()).unwrap_or(J2000_UNIX);
    (unix - J2000_UNIX) / 86_400.0
}

/// Parses `YYYY-MM-DDTHH:MM:SS` (UTC, the time part is optional) into days since J2000.
pub fn parse_time(s: &str) -> Option<f64> {
    let s = s.trim_end_matches('Z');
    let (date, time) = s.split_once(['T', ' ']).unwrap_or((s, "00:00:00"));
    let mut d = date.split('-').map(|p| p.parse::<i64>());
    let (y, m, day) = (d.next()?.ok()?, d.next()?.ok()?, d.next()?.ok()?);
    let mut t = time.split(':').map(|p| p.parse::<f64>());
    let (h, min, sec) = (t.next()?.ok()?, t.next().unwrap_or(Ok(0.0)).ok()?, t.next().unwrap_or(Ok(0.0)).ok()?);
    if d.next().is_some() || !(1..=12).contains(&m) || !(1..=31).contains(&day) || !(0.0..24.0).contains(&h) || !(0.0..60.0).contains(&min) || !(0.0..61.0).contains(&sec) {
        return None;
    }
    // Days from the civil date, proleptic Gregorian calendar.
    let (y2, m2) = if m <= 2 { (y - 1, m + 12) } else { (y, m) };
    let days = 365 * y2 + y2.div_euclid(4) - y2.div_euclid(100) + y2.div_euclid(400) + (153 * (m2 - 3) + 2) / 5 + day - 730_426;
    Some(days as f64 - 0.5 + (h * 3600.0 + min * 60.0 + sec) / 86_400.0)
}

/// Heliocentric position of the Earth in km, referred to the equinox of the date.
///
/// This is the low-precision solar position of the Astronomical Almanac, good to about 0.01° of
/// longitude. The viewer uses VSOP87; this is an independent check of it in the tests.
#[cfg(test)]
pub fn earth_position_km(days: f64) -> DVec3 {
    let g = (357.528 + 0.985_600_3 * days).to_radians();
    let sun_lon = 280.460 + 0.985_647_4 * days + 1.915 * g.sin() + 0.020 * (2.0 * g).sin();
    let r = (1.000_14 - 0.016_71 * g.cos() - 0.000_14 * (2.0 * g).cos()) * AU_KM;
    let lon = (sun_lon + 180.0).to_radians();
    let (xe, ye) = (r * lon.cos(), r * lon.sin());
    let eps = OBLIQUITY_DEG.to_radians();
    // Ecliptic to equatorial, then to the internal axes: (to RA 90°, north, to the equinox).
    DVec3::new(ye * eps.cos(), ye * eps.sin(), xe)
}

/// Rotation of the Earth around its axis: the Greenwich mean sidereal time, in radians.
/// The prime meridian of the textures is the model's +Z axis.
pub fn earth_rotation(days: f64) -> f64 {
    (280.460_618_37 + 360.985_647_366_29 * days).rem_euclid(360.0).to_radians()
}

/// Orientation of a camera that looks along `dir`, with `up` roughly up.
pub fn look_at(dir: DVec3, up: DVec3) -> DQuat {
    let back = -dir.normalize();
    let right = up.cross(back).normalize();
    DQuat::from_mat3(&DMat3::from_cols(right, back.cross(right), back)).normalize()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_times() {
        assert_eq!(parse_time("2000-01-01T12:00:00"), Some(0.0));
        assert_eq!(parse_time("2000-01-02"), Some(0.5));
        assert_eq!(parse_time("1999-12-31T12:00:00Z"), Some(-1.0));
        assert_eq!(parse_time("2026-11-14T12:00:00"), Some(9814.0));
        assert_eq!(parse_time("2024-02-29T12:00"), Some(8825.0));
        assert_eq!(parse_time("2026-13-01"), None);
        assert_eq!(parse_time("yesterday"), None);
    }

    #[test]
    fn earth_is_where_it_should_be() {
        // Around 4 January the Earth is at perihelion, about 0.983 AU from the Sun.
        let p = earth_position_km(parse_time("2026-01-04").unwrap());
        assert!((p.length() / AU_KM - 0.9833).abs() < 5e-4, "{}", p.length() / AU_KM);
        // At the March equinox the Sun is at the vernal equinox as seen from the Earth, so the
        // Earth is on the -Z side and close to the equatorial plane.
        let p = earth_position_km(parse_time("2026-03-20T14:46:00").unwrap()).normalize();
        assert!(p.z < -0.9999 && p.y.abs() < 1e-3, "{p:?}");
        // At the June solstice the Sun is 23.44° north, so the Earth is that far south.
        let p = earth_position_km(parse_time("2026-06-21T08:24:00").unwrap()).normalize();
        assert!((p.y.asin().to_degrees() + 23.44).abs() < 0.01, "{}", p.y.asin().to_degrees());
    }

    #[test]
    fn noon_is_over_greenwich() {
        // At 12:00 UTC the Sun is over the prime meridian, give or take the equation of time
        // (under 16 minutes, that is 4°).
        for date in ["2026-02-11T12:00:00", "2026-06-15T12:00:00", "2026-11-03T12:00:00"] {
            let days = parse_time(date).unwrap();
            let to_sun = -earth_position_km(days).normalize();
            let theta = earth_rotation(days);
            let greenwich = DVec3::new(theta.sin(), 0.0, theta.cos());
            let sun_ra = to_sun.x.atan2(to_sun.z);
            let diff = (sun_ra - theta + std::f64::consts::PI).rem_euclid(std::f64::consts::TAU) - std::f64::consts::PI;
            assert!(diff.to_degrees().abs() < 4.2, "{date}: {}", diff.to_degrees());
            assert!(greenwich.dot(to_sun) > 0.9);
        }
    }

    #[test]
    fn look_at_looks() {
        let q = look_at(DVec3::new(1.0, 2.0, 3.0), DVec3::Y);
        let f = q * DVec3::NEG_Z;
        assert!((f - DVec3::new(1.0, 2.0, 3.0).normalize()).length() < 1e-12);
        assert!((q * DVec3::Y).y > 0.0);
    }
}
