//! What a made-up star is like: its class, mass, size and temperature, from its brightness and
//! its kind. The numbers follow the main sequence of real stars; nothing here is stored.

use crate::procedural::{Kind, Star, CLASSES};

/// A star in words and numbers. Mass and radius are in units of the Sun's.
#[derive(Clone, Debug, PartialEq)]
pub struct Description {
    /// The spectral class, as astronomers write it: "G2 V" is a star like the Sun.
    pub class: String,
    pub mass: f64,
    pub radius: f64,
    /// Surface temperature in kelvins.
    pub temperature: f64,
    /// Visible light given off, in units of the Sun's.
    pub luminosity: f64,
}

impl Description {
    /// All the light and heat given off, in units of the Sun's: what warms the planets. Cool
    /// stars give most of it as heat that the eye does not see.
    pub fn power(&self) -> f64 {
        self.radius * self.radius * (self.temperature / 5_770.0).powi(4)
    }
}

/// The main sequence: absolute magnitude, mass, radius, temperature, from the hottest and
/// largest stars (O5: eleven times the Sun across) to the smallest (M9.5: a tenth of the Sun).
/// Values in between are interpolated. The rows are those of the table of Pecaut and Mamajek
/// (version 2022.04.16, see map/SOURCES.md).
const MAIN_SEQUENCE: [(f64, f64, f64, f64); 27] = [
    (-5.35, 43.0, 11.45, 41400.0), // O5V
    (-4.2, 20.2, 7.72, 33300.0), // O9V
    (-3.9, 17.7, 7.16, 31400.0), // B0V
    (-1.8, 7.3, 4.06, 20600.0), // B2V
    (-0.85, 4.7, 3.36, 15700.0), // B5V
    (0.0, 3.38, 2.86, 12300.0), // B8V
    (0.99, 2.18, 2.193, 9700.0), // A0V
    (2.01, 1.88, 1.785, 8100.0), // A5V
    (2.57, 1.61, 1.728, 7220.0), // F0V
    (3.37, 1.33, 1.473, 6550.0), // F5V
    (4.48, 1.06, 1.1, 5930.0), // G0V
    (4.8, 1.0, 1.012, 5770.0), // G2V
    (4.98, 0.98, 0.977, 5660.0), // G5V
    (5.78, 0.88, 0.813, 5270.0), // K0V
    (6.07, 0.82, 0.783, 5100.0), // K2V
    (7.28, 0.7, 0.701, 4440.0), // K5V
    (8.16, 0.64, 0.63, 4100.0), // K7V
    (8.8, 0.57, 0.588, 3850.0), // M0V
    (10.21, 0.44, 0.446, 3560.0), // M2V
    (11.15, 0.37, 0.361, 3430.0), // M3V
    (12.61, 0.23, 0.274, 3210.0), // M4V
    (14.15, 0.162, 0.196, 3060.0), // M5V
    (16.32, 0.102, 0.137, 2810.0), // M6V
    (17.7, 0.09, 0.12, 2680.0), // M7V
    (18.6, 0.085, 0.114, 2570.0), // M8V
    (19.4, 0.079, 0.102, 2380.0), // M9V
    (19.75, 0.078, 0.101, 2350.0), // M9.5V
];

fn along(table: &[(f64, f64, f64, f64)], m: f64) -> (f64, f64, f64) {
    let m = m.clamp(table[0].0, table[table.len() - 1].0);
    let i = table.windows(2).position(|w| m <= w[1].0).unwrap_or(table.len() - 2);
    let (a, b) = (table[i], table[i + 1]);
    let t = (m - a.0) / (b.0 - a.0);
    // Masses and sizes span decades: in between, their logarithms go evenly.
    let mix = |x: f64, y: f64| (x.ln() + (y.ln() - x.ln()) * t).exp();
    (mix(a.1, b.1), mix(a.2, b.2), mix(a.3, b.3))
}

/// The absolute magnitude of a main-sequence star of a mass, in solar masses.
pub fn abs_mag_of_mass(mass: f64) -> f64 {
    let t = &MAIN_SEQUENCE;
    let mass = mass.clamp(t[t.len() - 1].1, t[0].1);
    let i = t.windows(2).position(|w| mass >= w[1].1).unwrap_or(t.len() - 2);
    let (a, b) = (t[i], t[i + 1]);
    a.0 + (b.0 - a.0) * (mass.ln() - a.1.ln()) / (b.1.ln() - a.1.ln())
}

/// The class letter and its tenth for a main-sequence star of an absolute magnitude.
fn class_of(m: f64) -> (char, u32) {
    let i = CLASSES.iter().position(|c| m < c.1).unwrap_or(CLASSES.len() - 1);
    let from = if i == 0 { -7.0 } else { CLASSES[i - 1].1 };
    let tenth = ((m - from) / (CLASSES[i].1 - from) * 10.0).clamp(0.0, 9.0) as u32;
    (CLASSES[i].0, tenth)
}

/// Describes a star. `seed` is the star's own number, for what its brightness does not settle.
pub fn describe(star: &Star, seed: u64) -> Description {
    let m = star.abs_mag as f64;
    let luminosity = 10f64.powf((4.83 - m) / 2.5);
    // Two numbers in [0, 1) of the star's own.
    let draw = |k: u64| ((seed ^ k.wrapping_mul(0xD6E8_FEB8_6659_FD93)).wrapping_mul(0x9E37_79B9_7F4A_7C15) >> 11) as f64 / (1u64 << 53) as f64;
    match star.kind {
        Kind::MainSequence => {
            let (mass, radius, temperature) = along(&MAIN_SEQUENCE, m);
            let (letter, tenth) = class_of(m);
            Description { class: format!("{letter}{tenth} V"), mass, radius, temperature, luminosity }
        }
        Kind::Giant => {
            // An old star of one to three solar masses, swollen and cooled: orange, class G to K.
            let temperature = 4_000.0 + 1_200.0 * draw(1);
            let radius = (luminosity * 1.4).sqrt() * (5_770.0f64 / temperature).powi(2);
            let class = format!("{}{} III", if temperature > 4_900.0 { 'G' } else { 'K' }, (draw(2) * 10.0) as u32);
            Description { class, mass: 1.0 + 2.0 * draw(3), radius, temperature, luminosity }
        }
        Kind::WhiteDwarf => {
            // The core of a dead star, the size of the Earth; it cools and fades with age.
            let fade = ((m - 10.0) / 5.5).clamp(0.0, 1.0);
            Description { class: "white dwarf".into(), mass: 0.5 + 0.3 * draw(1), radius: 0.012, temperature: 25_000.0 * (5_000.0f64 / 25_000.0).powf(fade), luminosity }
        }
    }
}

impl Description {
    pub fn line(&self) -> String {
        format!("{}, {:.2} solar masses, {:.2} solar radii, {:.0} K, {} times the Sun's light", self.class, self.mass, self.radius, self.temperature, if self.luminosity < 0.01 { format!("{:.1e}", self.luminosity) } else { format!("{:.2}", self.luminosity) })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::procedural::StarId;
    use glam::DVec3;

    fn star(abs_mag: f32, kind: Kind) -> Star {
        Star { position: DVec3::ZERO, abs_mag, color: [1.0; 3], kind, id: StarId { level: 0, cell: [0; 3], index: 0 } }
    }

    #[test]
    fn the_sun_is_a_g2_star() {
        let sun = describe(&star(4.83, Kind::MainSequence), 1);
        assert_eq!(sun.class, "G2 V");
        assert!((sun.mass - 1.0).abs() < 0.02 && (sun.radius - 1.0).abs() < 0.02 && (sun.luminosity - 1.0).abs() < 1e-6);
    }

    #[test]
    fn brighter_main_sequence_stars_are_heavier_and_hotter() {
        let mut last = describe(&star(19.0, Kind::MainSequence), 1);
        assert!(last.class.starts_with('M') && last.mass < 0.1);
        for m in [12.0, 9.0, 6.0, 3.0, 0.0, -3.0, -5.5] {
            let d = describe(&star(m, Kind::MainSequence), 1);
            assert!(d.mass > last.mass && d.temperature > last.temperature && d.radius > last.radius, "{m}: {d:?}");
            last = d;
        }
        assert!(last.class.starts_with('O') && last.mass > 30.0);
    }

    #[test]
    fn mass_and_brightness_go_both_ways() {
        for m in [-5.0, 0.0, 4.83, 9.0, 15.0, 19.0] {
            let mass = describe(&star(m, Kind::MainSequence), 1).mass;
            assert!((abs_mag_of_mass(mass) - m as f64).abs() < 1e-4, "{m}");
        }
    }

    #[test]
    fn giants_are_large_and_white_dwarfs_small() {
        let giant = describe(&star(0.5, Kind::Giant), 42);
        assert!(giant.radius > 8.0 && giant.class.ends_with("III"), "{giant:?}");
        let dwarf = describe(&star(13.0, Kind::WhiteDwarf), 42);
        assert!(dwarf.radius < 0.02 && (0.5..0.8).contains(&dwarf.mass));
        // The same star is described the same way every time.
        assert_eq!(giant, describe(&star(0.5, Kind::Giant), 42));
    }
}
