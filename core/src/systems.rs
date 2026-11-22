//! Star systems: whether a made-up star is alone or has companions, and how they go around each
//! other. Everything follows from the star's own number, so the same star always has the same
//! system, and nothing is stored. The statistics are those of real stars (map/SOURCES.md).
//!
//! A system is a tree: a star, or two parts that orbit their common centre of mass, each of them
//! a star or a pair again. A triple is a close pair with a third star far out; orbits are plain
//! Kepler ellipses, and a system is only made if it would hold together.

use glam::{DQuat, DVec3};

use crate::procedural::{Kind, Rng, Star};
use crate::stars::{self, Description};

/// Stars a system has at most.
const MOST_STARS: usize = 7;
/// No companion is farther than this from its partner, in AU: about a tenth of a parsec.
const WIDEST_AU: f64 = 20_000.0;
const SOLAR_RADIUS_AU: f64 = 0.004_65;
/// The lightest thing that is a star, in solar masses.
const LIGHTEST_STAR: f64 = 0.078;

/// How often stars of a mass have companions: for primaries of at least `min_mass` solar masses,
/// the shares of systems with two, three, and four or more stars. The rest are single.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Multiplicity {
    pub min_mass: f64,
    pub double: f64,
    pub triple: f64,
    pub more: f64,
}

/// The table of real stars, lightest first: M dwarfs (Winters et al. 2019), sun-like stars
/// (Raghavan et al. 2010), then A, B and O stars (Offner et al. 2023: their shares of multiples
/// and of triples or more; how the latter split into three and more is a choice made here).
pub const MULTIPLICITY: [Multiplicity; 6] = [
    Multiplicity { min_mass: 0.0, double: 0.20, triple: 0.03, more: 0.003 },
    Multiplicity { min_mass: 0.57, double: 0.33, triple: 0.08, more: 0.03 },
    Multiplicity { min_mass: 1.61, double: 0.43, triple: 0.18, more: 0.07 },
    Multiplicity { min_mass: 3.0, double: 0.45, triple: 0.25, more: 0.11 },
    Multiplicity { min_mass: 8.0, double: 0.36, triple: 0.35, more: 0.22 },
    Multiplicity { min_mass: 17.7, double: 0.28, triple: 0.40, more: 0.28 },
];

/// A Kepler orbit of one part of a pair around the other.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Orbit {
    pub semi_major_au: f64,
    pub eccentricity: f64,
    pub period_years: f64,
    /// How the ellipse lies in space.
    plane: DQuat,
    /// Where on the orbit the pair is at the epoch (J2000): the mean anomaly, in radians.
    phase: f64,
}

impl Orbit {
    /// From the first part of the pair to the second, in AU, at a time in years since J2000.
    pub fn offset_au(&self, years: f64) -> DVec3 {
        let mean = (self.phase + std::f64::consts::TAU * years / self.period_years).rem_euclid(std::f64::consts::TAU);
        // Kepler's equation, by Newton's method.
        let e = self.eccentricity;
        let mut anomaly = if e < 0.8 { mean } else { std::f64::consts::PI };
        for _ in 0..30 {
            let step = (anomaly - e * anomaly.sin() - mean) / (1.0 - e * anomaly.cos());
            anomaly -= step;
            if step.abs() < 1e-13 {
                break;
            }
        }
        let a = self.semi_major_au;
        self.plane * DVec3::new(a * (anomaly.cos() - e), a * (1.0 - e * e).sqrt() * anomaly.sin(), 0.0)
    }
}

/// One star of a system.
#[derive(Clone, Debug, PartialEq)]
pub struct Member {
    pub about: Description,
    pub abs_mag: f64,
    pub kind: Kind,
}

#[derive(Clone, Debug, PartialEq)]
pub enum Part {
    Star(Member),
    /// Two parts and the orbit of the second around the first.
    Pair(Box<Part>, Box<Part>, Orbit),
}

impl Part {
    pub fn mass(&self) -> f64 {
        match self {
            Part::Star(m) => m.about.mass,
            Part::Pair(a, b, _) => a.mass() + b.mass(),
        }
    }

    pub fn stars(&self) -> usize {
        match self {
            Part::Star(_) => 1,
            Part::Pair(a, b, _) => a.stars() + b.stars(),
        }
    }

    /// Widest reach of the part from its centre of mass, in AU.
    fn reach_au(&self) -> f64 {
        match self {
            Part::Star(_) => 0.0,
            Part::Pair(a, b, o) => o.semi_major_au * (1.0 + o.eccentricity) + a.reach_au().max(b.reach_au()),
        }
    }

    fn place(&self, centre: DVec3, years: f64, out: &mut Vec<(DVec3, Member)>) {
        match self {
            Part::Star(m) => out.push((centre, m.clone())),
            Part::Pair(a, b, orbit) => {
                let r = orbit.offset_au(years);
                let (ma, mb) = (a.mass(), b.mass());
                a.place(centre - r * (mb / (ma + mb)), years, out);
                b.place(centre + r * (ma / (ma + mb)), years, out);
            }
        }
    }

    fn orbits<'a>(&'a self, out: &mut Vec<&'a Orbit>) {
        if let Part::Pair(a, b, o) = self {
            a.orbits(out);
            b.orbits(out);
            out.push(o);
        }
    }
}

/// A star with its companions.
#[derive(Clone, Debug, PartialEq)]
pub struct System {
    pub root: Part,
}

impl System {
    pub fn stars(&self) -> usize {
        self.root.stars()
    }

    /// Where the stars are, in AU from the centre of mass of the system, at a time in years since
    /// J2000. The first is the star the system was made from.
    pub fn places_au(&self, years: f64) -> Vec<(DVec3, Member)> {
        let mut out = Vec::new();
        self.root.place(DVec3::ZERO, years, &mut out);
        out
    }

    /// The orbits, innermost first.
    pub fn orbits(&self) -> Vec<&Orbit> {
        let mut out = Vec::new();
        self.root.orbits(&mut out);
        out
    }

    /// The absolute magnitude of all the stars together, as seen from too far to tell them apart.
    /// (From afar the viewer draws a system as its first star alone; this is for the tests.)
    #[cfg_attr(not(test), allow(dead_code))]
    pub fn abs_mag(&self) -> f64 {
        let light: f64 = self.places_au(0.0).iter().map(|(_, m)| 10f64.powf(-0.4 * m.abs_mag)).sum();
        -2.5 * light.log10()
    }

    /// One line for each star and each orbit.
    pub fn lines(&self) -> Vec<String> {
        let mut lines: Vec<String> = self.places_au(0.0).iter().map(|(_, m)| format!("star: {}", m.about.line())).collect();
        for o in self.orbits() {
            let period = if o.period_years < 1.0 { format!("{:.1} days", o.period_years * 365.25) } else { format!("{:.1} years", o.period_years) };
            lines.push(format!("orbit: {:.3} AU, eccentricity {:.2}, once in {period}", o.semi_major_au, o.eccentricity));
        }
        lines
    }
}

/// How far apart a third body must be for a pair to keep it (Mardling and Aarseth 2001): the
/// ratio of the outer semi-major axis to the inner one, for an outer orbit of eccentricity `e`
/// and an outer body of `q` times the mass of the inner pair.
pub fn stable_ratio(e: f64, q: f64) -> f64 {
    2.8 / (1.0 - e) * ((1.0 + q) * (1.0 + e) / (1.0 - e).sqrt()).powf(0.4)
}

fn normal(rng: &mut Rng) -> f64 {
    let (u, v) = (rng.next().max(1e-300), rng.next());
    (-2.0 * u.ln()).sqrt() * (std::f64::consts::TAU * v).cos()
}

/// The number of stars of a system whose first star has a mass.
fn count(table: &[Multiplicity], mass: f64, rng: &mut Rng) -> usize {
    let row = table.iter().rev().find(|r| mass >= r.min_mass).unwrap_or(&table[0]);
    let u = rng.next();
    if u < row.more {
        // Four, and each further star is rarer by a factor.
        let mut n = 4;
        while n < MOST_STARS && rng.next() < 0.25 {
            n += 1;
        }
        n
    } else if u < row.more + row.triple {
        3
    } else if u < row.more + row.triple + row.double {
        2
    } else {
        1
    }
}

/// A separation in AU for a companion of a star of a mass, and the sum of the two masses.
fn separation(primary: f64, total: f64, rng: &mut Rng) -> f64 {
    let au = if primary < 0.57 {
        // M dwarfs keep their companions close: most often at 20 AU.
        10f64.powf(20f64.log10() + 1.16 * normal(rng))
    } else {
        // Sun-like and heavier: periods spread over ten decades around 300 years.
        let years = 10f64.powf(5.03 + 2.28 * normal(rng)) / 365.25;
        (years * years * total).cbrt()
    };
    au.clamp(0.005, WIDEST_AU)
}

/// An orbit of a size and shape around a mass (in solar masses), lying any way in space.
pub fn orbit(semi_major_au: f64, eccentricity: f64, mass: f64, rng: &mut Rng) -> Orbit {
    // A plane with a random pole, turned by a random angle.
    let (z, turn) = (rng.next() * 2.0 - 1.0, rng.next() * std::f64::consts::TAU);
    let pole = DVec3::new((1.0 - z * z).sqrt() * turn.cos(), (1.0 - z * z).sqrt() * turn.sin(), z);
    let plane = DQuat::from_rotation_arc(DVec3::Z, pole) * DQuat::from_rotation_z(rng.next() * std::f64::consts::TAU);
    Orbit { semi_major_au, eccentricity, period_years: (semi_major_au.powi(3) / mass).sqrt(), plane, phase: rng.next() * std::f64::consts::TAU }
}

/// The system of a made-up star. `seed` is the star's own number.
pub fn system_of(star: &Star, seed: u64, table: &[Multiplicity]) -> System {
    let mut rng = Rng::seeded(seed ^ 0x5953_5445_4d53);
    let first = Member { about: stars::describe(star, seed), abs_mag: star.abs_mag as f64, kind: star.kind };
    let primary = first.about.mass;
    let wanted = count(table, primary, &mut rng);

    let mut root = Part::Star(first.clone());
    let mut size = first.about.radius * SOLAR_RADIUS_AU;
    // Companions are added from the inside out, each around all that is there already.
    let mut distances: Vec<f64> = (1..wanted).map(|_| separation(primary, primary * 1.5, &mut rng)).collect();
    distances.sort_by(f64::total_cmp);
    for mut a in distances {
        // Lighter than the first star, evenly spread, and still a star.
        let low = if primary < 0.57 { 0.5 } else { 0.1 };
        let mass = (primary * (low + (1.0 - low) * rng.next())).max(LIGHTEST_STAR).min(primary);
        let abs_mag = stars::abs_mag_of_mass(mass);
        let star = Star { abs_mag: abs_mag as f32, kind: Kind::MainSequence, ..*star };
        let member = Member { about: stars::describe(&star, rng.next().to_bits()), abs_mag, kind: Kind::MainSequence };
        let total = root.mass() + mass;
        // Close pairs have been made round by their tides.
        let mut e = if (a.powi(3) / total).sqrt() * 365.25 < 12.0 { 0.0 } else { 0.85 * rng.next() };
        // The stars must not touch at their closest, and what is inside must be left alone.
        let inner = root.reach_au();
        let fits = |a: f64, e: f64| a * (1.0 - e) > 3.0 * (size + member.about.radius * SOLAR_RADIUS_AU) && (inner == 0.0 || a >= inner * stable_ratio(e, mass / root.mass()));
        if !fits(a, e) {
            e = 0.0;
            let least = (3.0 * (size + member.about.radius * SOLAR_RADIUS_AU)).max(inner * stable_ratio(0.0, mass / root.mass())) * 1.05;
            a = a.max(least);
        }
        if a > WIDEST_AU {
            break; // no room left for another star
        }
        let o = orbit(a, e, total, &mut rng);
        size = size.max(a * (1.0 + e));
        root = Part::Pair(Box::new(root), Box::new(Part::Star(member)), o);
    }
    System { root }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::procedural::StarId;

    fn star(abs_mag: f32) -> Star {
        Star { position: DVec3::ZERO, abs_mag, color: [1.0; 3], kind: Kind::MainSequence, id: StarId { level: 0, cell: [0; 3], index: 0 } }
    }

    /// Shares of systems with one, two, three and more stars among `n` stars of a brightness.
    fn shares(abs_mag: f32, n: u64) -> [f64; 4] {
        let mut by_count = [0u32; 4];
        for seed in 0..n {
            let system = system_of(&star(abs_mag), seed.wrapping_mul(0x9E37_79B9_7F4A_7C15) | 1, &MULTIPLICITY);
            by_count[system.stars().min(4) - 1] += 1;
        }
        by_count.map(|c| c as f64 / n as f64)
    }

    #[test]
    fn stars_have_companions_as_often_as_real_ones() {
        // Sun-like stars: 56% single, 33% double, 8% triple, 3% more.
        let sun = shares(4.83, 40_000);
        for (got, expected) in sun.iter().zip([0.56, 0.33, 0.08, 0.03]) {
            assert!((got - expected).abs() < 0.012, "sun-like: {sun:?}");
        }
        // M dwarfs: 76.7% single, 20% double, 3% triple, 0.3% more.
        let dwarf = shares(12.0, 40_000);
        for (got, expected) in dwarf.iter().zip([0.767, 0.20, 0.03, 0.003]) {
            assert!((got - expected).abs() < 0.012, "M dwarfs: {dwarf:?}");
        }
        // O stars are almost never alone.
        assert!(shares(-5.0, 5_000)[0] < 0.08);
    }

    #[test]
    fn systems_hold_together() {
        for seed in 1..20_000u64 {
            let system = system_of(&star(if seed % 3 == 0 { -2.0 } else { 4.0 }), seed.wrapping_mul(0xD6E8_FEB8_6659_FD93) | 1, &MULTIPLICITY);
            let orbits = system.orbits();
            for pair in orbits.windows(2) {
                let (inner, outer) = (pair[0], pair[1]);
                // The outer orbit clears the inner one by the margin of the stability rule (for
                // any mass ratio the rule asks for at least 2.8).
                assert!(outer.semi_major_au * (1.0 - outer.eccentricity) > 2.8 * inner.semi_major_au * (1.0 + inner.eccentricity) * 0.99, "{:?}", system.lines());
            }
            for o in &orbits {
                assert!(o.semi_major_au <= WIDEST_AU && o.eccentricity < 0.85 && o.period_years > 0.0);
                // The stars of a pair never touch.
                assert!(o.semi_major_au * (1.0 - o.eccentricity) > 2.0 * SOLAR_RADIUS_AU * 0.1, "{:?}", system.lines());
            }
            assert!(system.stars() <= MOST_STARS);
        }
    }

    #[test]
    fn the_same_star_has_the_same_system() {
        let (a, b) = (system_of(&star(3.0), 99, &MULTIPLICITY), system_of(&star(3.0), 99, &MULTIPLICITY));
        assert_eq!(a, b);
        assert_eq!(a.places_au(12.5), b.places_au(12.5));
    }

    #[test]
    fn orbits_are_ellipses_that_close() {
        let mut rng = Rng::seeded(5);
        let o = orbit(10.0, 0.6, 2.0, &mut rng);
        assert!((o.period_years - (1000.0f64 / 2.0).sqrt()).abs() < 1e-9);
        // Back at the same place after one turn, and between the nearest and the farthest point.
        assert!((o.offset_au(3.0) - o.offset_au(3.0 + o.period_years)).length() < 1e-6);
        for k in 0..200 {
            let r = o.offset_au(k as f64 * o.period_years / 200.0).length();
            assert!(r > 10.0 * 0.4 - 1e-9 && r < 10.0 * 1.6 + 1e-9, "{r}");
        }
        // The companion adds its light: two equal stars are 0.75 magnitudes brighter than one.
        let double = System { root: Part::Pair(Box::new(Part::Star(member(5.0))), Box::new(Part::Star(member(5.0))), o) };
        assert!((double.abs_mag() - (5.0 - 0.7526)).abs() < 1e-3);
    }

    fn member(abs_mag: f64) -> Member {
        Member { about: stars::describe(&star(abs_mag as f32), 1), abs_mag, kind: Kind::MainSequence }
    }
}
