//! The planets of a made-up star system. Like everything about a made-up star, they follow from
//! the star's own number: the same star has the same planets on every visit, and nothing is
//! stored. How many there are, how large and how far apart follows what is known of real planets
//! (map/SOURCES.md); where that leaves a choice, the choice made here is said in the comments.

use glam::DVec3;

use crate::procedural::Rng;
use crate::systems::{orbit, Orbit, Part, System};

pub const EARTH_RADIUS_AU: f64 = 4.2635e-5;
const EARTH_MASSES_PER_SUN: f64 = 332_946.0;
const JUPITER_MASS: f64 = 317.8;
/// The habitable zone of the Sun, in AU: between these distances a rocky planet can keep liquid
/// water on its surface. Nearer, the water boils away (as on Venus); farther, it freezes (as on
/// Mars). These are the limits of Kopparapu et al. 2013 (water loss, maximum greenhouse); for
/// another star they scale with the square root of its power.
pub const HABITABLE_ZONE: (f64, f64) = (0.99, 1.70);
/// Planets a star has at most.
const MOST_PLANETS: usize = 12;
/// Neighbours are at least this many mutual Hill radii apart (real ones: nearly all over 10,
/// most often 20), and their periods differ by at least this factor.
const LEAST_HILL_RADII: f64 = 8.0;
const LEAST_PERIOD_RATIO: f64 = 1.2;

/// How many planets stars have, a choice the map can change: the mean number per star (the
/// surveys give 1.6 per star between 0.5 and 10 AU, and 2.5 small ones within 200 days for M
/// dwarfs, each for a part of the range only), and the share of stars of one solar mass and the
/// Sun's metals with a giant planet (it grows with the mass and the metals of the star).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Abundance {
    pub mean_planets: f64,
    pub giants_at_one_solar_mass: f64,
    /// Of the stars with a giant, the share where it is a hot Jupiter.
    pub hot_share: f64,
}

pub const ABUNDANCE: Abundance = Abundance { mean_planets: 4.5, giants_at_one_solar_mass: 0.07, hot_share: 0.15 };

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Kind {
    /// Rock so hot that its surface is molten.
    Lava,
    Rocky,
    /// Rock with liquid water: in the zone where a planet like the Earth keeps its oceans.
    Temperate,
    /// Rock under ice.
    Icy,
    /// Between the Earth and Neptune in size, with a thick air.
    SubNeptune,
    IceGiant,
    GasGiant,
}

impl Kind {
    pub fn name(self) -> &'static str {
        match self {
            Kind::Lava => "lava planet",
            Kind::Rocky => "rocky planet",
            Kind::Temperate => "temperate rocky planet",
            Kind::Icy => "icy planet",
            Kind::SubNeptune => "sub-Neptune",
            Kind::IceGiant => "ice giant",
            Kind::GasGiant => "gas giant",
        }
    }

    /// A number for the kind, for the shader that draws its surface.
    pub fn code(self) -> u8 {
        match self {
            Kind::Lava => 1,
            Kind::Rocky => 2,
            Kind::Temperate => 3,
            Kind::Icy => 4,
            Kind::SubNeptune => 5,
            Kind::IceGiant => 6,
            Kind::GasGiant => 7,
        }
    }

    /// The color it is drawn in.
    pub fn color(self) -> [f32; 3] {
        match self {
            Kind::Lava => [0.95, 0.35, 0.15],
            Kind::Rocky => [0.62, 0.55, 0.48],
            Kind::Temperate => [0.35, 0.55, 0.85],
            Kind::Icy => [0.85, 0.92, 1.0],
            Kind::SubNeptune => [0.55, 0.75, 0.8],
            Kind::IceGiant => [0.45, 0.65, 0.95],
            Kind::GasGiant => [0.85, 0.72, 0.55],
        }
    }
}

/// What a planet goes around.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Host {
    /// One star of the system, by its place in `System::places_au`.
    Star(usize),
    /// All the stars of the system together.
    All,
}

#[derive(Clone, Debug, PartialEq)]
pub struct Planet {
    pub kind: Kind,
    pub host: Host,
    /// In Earth masses and Earth radii.
    pub mass: f64,
    pub radius: f64,
    pub orbit: Orbit,
    /// Temperature of a dark body at its distance, in kelvins.
    pub temperature: f64,
    pub moons: u32,
    pub rings: bool,
}

impl Planet {
    pub fn line(&self) -> String {
        let size = if self.mass > 0.3 * JUPITER_MASS { format!("{:.2} Jupiter masses", self.mass / JUPITER_MASS) } else { format!("{:.2} Earth masses", self.mass) };
        let period = if self.orbit.period_years < 1.0 { format!("{:.1} days", self.orbit.period_years * 365.25) } else { format!("{:.1} years", self.orbit.period_years) };
        let around = match self.host {
            Host::Star(i) => format!("star {}", i + 1),
            Host::All => "all the stars".to_owned(),
        };
        let extras = match (self.moons, self.rings) {
            (0, false) => String::new(),
            (0, true) => ", rings".to_owned(),
            (n, false) => format!(", {n} moon{}", if n == 1 { "" } else { "s" }),
            (n, true) => format!(", {n} moon{}, rings", if n == 1 { "" } else { "s" }),
        };
        format!("{}: {size}, {:.2} Earth radii, {:.3} AU from {around}, once in {period}, {:.0} K{extras}", self.kind.name(), self.radius, self.orbit.semi_major_au, self.temperature)
    }

    /// Its absolute magnitude in the light it reflects, given that of what it goes around.
    pub fn abs_mag(&self, host_abs_mag: f64) -> f64 {
        // A third of the light that falls on its disc, half of that towards a viewer to the side.
        let share = 0.3 * 0.5 * (self.radius * EARTH_RADIUS_AU / (2.0 * self.orbit.semi_major_au)).powi(2);
        host_abs_mag - 2.5 * share.log10()
    }
}

fn poisson(mean: f64, rng: &mut Rng) -> usize {
    let (limit, mut product, mut n) = ((-mean).exp(), rng.next(), 0);
    while product > limit {
        product *= rng.next();
        n += 1;
    }
    n
}

fn normal(rng: &mut Rng) -> f64 {
    let (u, v) = (rng.next().max(1e-300), rng.next());
    (-2.0 * u.ln()).sqrt() * (std::f64::consts::TAU * v).cos()
}

/// Mass in Earth masses of a planet of a radius in Earth radii that is not a giant: rock up to
/// 1.5 Earth radii, then more and more air.
fn mass_of_radius(radius: f64) -> f64 {
    if radius < 1.5 { radius.powf(3.7) } else { 2.7 * radius.powf(1.3) }
}

/// The radius of the Hill sphere of two neighbours together, in AU.
fn mutual_hill(a1: f64, m1: f64, a2: f64, m2: f64, star: f64) -> f64 {
    (a1 + a2) / 2.0 * ((m1 + m2) / (3.0 * star * EARTH_MASSES_PER_SUN)).cbrt()
}

/// Where planets can go around a host: the widest and the narrowest orbit, in AU. A companion
/// star clears a zone around its own orbit (Holman and Wiegert 1999).
fn room(part: &Part, index: &mut usize, outer: Option<(f64, f64, f64)>, out: &mut Vec<(Host, f64, f64, f64, f64, f64)>) {
    match part {
        Part::Star(member) => {
            // Not inside the star's own glare, and not farther than its companion allows.
            let near = (member.about.radius * 0.00465 * 3.0).max(0.01);
            let far = outer.map_or(60.0 * member.about.mass.sqrt(), |(a, e, mu)| a * (0.464 - 0.380 * mu - 0.631 * e + 0.586 * mu * e + 0.150 * e * e - 0.198 * mu * e * e));
            out.push((Host::Star(*index), near, far, member.about.mass, member.about.power(), member.abs_mag));
            *index += 1;
        }
        Part::Pair(a, b, o) => {
            let (ma, mb) = (a.mass(), b.mass());
            // Each side keeps what the other leaves it; the nearer limit of the two orbits counts.
            let tighter = |own: (f64, f64, f64)| outer.map_or(own, |up| if up.0 < own.0 { up } else { own });
            room(a, index, Some(tighter((o.semi_major_au, o.eccentricity, mb / (ma + mb)))), out);
            room(b, index, Some(tighter((o.semi_major_au, o.eccentricity, ma / (ma + mb)))), out);
        }
    }
}

/// The planets of a system. `seed` is the number of the star the system was made from.
pub fn planets_of(system: &System, seed: u64, abundance: &Abundance) -> Vec<Planet> {
    let mut rng = Rng::seeded(seed ^ 0x504c_414e_4554);
    let mut hosts = Vec::new();
    room(&system.root, &mut 0, None, &mut hosts);
    // Around a close pair as a whole, outside the zone the pair stirs up.
    if let Part::Pair(a, b, o) = &system.root {
        if system.stars() == 2 && o.semi_major_au < 3.0 {
            let (e, mu) = (o.eccentricity, b.mass() / (a.mass() + b.mass()));
            let near = o.semi_major_au * (1.60 + 5.10 * e - 2.22 * e * e + 4.12 * mu - 4.27 * e * mu - 5.09 * mu * mu + 4.61 * e * e * mu * mu);
            let light: f64 = hosts.iter().map(|h| h.4).sum();
            let abs_mag = -2.5 * hosts.iter().map(|h| 10f64.powf(-0.4 * h.5)).sum::<f64>().log10();
            hosts.push((Host::All, near, 60.0 * system.root.mass().sqrt(), system.root.mass(), light, abs_mag));
        }
    }

    let mut planets = Vec::new();
    for (host, near, far, star_mass, light, _) in hosts {
        if far <= near * 1.5 {
            continue;
        }
        // The metals of the star, compared with the Sun's, on a logarithmic scale.
        let metals = -0.05 + 0.2 * normal(&mut rng);
        let giant = rng.next() < (abundance.giants_at_one_solar_mass * star_mass * 10f64.powf(1.2 * metals)).min(0.6);
        let hot = giant && rng.next() < abundance.hot_share;
        let count = (poisson(abundance.mean_planets, &mut rng) + giant as usize).min(MOST_PLANETS);
        if count == 0 {
            continue;
        }
        // Which of the planets is the giant: the innermost if it is a hot one, else one of the outer ones.
        let giant_at = if !giant { usize::MAX } else if hot { 0 } else { count / 2 + ((count - count / 2) as f64 * rng.next()) as usize };

        // Orbits from the inside out. Planets of small stars huddle closer.
        let mut a = if hot { 0.02 + 0.06 * rng.next() } else { 0.04 * star_mass.sqrt() * 10f64.powf(0.9 * rng.next()) }.max(near);
        let mut last: Option<(f64, f64)> = None;
        for i in 0..count {
            let is_giant = i == giant_at.min(count - 1) && giant;
            let (radius, mass) = if is_giant {
                // Between a third of Jupiter and ten Jupiters, evenly on a logarithmic scale.
                let mass = JUPITER_MASS * 10f64.powf(-0.5 + 1.5 * rng.next());
                (11.2 * (mass / JUPITER_MASS).powf(0.05).min(1.2), mass)
            } else {
                // Small ones are the many: evenly on a logarithmic scale from half an Earth to
                // three Earths across, and now and then one the size of Neptune.
                let radius = if rng.next() < 0.08 { 3.0 + 3.0 * rng.next() } else { 0.5 * 6f64.powf(rng.next()) };
                (radius, mass_of_radius(radius))
            };
            // The next orbit: a period two times longer, more or less, and well clear of the last.
            if let Some((a_last, m_last)) = last {
                a = a_last * (LEAST_PERIOD_RATIO.max(2.0 * 10f64.powf(0.18 * normal(&mut rng)))).powf(2.0 / 3.0);
                for _ in 0..40 {
                    if (a - a_last) >= LEAST_HILL_RADII * mutual_hill(a_last, m_last, a, mass, star_mass) {
                        break;
                    }
                    a *= 1.1;
                }
            }
            if a > far {
                break;
            }
            // Orbits in systems of several planets are nearly round.
            let e = (0.05 * (-2.0 * rng.next().max(1e-12).ln()).sqrt()).min(0.3);
            let temperature = 255.0 * light.powf(0.25) / a.sqrt();
            let kind = if radius >= 7.0 {
                Kind::GasGiant
            } else if radius >= 4.0 {
                Kind::IceGiant
            } else if radius >= 1.6 {
                Kind::SubNeptune
            } else if temperature > 1000.0 {
                Kind::Lava
            } else if (HABITABLE_ZONE.0 * light.sqrt()..HABITABLE_ZONE.1 * light.sqrt()).contains(&a) {
                Kind::Temperate
            } else if temperature < 150.0 {
                Kind::Icy
            } else {
                Kind::Rocky
            };
            // Moons and rings of planets around other stars have not been seen: these numbers
            // are made after the solar system, where giants have many moons and often rings.
            let (moons, rings) = match kind {
                Kind::GasGiant | Kind::IceGiant => (2 + (rng.next() * 10.0) as u32, rng.next() < 0.5),
                Kind::SubNeptune => ((rng.next() * 3.0) as u32, rng.next() < 0.05),
                _ if temperature > 600.0 => (0, false),
                _ => ((rng.next() * rng.next() * 3.0) as u32, false),
            };
            planets.push(Planet { kind, host, mass, radius, orbit: orbit(a, e, star_mass, &mut rng), temperature, moons, rings });
            last = Some((a, mass));
        }
    }
    planets
}

/// Where the planets are, in AU from the centre of mass of the system, at a time in years since
/// J2000, with the absolute magnitude each has in reflected light.
pub fn places_au(system: &System, planets: &[Planet], years: f64) -> Vec<(DVec3, f64)> {
    let stars = system.places_au(years);
    let all_light = -2.5 * stars.iter().map(|(_, m)| 10f64.powf(-0.4 * m.abs_mag)).sum::<f64>().log10();
    planets
        .iter()
        .map(|p| {
            let (centre, light) = match p.host {
                Host::Star(i) => (stars[i].0, stars[i].1.abs_mag),
                Host::All => (DVec3::ZERO, all_light),
            };
            (centre + p.orbit.offset_au(years), p.abs_mag(light))
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::procedural::{Kind as StarKind, Star, StarId};
    use crate::systems::{system_of, MULTIPLICITY};

    fn systems(abs_mag: f32, n: u64) -> impl Iterator<Item = (System, Vec<Planet>)> {
        (1..=n).map(move |i| {
            let seed = i.wrapping_mul(0x9E37_79B9_7F4A_7C15) | 1;
            let star = Star { position: DVec3::ZERO, abs_mag, color: [1.0; 3], kind: StarKind::MainSequence, id: StarId { level: 0, cell: [0; 3], index: 0 } };
            let system = system_of(&star, seed, &MULTIPLICITY);
            let planets = planets_of(&system, seed, &ABUNDANCE);
            (system, planets)
        })
    }

    #[test]
    fn planets_are_as_common_as_the_surveys_say() {
        let n = 20_000;
        // Sun-like stars, single ones: giants at about one star in fourteen for the Sun's mass
        // (the surveys: 10.5% for sun-like stars up to twice as heavy), hot Jupiters at about 1%.
        let (mut stars, mut planets, mut giants, mut hot, mut temperate) = (0, 0, 0, 0, 0);
        for (system, list) in systems(4.83, n) {
            if system.stars() > 1 {
                continue;
            }
            stars += 1;
            planets += list.len();
            giants += list.iter().any(|p| p.kind == Kind::GasGiant) as u32;
            hot += list.iter().any(|p| p.kind == Kind::GasGiant && p.orbit.semi_major_au < 0.1) as u32;
            temperate += list.iter().any(|p| p.kind == Kind::Temperate) as u32;
        }
        let share = |c: u32| c as f64 / stars as f64;
        assert!((3.5..5.0).contains(&(planets as f64 / stars as f64)), "{} planets a star", planets as f64 / stars as f64);
        assert!((0.05..0.11).contains(&share(giants)), "giants: {}", share(giants));
        assert!((0.005..0.02).contains(&share(hot)), "hot Jupiters: {}", share(hot));
        // Earth-size planets where water is liquid: the surveys say 0.37 to 0.60 a star, with
        // error bars as large as the numbers; here somewhat fewer.
        assert!((0.15..0.6).contains(&share(temperate)), "temperate: {}", share(temperate));

        // M dwarfs: few giants, and their small planets close in.
        let (mut stars, mut giants, mut close) = (0, 0, 0);
        for (system, list) in systems(11.0, n) {
            if system.stars() > 1 {
                continue;
            }
            stars += 1;
            giants += list.iter().any(|p| p.kind == Kind::GasGiant) as u32;
            close += list.iter().filter(|p| p.orbit.period_years * 365.25 < 200.0 && (1.0..4.0).contains(&p.radius)).count();
        }
        assert!(giants as f64 / (stars as f64) < 0.04, "giants of M dwarfs: {}", giants as f64 / stars as f64);
        // The surveys: 2.5 planets of one to four Earth radii within 200 days.
        let close = close as f64 / stars as f64;
        assert!((1.2..3.2).contains(&close), "close small planets of M dwarfs: {close}");
    }

    #[test]
    fn planets_keep_clear_of_each_other_and_of_companion_stars() {
        for (system, list) in systems(4.0, 10_000).chain(systems(-1.0, 3_000)) {
            assert!(list.len() <= MOST_PLANETS * system.stars() + MOST_PLANETS);
            let stars = system.places_au(0.0);
            let orbits = system.orbits();
            for host in (0..stars.len()).map(Host::Star).chain([Host::All]) {
                let of_host: Vec<&Planet> = list.iter().filter(|p| p.host == host).collect();
                let mass = match host {
                    Host::Star(i) => stars[i].1.about.mass,
                    Host::All => system.root.mass(),
                };
                for pair in of_host.windows(2) {
                    let (inner, outer) = (pair[0], pair[1]);
                    assert!(outer.orbit.period_years / inner.orbit.period_years >= LEAST_PERIOD_RATIO * 0.999, "{} / {}", inner.line(), outer.line());
                    let hill = mutual_hill(inner.orbit.semi_major_au, inner.mass, outer.orbit.semi_major_au, outer.mass, mass);
                    assert!(outer.orbit.semi_major_au - inner.orbit.semi_major_au >= LEAST_HILL_RADII * hill * 0.999);
                }
                // A planet of one star of a pair stays well inside the pair's orbit; one around
                // a pair as a whole, well outside.
                if let (Some(tightest), Some(p)) = (orbits.first(), of_host.last()) {
                    match host {
                        Host::Star(_) => assert!(p.orbit.semi_major_au < 0.47 * orbits.iter().map(|o| o.semi_major_au).fold(f64::INFINITY, f64::min).max(tightest.semi_major_au) || stars.len() > 2, "{}", p.line()),
                        Host::All => assert!(of_host[0].orbit.semi_major_au > 1.6 * tightest.semi_major_au),
                    }
                }
            }
        }
    }

    #[test]
    fn the_same_star_has_the_same_planets() {
        let (a, b): (Vec<_>, Vec<_>) = (systems(5.0, 50).collect(), systems(5.0, 50).collect());
        assert_eq!(a, b);
        // And they move: a planet is elsewhere a quarter of its year later.
        let (system, planets) = a.into_iter().find(|(_, p)| !p.is_empty()).expect("a star with planets");
        let quarter = planets[0].orbit.period_years / 4.0;
        let (now, later) = (places_au(&system, &planets, 0.0), places_au(&system, &planets, quarter));
        assert!((now[0].0 - later[0].0).length() > 0.5 * planets[0].orbit.semi_major_au);
        // A planet shines by reflected light only: far fainter than its star.
        assert!(now[0].1 > system.abs_mag() + 15.0);
    }
}
