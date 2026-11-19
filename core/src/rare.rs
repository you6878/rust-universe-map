//! Objects that are not ordinary stars: the dead ones (neutron stars, black holes), the failed
//! ones (brown dwarfs), planets without a star, and the very rare. Each kind is a layer of its
//! own over the galaxy, spread like the stars are, with as many in all as the map says.
//!
//! A kind lives in cubes of a size that suits how rare it is: about four objects to a cube near
//! the Sun. That is what makes "the nearest black hole" a question with a quick answer: a few
//! cubes around the camera are computed, not billions of stars looked through. As with the
//! stars, a cube always gives the same objects, and nothing is stored.

use glam::DVec3;

use crate::procedural::{Model, Rng};

/// A kind of rare object: its name and how many the galaxy has.
#[derive(Clone, Debug, PartialEq)]
pub struct Class {
    pub name: String,
    pub total: f64,
}

/// The kinds of the default map, with the numbers of map/SOURCES.md. Brown dwarfs and
/// free-floating planets are given per star there: one for six stars, and 21 a star.
pub fn default_classes(stars: f64) -> Vec<Class> {
    let class = |name: &str, total: f64| Class { name: name.to_owned(), total };
    vec![
        class("neutron star", 1e9),
        class("black hole", 1e8),
        class("radio pulsar", 3e4),
        class("magnetar", 30.0),
        class("x-ray binary", 350.0),
        class("hypervelocity star", 1000.0),
        class("brown dwarf", stars / 6.0),
        class("free-floating planet", stars * 21.0),
    ]
}

/// The address of a rare object: its kind, its cube and its place in the list of the cube.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct RareId {
    pub class: u16,
    pub cell: [i64; 3],
    pub index: u32,
}

#[derive(Clone, Debug, PartialEq)]
pub struct Rare {
    pub id: RareId,
    /// Galactic, parsecs.
    pub position: DVec3,
    /// What it is like, in a few words and numbers.
    pub about: String,
}

/// The side of the cubes of a kind, in parsecs: a power of two times two parsecs, such that a
/// cube near the Sun holds about four.
pub fn side(model: &Model, systems: f64, class: &Class) -> f64 {
    let density = class.total * model.density(DVec3::ZERO) / systems;
    let wanted = (4.0 / density.max(1e-300)).cbrt();
    2.0 * 2f64.powi((wanted / 2.0).log2().round().clamp(0.0, 14.0) as i32)
}

fn poisson(mean: f64, rng: &mut Rng) -> usize {
    if mean > 30.0 {
        // Many: the bell curve is close enough.
        let (u, v) = (rng.next().max(1e-300), rng.next());
        return (mean + mean.sqrt() * (-2.0 * u.ln()).sqrt() * (std::f64::consts::TAU * v).cos()).round().max(0.0) as usize;
    }
    let (limit, mut product, mut n) = ((-mean).exp(), rng.next(), 0);
    while product > limit {
        product *= rng.next();
        n += 1;
    }
    n
}

/// A few words on one object of a kind, from draws of its own.
fn describe(name: &str, rng: &mut Rng) -> String {
    let mut between = |low: f64, high: f64| low * (high / low).powf(rng.next());
    match name {
        "neutron star" => format!("neutron star, {:.2} solar masses, 12 km across", 1.1 + 0.9 * between(0.1, 1.0)),
        "black hole" => format!("black hole, {:.1} solar masses", between(5.0, 30.0)),
        "radio pulsar" => format!("radio pulsar, turns once in {:.3} s", between(0.03, 3.0)),
        "magnetar" => format!("magnetar, turns once in {:.1} s, magnetic field {:.0e} tesla", between(2.0, 12.0), between(1e10, 1e11)),
        "x-ray binary" => format!("x-ray binary: a {} with a star of {:.1} solar masses feeding it", if between(0.1, 1.0) < 0.5 { "neutron star" } else { "black hole" }, between(0.3, 20.0)),
        "hypervelocity star" => format!("hypervelocity star, leaving the galaxy at {:.0} km/s", between(500.0, 1200.0)),
        "brown dwarf" => format!("brown dwarf, {:.0} Jupiter masses, {:.0} K", between(13.0, 80.0), between(300.0, 2200.0)),
        "free-floating planet" => format!("free-floating planet, {:.1} Earth masses", between(0.33, 6660.0f64.powf(0.35))),
        other => other.to_owned(),
    }
}

/// The objects of one kind in one cube. `systems` is the number of systems of the whole model.
pub fn generate(model: &Model, systems: f64, classes: &[Class], class: usize, cell: [i64; 3]) -> Vec<Rare> {
    let kind = &classes[class];
    let s = side(model, systems, kind);
    let corner = DVec3::new(cell[0] as f64, cell[1] as f64, cell[2] as f64) * s;
    let mut rng = Rng::seeded(model.seed ^ (class as u64 + 1).wrapping_mul(0x5241_5245) ^ (cell[0] as u64).wrapping_mul(0x9E37_79B9) ^ (cell[1] as u64).wrapping_mul(0x85EB_CA6B_0001) ^ (cell[2] as u64).wrapping_mul(0xC2B2_AE35_0000_0001));
    // How much of the galaxy is in the cube: the density at points spread over it.
    let n = if s > 400.0 { 8 } else { 3 };
    let (mut sum, mut max) = (0.0, 0.0f64);
    for i in 0..n {
        for j in 0..n {
            for k in 0..n {
                let d = model.density(corner + DVec3::new(i as f64 + rng.next(), j as f64 + rng.next(), k as f64 + rng.next()) / n as f64 * s);
                sum += d;
                max = max.max(d);
            }
        }
    }
    let mean = kind.total * (sum / (n * n * n) as f64 * s * s * s) / systems;
    let count = poisson(mean, &mut rng);
    let mut out = Vec::with_capacity(count);
    let mut tries = 0;
    while out.len() < count && tries < count * 60 + 200 {
        tries += 1;
        let position = corner + DVec3::new(rng.next(), rng.next(), rng.next()) * s;
        if rng.next() * max * 1.5 > model.density(position) {
            continue;
        }
        let id = RareId { class: class as u16, cell, index: out.len() as u32 };
        out.push(Rare { id, position, about: describe(&kind.name, &mut rng) });
    }
    out
}

/// The nearest objects of a kind to a place, nearest first, with their distances in parsecs.
/// Cubes are looked through in shells around the place, up to `shells` cubes away; if the kind is
/// so rare that none is found there, the list is empty.
pub fn nearest(model: &Model, systems: f64, classes: &[Class], class: usize, place: DVec3, count: usize, shells: i64) -> Vec<(Rare, f64)> {
    let s = side(model, systems, &classes[class]);
    let home = (place / s).floor();
    let home = [home.x as i64, home.y as i64, home.z as i64];
    let mut found: Vec<(Rare, f64)> = Vec::new();
    for shell in 0..=shells {
        for dx in -shell..=shell {
            for dy in -shell..=shell {
                for dz in -shell..=shell {
                    if dx.abs().max(dy.abs()).max(dz.abs()) != shell {
                        continue; // inner cubes were looked through before
                    }
                    for rare in generate(model, systems, classes, class, [home[0] + dx, home[1] + dy, home[2] + dz]) {
                        let d = rare.position.distance(place);
                        found.push((rare, d));
                    }
                }
            }
        }
        found.sort_by(|a, b| a.1.total_cmp(&b.1));
        // Anything in a cube not looked through yet is at least this far.
        if found.len() >= count && found[count - 1].1 <= shell as f64 * s {
            break;
        }
    }
    found.truncate(count);
    found
}

#[cfg(test)]
mod tests {
    use super::*;

    fn setup() -> (Model, f64, Vec<Class>) {
        let model = Model::default();
        let systems = model.total();
        let classes = default_classes(systems * 1.35);
        (model, systems, classes)
    }

    #[test]
    fn the_galaxy_has_as_many_as_the_map_says() {
        let (model, systems, classes) = setup();
        // Hypervelocity stars: a thousand in all. Count every cube that can hold any.
        let class = classes.iter().position(|c| c.name == "hypervelocity star").unwrap();
        let s = side(&model, systems, &classes[class]);
        let reach = (22_000.0 / s).ceil() as i64;
        let mut total = 0;
        for x in -reach..reach {
            for y in -(8_000.0 / s).ceil() as i64..(8_000.0 / s).ceil() as i64 {
                for z in -reach..reach {
                    total += generate(&model, systems, &classes, class, [x, y, z + (8_000.0 / s) as i64]).len();
                }
            }
        }
        assert!((850..1150).contains(&total), "{total} hypervelocity stars");
    }

    #[test]
    fn cubes_suit_how_rare_a_kind_is() {
        let (model, systems, classes) = setup();
        let side_of = |name: &str| side(&model, systems, classes.iter().find(|c| c.name == name).unwrap());
        // Planets without a star are everywhere; black holes take tens of parsecs to find; a
        // magnetar, thousands.
        assert_eq!(side_of("free-floating planet"), 2.0);
        assert!((16.0..=128.0).contains(&side_of("black hole")), "{}", side_of("black hole"));
        assert!(side_of("magnetar") >= 4096.0);
    }

    #[test]
    fn the_nearest_is_found_and_is_the_nearest() {
        let (model, systems, classes) = setup();
        for name in ["black hole", "brown dwarf", "neutron star"] {
            let class = classes.iter().position(|c| c.name == name).unwrap();
            let place = DVec3::new(100.0, 20.0, 3000.0);
            let found = nearest(&model, systems, &classes, class, place, 3, 6);
            assert_eq!(found.len(), 3, "{name}");
            assert!(found[0].1 <= found[1].1 && found[1].1 <= found[2].1);
            // Looking through more cubes finds nothing nearer.
            let wider = nearest(&model, systems, &classes, class, place, 1, 8);
            assert_eq!(wider[0].0, found[0].0, "{name}");
            // The same again, the same answer.
            assert_eq!(nearest(&model, systems, &classes, class, place, 3, 6), found);
            // And the cube gives the object back from its address.
            let again = generate(&model, systems, &classes, class, found[0].0.id.cell);
            assert_eq!(again[found[0].0.id.index as usize], found[0].0);
        }
    }
}
