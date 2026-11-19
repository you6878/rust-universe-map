//! Finding where to go by name: the bodies of the Solar System, the stars of the catalogs that
//! have a name, and the far side of the Milky Way.

use glam::DVec3;

use crate::Dest;

/// A star of the catalogs with a name, kept so that it can be travelled to.
pub struct NamedStar {
    /// All its names, the first the one to show.
    pub names: Vec<String>,
    /// The names as they are compared: lower case, without spaces.
    keys: Vec<String>,
    /// Equatorial position in parsecs, from the Sun.
    pub position_pc: DVec3,
    /// Magnitude seen from the Earth: the brighter (smaller) comes first among equal matches.
    pub app_mag: f32,
}

impl NamedStar {
    pub fn new(names: Vec<String>, position_pc: DVec3, app_mag: f32) -> NamedStar {
        let keys = names.iter().map(|n| key(n)).collect();
        NamedStar { names, keys, position_pc, app_mag }
    }
}

/// Whether the names of a star are real names, and not only its number in the file.
pub fn has_name(id: i64, names: &[String]) -> bool {
    names.iter().any(|n| !n.is_empty() && *n != id.to_string())
}

/// How a name is compared: case and spaces do not count, so "hip32349" finds "HIP 32349".
fn key(name: &str) -> String {
    name.chars().filter(|c| !c.is_whitespace()).flat_map(char::to_lowercase).collect()
}

/// How well a name matches what was typed: 0 the same, 1 it starts so, 2 it has it inside.
fn closeness(name_key: &str, query_key: &str) -> Option<u8> {
    if name_key == query_key {
        Some(0)
    } else if name_key.starts_with(query_key) {
        Some(1)
    } else if name_key.contains(query_key) {
        Some(2)
    } else {
        None
    }
}

/// Something found, and the name it was found by.
#[derive(Clone, Debug)]
pub struct Hit {
    pub dest: Dest,
    pub name: String,
}

const FAR_SIDE: [&str; 3] = ["far side of the Milky Way", "farside", "galaxy"];

/// What answers to `query`, the best first, at most `limit`: exact names before names that
/// start so before names that have it inside; then bodies before stars, and bright stars before
/// faint ones.
pub fn search(query: &str, bodies: &[String], stars: &[NamedStar], limit: usize) -> Vec<Hit> {
    let q = key(query);
    if q.is_empty() {
        return Vec::new();
    }
    // (closeness, group, brightness, order) and what was found.
    let mut found: Vec<((u8, u8, f32, usize), Hit)> = Vec::new();
    for (i, name) in bodies.iter().enumerate() {
        if let Some(c) = closeness(&key(name), &q) {
            found.push(((c, 0, 0.0, i), Hit { dest: Dest::Body(i), name: name.clone() }));
        }
    }
    if let Some(c) = FAR_SIDE.iter().filter_map(|n| closeness(&key(n), &q)).min() {
        found.push(((c, 1, 0.0, 0), Hit { dest: Dest::FarSide, name: FAR_SIDE[0].into() }));
    }
    for (i, star) in stars.iter().enumerate() {
        // The best of its names, and that name to show when it is not the first.
        if let Some((c, n)) = star.keys.iter().enumerate().filter_map(|(n, k)| closeness(k, &q).map(|c| (c, n))).min() {
            let name = if n == 0 { star.names[0].clone() } else { format!("{} ({})", star.names[0], star.names[n]) };
            let mag = if star.app_mag.is_finite() { star.app_mag } else { 99.0 };
            found.push(((c, 2, mag, i), Hit { dest: Dest::Star(i), name }));
        }
    }
    found.sort_by(|a, b| a.0 .0.cmp(&b.0 .0).then(a.0 .1.cmp(&b.0 .1)).then(a.0 .2.total_cmp(&b.0 .2)).then(a.0 .3.cmp(&b.0 .3)));
    found.into_iter().take(limit).map(|(_, hit)| hit).collect()
}

/// The one thing a name means, for `--travel`: the best match, if it is the whole name.
pub fn exact(query: &str, bodies: &[String], stars: &[NamedStar]) -> Option<Dest> {
    let q = key(query);
    let best = search(query, bodies, stars, 1).into_iter().next()?;
    let names: Vec<String> = match best.dest {
        Dest::Body(i) => vec![bodies[i].clone()],
        Dest::Star(i) => stars[i].names.clone(),
        Dest::FarSide => FAR_SIDE.iter().map(|s| s.to_string()).collect(),
    };
    names.iter().any(|n| key(n) == q).then_some(best.dest)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn bodies() -> Vec<String> {
        ["Sun", "Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune", "Moon"].map(String::from).to_vec()
    }

    fn stars() -> Vec<NamedStar> {
        let star = |names: &[&str], mag: f32| NamedStar::new(names.iter().map(|s| s.to_string()).collect(), DVec3::X, mag);
        vec![star(&["Sirius", "HIP 32349"], -1.46), star(&["Sirrah", "Alpheratz", "HIP 677"], 2.06), star(&["Vega", "HIP 91262"], 0.03), star(&["Marsic", "HIP 80463"], 5.4), star(&["Rigel", "HIP 24436"], 0.13)]
    }

    fn names(hits: &[Hit]) -> Vec<&str> {
        hits.iter().map(|h| h.name.as_str()).collect()
    }

    #[test]
    fn a_whole_name_comes_before_names_that_start_with_it() {
        let hits = search("mars", &bodies(), &stars(), 5);
        assert_eq!(names(&hits), ["Mars", "Marsic"]);
        assert!(matches!(hits[0].dest, Dest::Body(4)));
        assert!(matches!(hits[1].dest, Dest::Star(3)));
    }

    #[test]
    fn case_and_spaces_do_not_count_and_other_names_are_found_too() {
        let hits = search("hip32349", &bodies(), &stars(), 5);
        assert_eq!(names(&hits), ["Sirius (HIP 32349)"]);
        assert!(matches!(hits[0].dest, Dest::Star(0)));
        assert!(matches!(search("ALPHERATZ", &bodies(), &stars(), 5)[0].dest, Dest::Star(1)));
    }

    #[test]
    fn among_equal_matches_the_brighter_star_comes_first() {
        // Both start with "sir": Sirius (-1.46) is brighter than Sirrah (2.06).
        assert_eq!(names(&search("Sir", &bodies(), &stars(), 5)), ["Sirius", "Sirrah"]);
        // "ga" is inside Vega and the galaxy: the galaxy is no star, so it comes first.
        assert_eq!(names(&search("ga", &bodies(), &stars(), 5)), ["far side of the Milky Way", "Vega"]);
    }

    #[test]
    fn the_list_is_cut_and_nothing_typed_finds_nothing() {
        assert_eq!(search("r", &bodies(), &stars(), 3).len(), 3);
        assert!(search("  ", &bodies(), &stars(), 5).is_empty());
        assert!(search("zzz", &bodies(), &stars(), 5).is_empty());
    }

    #[test]
    fn travel_on_the_command_line_needs_the_whole_name() {
        assert!(matches!(exact("sirius", &bodies(), &stars()), Some(Dest::Star(0))));
        assert!(matches!(exact("HIP 91262", &bodies(), &stars()), Some(Dest::Star(2))));
        assert!(matches!(exact("farside", &bodies(), &stars()), Some(Dest::FarSide)));
        assert!(matches!(exact("earth", &bodies(), &stars()), Some(Dest::Body(3))));
        assert!(exact("Siri", &bodies(), &stars()).is_none());
    }

    #[test]
    fn a_number_alone_is_no_name() {
        assert!(!has_name(42, &["42".to_string()]));
        assert!(has_name(42, &["42".to_string(), "Vega".to_string()]));
    }
}
