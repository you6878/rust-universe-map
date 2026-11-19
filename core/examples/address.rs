//! What is at an address, computed from the address alone:
//!
//!   cargo run --release --example address -- "Andromeda Galaxy/1"
//!   cargo run --release --example address -- "Andromeda Galaxy/123456789012/2"
//!   cargo run --release --example address -- "Milky Way/280"
//!
//! Reads map/universe-computed.json (the layers that need no other data) and prints the star,
//! its companions and its planets. Nothing is looked up: the same address gives the same star
//! on every machine.

use std::path::Path;
use std::time::Instant;

use universe_map::map::Map;
use universe_map::places::{Atlas, StarRef};

fn main() {
    let address = std::env::args().nth(1).unwrap_or_else(|| "Andromeda Galaxy/1".into());
    let file = Path::new(env!("CARGO_MANIFEST_DIR")).join("../map/universe-computed.json");
    let map = match Map::load(Some(&file), None) {
        Ok(m) => m,
        Err(e) => {
            eprintln!("{e}");
            std::process::exit(1);
        }
    };
    let Some(model) = map.star_field.clone() else {
        eprintln!("the map has no starField");
        std::process::exit(1);
    };
    let atlas = Atlas::new(model);
    println!("galaxies with stars of their own: {}", atlas.galaxies.len());
    for g in atlas.galaxies.iter().take(6) {
        println!("  {:28} {:>16} star systems, {:>12.0} light years away", g.name, g.systems, g.light_years());
    }
    let systems: f64 = atlas.galaxies.iter().map(|g| g.systems as f64).sum();
    println!("  ... in all {:.3e} star systems, none of them stored", systems);
    let a = match atlas.parse(&address) {
        Ok(a) => a,
        Err(e) => {
            eprintln!("{address}: {e}");
            std::process::exit(2);
        }
    };
    let StarRef::Number(n) = a.star else {
        eprintln!("{address}: a named star needs the star catalog; give a number");
        std::process::exit(2);
    };
    let clock = Instant::now();
    let star = atlas.star(a.galaxy, n);
    let took = clock.elapsed();
    println!("\n{}  (computed in {:.3} ms)", atlas.canonical(&a), took.as_secs_f64() * 1e3);
    println!("  {:.0} light years from the Sun", star.position_pc.length() * 3.261_563_777);
    println!("  star: {}", star.about.line());
    println!("  planets: {}", star.planets.len());
    for (i, p) in star.planets.iter().enumerate() {
        println!("    {}: {}", i + 1, p.line());
    }
}
