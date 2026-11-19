//! universe-map: a computed map of the universe, with no engine and no renderer.
//!
//! `map` reads `map/universe.json`; `places` finds a star or a planet by its address
//! ("Andromeda Galaxy/1/2"); `procedural`, `stars`, `systems`, `planets` and `rare` compute what
//! is there from the address alone; `fitted` makes the stars of a galaxy from its numbers;
//! `ground` computes the ground of a planet; `astro`, `vsop87` and `bodies` place the Solar
//! System; `travel` and `journey` find and fly to a place. Nothing here draws.
#![allow(dead_code)]

pub mod astro;
pub mod billboards;
pub mod bodies;
pub mod fitted;
pub mod ground;
pub mod journey;
pub mod map;
pub mod photo_volume;
pub mod places;
pub mod planets;
pub mod procedural;
pub mod rare;
pub mod stars;
pub mod survey;
pub mod systems;
pub mod travel;
pub mod vsop87;

/// What the modules ask of a renderer: only the turn from galactic to equatorial axes.
pub mod render {
    pub use crate::places::galactic_to_equatorial;
}

/// Where a trip goes: a body of the Solar System, a named star of the catalog (their indices
/// in the lists given to `travel::search`), or the far side of the Milky Way.
#[derive(Clone, Copy, Debug)]
pub enum Dest {
    Body(usize),
    Star(usize),
    FarSide,
}
