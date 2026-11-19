//! A galaxy not made yet, as the few numbers it is made from: its code.
//!
//! A galaxy of a survey is not stored as stars, nor as pictures. It is where it is and its redshift
//! (what the survey gives) and what follows from its number as a seed: its kind, its size, its
//! light, how round and how thick it is, how it is turned. That is the whole galaxy until someone
//! goes there. By analogy with the latent vector a generative model decodes a picture from, the
//! code is what the galaxy is decoded from:
//!
//!   code  ->  three prompts (prompt_y, prompt_x, prompt_z)
//!         ->  an image model draws its three views (the caller's part: this crate asks nobody)
//!         ->  xyz_make::make puts them on one grid
//!         ->  photo_volume and procedural compute its stars where the three views overlap
//!
//! The same number gives the same code on every machine. The pictures an image model draws from it
//! are not the same twice, so a galaxy made once is to be kept and shared, not made again.
//!
//! The SDSS file (Gaia Sky's catalog-sdss-12, 327,835 galaxies) has only where each one is, so its
//! kind follows its redshift as SDSS chose its targets: beyond z 0.35 its luminous red galaxies and
//! quasars, whose galaxies are giant ellipticals by now (each is shown where and as it is now);
//! nearer, its main sample, of spirals, ellipticals and irregulars about as common as in the local
//! universe. The kind, the size and the light of one galaxy are therefore drawn, not measured.

use glam::{DMat3, DVec3};

use crate::survey::SurveyGalaxy;
use crate::xyz_make::Spec;

#[derive(Clone, Debug)]
pub struct GalaxyCode {
    /// "SDSS <number>"
    pub name: String,
    /// The same as a file name: "SDSS_<number>".
    pub stem: String,
    /// Equatorial internal axes, parsecs from the Sun, where it is now (comoving).
    pub position_pc: DVec3,
    pub z: f64,
    pub seed: u64,
    /// "spiral", "elliptical", "irregular" or "giant elliptical"
    pub kind: &'static str,
    /// Half its width, kiloparsecs.
    pub half_kpc: f64,
    pub abs_mag: f64,
    /// How many stars it has, from its light (the Milky Way's 2.7e11 at -20.8).
    pub stars: f64,
    /// Short over long axis, seen from above.
    pub axis_ratio: f64,
    /// Height over width, seen from the side.
    pub thickness: f64,
}

fn mix(mut x: u64) -> u64 {
    x = x.wrapping_add(0x9E37_79B9_7F4A_7C15);
    x = (x ^ (x >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    x = (x ^ (x >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    x ^ (x >> 31)
}

impl GalaxyCode {
    pub fn new(g: &SurveyGalaxy) -> GalaxyCode {
        let seed = mix(g.id as u64);
        let u = |k: u64| (mix(seed ^ k) >> 11) as f64 / (1u64 << 53) as f64;
        let between = |k: u64, lo: f64, hi: f64| lo + (hi - lo) * u(k);
        let kind = if g.z >= 0.35 {
            "giant elliptical"
        } else {
            match u(1) {
                x if x < 0.6 => "spiral",
                x if x < 0.9 => "elliptical",
                _ => "irregular",
            }
        };
        let (half_kpc, abs_mag, axis_ratio, thickness) = match kind {
            "giant elliptical" => (between(2, 60.0, 80.0), between(3, -23.3, -22.5), between(4, 0.7, 0.95), between(5, 0.5, 0.8)),
            "elliptical" => (between(2, 30.0, 50.0), between(3, -22.0, -20.5), between(4, 0.6, 0.95), between(5, 0.5, 0.8)),
            "spiral" => (between(2, 25.0, 40.0), between(3, -21.5, -20.0), 0.95, 0.12),
            _ => (between(2, 8.0, 15.0), between(3, -18.5, -17.0), between(4, 0.6, 0.9), between(5, 0.4, 0.6)),
        };
        GalaxyCode {
            name: format!("SDSS {}", g.id),
            stem: format!("SDSS_{}", g.id),
            position_pc: g.position_pc,
            z: g.z,
            seed,
            kind,
            half_kpc,
            abs_mag,
            stars: 2.7e11 * 10f64.powf(-0.4 * (abs_mag + 20.8)),
            axis_ratio,
            thickness,
        }
    }

    /// The galaxy of a list by its name, "SDSS <number>".
    pub fn find(list: &[SurveyGalaxy], name: &str) -> Option<GalaxyCode> {
        let name = name.trim();
        let id: i64 = name.strip_prefix("SDSS").or_else(|| name.strip_prefix("sdss"))?.trim().trim_start_matches(['_', '#']).parse().ok()?;
        list.iter().find(|g| g.id == id).map(GalaxyCode::new)
    }

    /// Its axes (galactic rows: long, pole, short), turned as its seed says.
    pub fn axes(&self) -> DMat3 {
        let u = |k: u64| (mix(self.seed ^ k) >> 11) as f64 / (1u64 << 53) as f64;
        let q = glam::DQuat::from_xyzw(u(20) - 0.5, u(21) - 0.5, u(22) - 0.5, u(23) - 0.5).normalize();
        DMat3::from_quat(q)
    }

    pub fn about(&self) -> String {
        format!("made from its code: SDSS DR12 galaxy {} at z {:.3}, as it is now ({}); its pictures by an image model", self.name, self.z, self.kind)
    }

    /// What the image model is asked for: its view from above...
    pub fn prompt_y(&self) -> String {
        let body = match self.kind {
            "giant elliptical" => format!("A giant elliptical galaxy at the centre of a galaxy cluster, seen straight along its shortest axis: a smooth, featureless warm golden-white glow, axis ratio about {:.1}, a bright compact core fading smoothly outward (de Vaucouleurs profile), a faint extended halo, a few faint globular clusters as tiny points near it", self.axis_ratio),
            "elliptical" => format!("An elliptical galaxy seen straight along its shortest axis: a smooth, featureless golden glow, axis ratio about {:.1}, a bright core fading smoothly outward", self.axis_ratio),
            "spiral" => "A spiral galaxy seen exactly face-on: a yellowish bulge, blue spiral arms with pink star-forming regions and brown dust lanes along their inner edges".to_string(),
            _ => "An irregular galaxy seen from above: clumpy, blue, knots of young stars and pink gas, no clear arms".to_string(),
        };
        format!("{body}. {}", Self::COMMON)
    }

    /// ... and from the side along its long axis, drawn from the view from above (given to the
    /// image model as its input picture) so that it is the same galaxy...
    pub fn prompt_x(&self) -> String {
        self.prompt_side("its longest axis lies left to right")
    }

    /// ... and from the side along its short axis.
    pub fn prompt_z(&self) -> String {
        self.prompt_side("its intermediate axis lies left to right")
    }

    fn prompt_side(&self, how: &str) -> String {
        let side = match self.kind {
            "giant elliptical" | "elliptical" => format!("a triaxial elliptical, not a disc: no dust lane, no arms, no bright line along its middle; height over width about {:.1}", self.thickness),
            "spiral" => "a thin disc seen exactly edge-on with a dust lane along its middle and the bulge rising above and below it".to_string(),
            _ => format!("a thick clumpy body seen from the side, height over width about {:.1}", self.thickness),
        };
        format!("This is a galaxy seen from above. Draw the SAME galaxy seen from the side, exactly perpendicular to the first view ({how}): the same colours, the same core and light profile, the same extent left to right, as {side}. The same scale as the input. {}", Self::COMMON)
    }

    const COMMON: &'static str = "Astronomical photograph, exactly centred, the galaxy filling about 80% of the frame, pure black background, no foreground stars, no other galaxies, no text, no labels.";

    /// What xyz_make::make needs to make it from its three pictures (y.png, x.png, z.png).
    pub fn spec(&self) -> Spec {
        Spec {
            name: self.name.clone(),
            stem: self.stem.clone(),
            centre_pc: crate::places::galactic_to_equatorial().transpose() * self.position_pc,
            axes: self.axes(),
            half_kpc: self.half_kpc,
            stars: self.stars,
            abs_mag: self.abs_mag,
            about: self.about(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn galaxy(id: i64, z: f64) -> SurveyGalaxy {
        SurveyGalaxy { id, position_pc: DVec3::new(1.0e8, 2.0e8, 3.0e8), z }
    }

    #[test]
    fn the_same_number_is_the_same_galaxy() {
        let (a, b) = (GalaxyCode::new(&galaxy(75, 0.1)), GalaxyCode::new(&galaxy(75, 0.1)));
        assert_eq!(a.kind, b.kind);
        assert_eq!(a.half_kpc, b.half_kpc);
        assert_eq!(a.abs_mag, b.abs_mag);
        assert_eq!(a.axes(), b.axes());
        assert_eq!(a.prompt_y(), b.prompt_y());
    }

    #[test]
    fn galaxies_differ_and_their_kinds_are_as_common_as_said() {
        let codes: Vec<GalaxyCode> = (1..=10_000).map(|id| GalaxyCode::new(&galaxy(id, 0.1))).collect();
        let share = |kind: &str| codes.iter().filter(|c| c.kind == kind).count() as f64 / codes.len() as f64;
        assert!((share("spiral") - 0.6).abs() < 0.03, "spirals {}", share("spiral"));
        assert!((share("elliptical") - 0.3).abs() < 0.03, "ellipticals {}", share("elliptical"));
        assert!((share("irregular") - 0.1).abs() < 0.03, "irregulars {}", share("irregular"));
        let mut sizes: Vec<u64> = codes.iter().map(|c| c.half_kpc.to_bits()).collect();
        sizes.sort_unstable();
        sizes.dedup();
        assert!(sizes.len() > 9_990, "only {} sizes among 10000 galaxies", sizes.len());
        assert!(codes.iter().all(|c| c.stars > 1e9 && c.stars < 1e13));
    }

    #[test]
    fn far_galaxies_are_giant_ellipticals() {
        assert!((1..=200).all(|id| GalaxyCode::new(&galaxy(id, 0.5)).kind == "giant elliptical"));
    }

    #[test]
    fn found_by_name() {
        let list = vec![galaxy(7, 0.1), galaxy(75, 0.2)];
        assert_eq!(GalaxyCode::find(&list, "SDSS 75").unwrap().name, "SDSS 75");
        assert_eq!(GalaxyCode::find(&list, "sdss_7").unwrap().stem, "SDSS_7");
        assert!(GalaxyCode::find(&list, "SDSS 8").is_none());
        assert!(GalaxyCode::find(&list, "Andromeda").is_none());
    }

    #[test]
    fn its_axes_are_a_turn() {
        let m = GalaxyCode::new(&galaxy(75, 0.1)).axes();
        assert!((m.determinant() - 1.0).abs() < 1e-9);
        assert!((m * m.transpose() - DMat3::IDENTITY).to_cols_array().iter().all(|v| v.abs() < 1e-9));
    }
}
