//! A galaxy of the SDSS survey as its code, and what an image model is asked to draw it from:
//!
//!   cargo run --release --example galaxy_code -- "SDSS 75"
//!   cargo run --release --example galaxy_code -- "SDSS 75" path/to/sdss_dr12.bin
//!   cargo run --release --example galaxy_code -- "SDSS 75" path/to/sdss_dr12.bin pictures/
//!
//! The list of galaxies is Gaia Sky's catalog-sdss-12 (11.8 MB, 327,835 galaxies; not in this
//! repository): by default ~/.gaiasky/data/catalog-sdss-12/sdss/sdss_dr12.bin.
//!
//! With a folder holding y.png, x.png and z.png (the three views an image model drew from the
//! prompts printed here; x and z are drawn with y as the input picture), the galaxy is made there:
//! its three views on one grid, its disc's pictures and its entry (galaxy.json). This example asks
//! no image model itself.

use std::path::PathBuf;

use universe_map::galaxy_code::GalaxyCode;
use universe_map::{survey, xyz_make};

fn main() {
    let mut args = std::env::args().skip(1);
    let name = args.next().unwrap_or_else(|| "SDSS 75".into());
    let file = args.next().map(PathBuf::from).unwrap_or_else(|| {
        PathBuf::from(std::env::var("HOME").unwrap_or_default()).join(".gaiasky/data/catalog-sdss-12/sdss/sdss_dr12.bin")
    });
    let pictures = args.next().map(PathBuf::from);
    let list = match survey::read(&file) {
        Ok(l) => l,
        Err(e) => {
            eprintln!("{e}\nthe list of galaxies is Gaia Sky's catalog-sdss-12 (see the top of this file)");
            std::process::exit(1);
        }
    };
    println!("{} galaxies in {} ({:.1} MB, {:.0} bytes a galaxy)", list.len(), file.display(), size_mb(&file), size_mb(&file) * 1e6 / list.len() as f64);
    let Some(code) = GalaxyCode::find(&list, &name) else {
        eprintln!("{name}: no such galaxy (\"SDSS <number>\")");
        std::process::exit(2);
    };
    println!("\n{}", code.name);
    println!("  {:.0} million light years away (now), z {:.3}", code.position_pc.length() * 3.261_563_777 / 1e6, code.z);
    println!("  {}, {:.0} kpc across, absolute magnitude {:.1}, about {:.2e} stars", code.kind, code.half_kpc * 2.0, code.abs_mag, code.stars);
    println!("  axis ratio {:.2} from above, height over width {:.2} from the side", code.axis_ratio, code.thickness);
    println!("\nfrom above (y.png):\n  {}", code.prompt_y());
    println!("\nfrom the side, with y.png as the input (x.png):\n  {}", code.prompt_x());
    println!("\nfrom the other side, with y.png as the input (z.png):\n  {}", code.prompt_z());
    if let Some(dir) = pictures {
        match xyz_make::make(&dir, &code.spec(), 1024) {
            Ok(said) => println!("\nmade in {}: {said}", dir.display()),
            Err(e) => {
                eprintln!("\n{e}");
                std::process::exit(3);
            }
        }
    }
}

fn size_mb(file: &std::path::Path) -> f64 {
    std::fs::metadata(file).map(|m| m.len() as f64 / 1e6).unwrap_or(0.0)
}
