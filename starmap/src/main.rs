//! Prints a summary of a Gaia Sky star catalog file. Usage: starmap <particles.bin> [star name]

use std::env;
use std::process::ExitCode;
use std::time::Instant;

fn main() -> ExitCode {
    let args: Vec<String> = env::args().collect();
    let Some(path) = args.get(1) else {
        eprintln!("usage: starmap <particles.bin> [star name]");
        return ExitCode::from(2);
    };
    let start = Instant::now();
    let cat = match starmap::load(path) {
        Ok(c) => c,
        Err(e) => {
            eprintln!("{path}: {e}");
            return ExitCode::FAILURE;
        }
    };
    println!("version {}, {} stars, read in {:.0?}", cat.version, cat.stars.len(), start.elapsed());

    let show = |s: &starmap::Star| {
        let [r, g, b, _] = s.rgba();
        println!(
            "{:<22} mag {:>6.2}  {:>10.3} pc  {:>6.0} K  rgb({r},{g},{b})  x={:.4e} y={:.4e} z={:.4e}",
            s.names[0], s.app_mag, s.distance_pc(), s.t_eff, s.x, s.y, s.z
        );
    };

    if let Some(name) = args.get(2) {
        let found: Vec<_> = cat.stars.iter().filter(|s| s.names.iter().any(|n| n.eq_ignore_ascii_case(name))).collect();
        if found.is_empty() {
            eprintln!("no star named {name}");
            return ExitCode::FAILURE;
        }
        for s in found {
            show(s);
            println!("  names: {}", s.names.join(", "));
        }
    } else {
        let mut brightest: Vec<_> = cat.stars.iter().filter(|s| s.app_mag.is_finite()).collect();
        brightest.sort_by(|a, b| a.app_mag.total_cmp(&b.app_mag));
        println!("brightest:");
        for s in brightest.iter().take(10) {
            show(s);
        }
    }
    ExitCode::SUCCESS
}
