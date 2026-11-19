//! The galaxies of a redshift survey: SDSS DR12 as Gaia Sky ships it
//! (`$data/catalog-sdss-12/sdss/sdss_dr12.bin`, 327 835 galaxies), for drawing the large-scale
//! structure of the universe: the walls and filaments of galaxies round empty voids.
//!
//! The file (gaiasky.data.group.BinaryPointDataProvider, big endian): an int count, a bool
//! "extra" (false here), then per galaxy a long id, an int name length, that many UTF-16 chars,
//! and doubles right ascension (deg), declination (deg) and distance (pc).
//!
//! Its distances are the redshift through Hubble's law, d = cz / H0: they run to 29.6 Gpc, more
//! than the radius of the observable universe (about 14.3 Gpc comoving), for quasars at z ~ 7.
//! Here each is turned back into its redshift and then into the comoving distance of a flat
//! Lambda-CDM universe, which is where the galaxy is now. The H0 the file was made with is not
//! written in it; FILE_H0 is assumed (the usual round value). Below z ~ 0.1, where most of the
//! survey's galaxies are, the two distances differ by a few percent.

use glam::DVec3;

/// Speed of light, km/s.
const C_KM_S: f64 = 299_792.458;
/// The Hubble constant the file's distances were made with (assumed), km/s/Mpc.
pub const FILE_H0: f64 = 70.0;
/// The universe the comoving distances are computed in: Planck 2018 (H0 67.4, matter 0.315),
/// flat, radiation left out (it matters only far beyond these redshifts).
pub const H0: f64 = 67.4;
pub const OMEGA_M: f64 = 0.315;

/// One galaxy of the survey: where it is now (equatorial internal axes, as the catalog's stars:
/// x = cos dec sin ra, y = sin dec, z = cos dec cos ra; parsecs from the Sun) and its redshift.
#[derive(Clone, Copy, Debug)]
pub struct SurveyGalaxy {
    pub position_pc: DVec3,
    pub z: f64,
}

/// The comoving distance to redshift `z`, parsecs: (c / H0) times the integral of 1 / E(z) from 0
/// to z, E = sqrt(matter (1 + z)^3 + (1 - matter)), by Simpson's rule.
pub fn comoving_pc(z: f64) -> f64 {
    if z <= 0.0 {
        return 0.0;
    }
    let e = |x: f64| 1.0 / (OMEGA_M * (1.0 + x).powi(3) + (1.0 - OMEGA_M)).sqrt();
    let n = 200; // even
    let h = z / n as f64;
    let mut s = e(0.0) + e(z);
    for i in 1..n {
        s += e(i as f64 * h) * if i % 2 == 1 { 4.0 } else { 2.0 };
    }
    C_KM_S / H0 * (s * h / 3.0) * 1e6
}

/// Reads the survey file: every galaxy, with its distance made comoving (see above).
pub fn read(path: &std::path::Path) -> Result<Vec<SurveyGalaxy>, String> {
    let raw = std::fs::read(path).map_err(|e| format!("{}: {e}", path.display()))?;
    let be_i32 = |o: usize| -> Result<i32, String> {
        raw.get(o..o + 4).map(|b| i32::from_be_bytes(b.try_into().unwrap())).ok_or_else(|| "cut short".to_string())
    };
    let be_f64 = |o: usize| -> Result<f64, String> {
        raw.get(o..o + 8).map(|b| f64::from_be_bytes(b.try_into().unwrap())).ok_or_else(|| "cut short".to_string())
    };
    let n = be_i32(0)?;
    if n < 0 || raw.len() < 5 {
        return Err("not a particle file".into());
    }
    if raw[4] != 0 {
        return Err("the file has extra fields (magnitudes, colours): not read here".into());
    }
    let mut out = Vec::with_capacity(n as usize);
    let mut o = 5;
    for _ in 0..n {
        let name_len = be_i32(o + 8)?;
        if name_len < 0 {
            return Err("bad name length".into());
        }
        o += 12 + 2 * name_len as usize;
        let (ra, dec, dist_pc) = (be_f64(o)?.to_radians(), be_f64(o + 8)?.to_radians(), be_f64(o + 16)?);
        o += 24;
        let z = dist_pc / 1e6 * FILE_H0 / C_KM_S;
        let d = comoving_pc(z);
        out.push(SurveyGalaxy { position_pc: DVec3::new(dec.cos() * ra.sin(), dec.sin(), dec.cos() * ra.cos()) * d, z });
    }
    if o != raw.len() {
        return Err(format!("{} bytes left over after {n} galaxies", raw.len() - o));
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn comoving_distances_are_the_known_ones() {
        // flat LCDM, H0 67.4, matter 0.315, by scipy's adaptive quadrature of the same integral:
        // z 0.1 -> 434.11 Mpc, z 1 -> 3401.26 Mpc, z 7 -> 8810.0 Mpc
        let mpc = |z: f64| comoving_pc(z) / 1e6;
        assert!((mpc(0.1) - 434.11).abs() < 0.05, "{}", mpc(0.1));
        assert!((mpc(1.0) - 3401.26).abs() < 0.5, "{}", mpc(1.0));
        assert!((mpc(7.0) - 8810.0).abs() < 2.0, "{}", mpc(7.0));
        // near: Hubble's law
        assert!((mpc(0.01) - 0.01 * C_KM_S / H0).abs() / mpc(0.01) < 0.01);
    }

    /// The file as Gaia Sky ships it, when it is in ~/.gaiasky/data (skipped otherwise).
    #[test]
    fn reads_the_sdss_file() {
        let Some(home) = std::env::var_os("HOME") else { return };
        let path = std::path::Path::new(&home).join(".gaiasky/data/catalog-sdss-12/sdss/sdss_dr12.bin");
        if !path.exists() {
            eprintln!("no {}: skipped", path.display());
            return;
        }
        let g = read(&path).unwrap();
        assert_eq!(g.len(), 327_835);
        let mut d: Vec<f64> = g.iter().map(|x| x.position_pc.length() / 1e6).collect();
        d.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let median = d[d.len() / 2];
        let max = d[d.len() - 1];
        eprintln!("SDSS: {} galaxies, comoving Mpc median {median:.0}, max {max:.0}", g.len());
        // the file's own median is 344 Mpc (z 0.080): comoving a little more with H0 67.4
        assert!(median > 330.0 && median < 370.0, "{median}");
        // nothing beyond the observable universe any more (the file had 29.6 Gpc)
        assert!(max < 14_300.0, "{max}");
        // directions on the unit sphere times the distance
        assert!(g.iter().all(|x| x.position_pc.is_finite()));
    }
}
