//! Reader for the Gaia Sky binary star catalog format (`particles_*.bin`).
//!
//! This mirrors `gaiasky.data.group.BinaryDataProvider` and `BinaryVersion3`/`BinaryVersion4`.
//! All numbers in the file are big-endian, as written by Java's `DataOutputStream`.

use std::fmt;
use std::fs;
use std::io;
use std::path::Path;

/// Internal length units per parsec. One internal unit is 1e9 m.
pub const PC_TO_U: f64 = 3.08567758149137e13 * 1000.0 * 1e-9;

/// One star of the catalog. Positions and velocities are in internal units, as in the Java code.
#[derive(Debug, Clone, PartialEq)]
pub struct Star {
    pub id: i64,
    /// All the names of the star. If the file has none, the id is used as the only name.
    pub names: Vec<String>,
    /// Equatorial cartesian position, in internal units.
    pub x: f64,
    pub y: f64,
    pub z: f64,
    /// Velocity vector, in internal units per year.
    pub vx: f32,
    pub vy: f32,
    pub vz: f32,
    /// Proper motions in mas/yr and radial velocity in km/s.
    pub mu_alpha: f32,
    pub mu_delta: f32,
    pub rad_vel: f32,
    pub app_mag: f32,
    pub abs_mag: f32,
    /// Color packed in a float, ABGR bits, as libGDX does.
    pub color: f32,
    pub size: f32,
    /// Hipparcos number, taken from the names. -1 if the star has none.
    pub hip: i32,
    /// Effective temperature in K.
    pub t_eff: f32,
    /// Surface gravity and metallicity. NaN in files of version 3.
    pub log_g: f32,
    pub mh: f32,
}

impl Star {
    /// Distance from the origin (the Sun), in parsecs.
    pub fn distance_pc(&self) -> f64 {
        (self.x * self.x + self.y * self.y + self.z * self.z).sqrt() / PC_TO_U
    }

    /// The color as red, green, blue and alpha in 0..=255.
    pub fn rgba(&self) -> [u8; 4] {
        let bits = self.color.to_bits();
        [bits as u8, (bits >> 8) as u8, (bits >> 16) as u8, (bits >> 24) as u8]
    }
}

#[derive(Debug)]
pub enum Error {
    Io(io::Error),
    /// The file ended in the middle of a value, at this byte offset.
    Truncated(usize),
    UnsupportedVersion(i32),
    BadCount(i32),
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Error::Io(e) => write!(f, "{e}"),
            Error::Truncated(at) => write!(f, "file ends unexpectedly at byte {at}"),
            Error::UnsupportedVersion(v) => write!(f, "unsupported catalog version {v} (3 and 4 are supported)"),
            Error::BadCount(n) => write!(f, "invalid star count {n}"),
        }
    }
}

impl std::error::Error for Error {}

impl From<io::Error> for Error {
    fn from(e: io::Error) -> Self {
        Error::Io(e)
    }
}

/// A star catalog file: its format version and its stars.
#[derive(Debug)]
pub struct Catalog {
    pub version: i32,
    pub stars: Vec<Star>,
}

struct Cursor<'a> {
    buf: &'a [u8],
    pos: usize,
}

impl<'a> Cursor<'a> {
    fn take<const N: usize>(&mut self) -> Result<[u8; N], Error> {
        let end = self.pos + N;
        let bytes = self.buf.get(self.pos..end).ok_or(Error::Truncated(self.pos))?;
        self.pos = end;
        Ok(bytes.try_into().unwrap())
    }
    fn i32(&mut self) -> Result<i32, Error> {
        Ok(i32::from_be_bytes(self.take()?))
    }
    fn i64(&mut self) -> Result<i64, Error> {
        Ok(i64::from_be_bytes(self.take()?))
    }
    fn f32(&mut self) -> Result<f32, Error> {
        Ok(f32::from_be_bytes(self.take()?))
    }
    fn f64(&mut self) -> Result<f64, Error> {
        Ok(f64::from_be_bytes(self.take()?))
    }
    /// A string of `len` UTF-16 code units, as written by `DataOutputStream.writeChars`.
    fn utf16(&mut self, len: usize) -> Result<String, Error> {
        let end = self.pos + len * 2;
        let bytes = self.buf.get(self.pos..end).ok_or(Error::Truncated(self.pos))?;
        self.pos = end;
        let units = bytes.chunks_exact(2).map(|c| u16::from_be_bytes([c[0], c[1]]));
        Ok(char::decode_utf16(units).map(|r| r.unwrap_or(char::REPLACEMENT_CHARACTER)).collect())
    }
}

/// Reads a catalog from the bytes of a particle file.
pub fn parse(buf: &[u8]) -> Result<Catalog, Error> {
    let mut c = Cursor { buf, pos: 0 };
    // Files of version 2 and later start with a negative token and the version number.
    let token = c.i32()?;
    if token >= 0 {
        return Err(Error::UnsupportedVersion(1));
    }
    let version = c.i32()?;
    // Version 3 has 11 floats per star. Version 4 adds log g and metallicity.
    let n_floats = match version {
        3 => 11,
        4 => 13,
        v => return Err(Error::UnsupportedVersion(v)),
    };
    let count = c.i32()?;
    if count < 0 {
        return Err(Error::BadCount(count));
    }

    let mut stars = Vec::with_capacity(count as usize);
    for _ in 0..count {
        let (x, y, z) = (c.f64()?, c.f64()?, c.f64()?);
        let mut f = [f32::NAN; 13];
        for v in f.iter_mut().take(n_floats) {
            *v = c.f32()?;
        }
        let id = c.i64()?;
        let name_len = c.i32()?;
        let names: Vec<String> = if name_len <= 0 {
            vec![id.to_string()]
        } else {
            c.utf16(name_len as usize)?.split('|').map(str::to_owned).collect()
        };
        let hip = names
            .iter()
            .find(|n| n.starts_with("HIP "))
            .and_then(|n| n.trim()[4..].parse().ok())
            .unwrap_or(-1);
        stars.push(Star {
            id,
            names,
            x,
            y,
            z,
            vx: f[0],
            vy: f[1],
            vz: f[2],
            mu_alpha: f[3],
            mu_delta: f[4],
            rad_vel: f[5],
            app_mag: f[6],
            abs_mag: f[7],
            color: f[8],
            size: f[9],
            hip,
            t_eff: f[10],
            log_g: f[11],
            mh: f[12],
        });
    }
    Ok(Catalog { version, stars })
}

/// Reads a catalog from a particle file.
pub fn load(path: impl AsRef<Path>) -> Result<Catalog, Error> {
    parse(&fs::read(path)?)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample(version: i32) -> Vec<u8> {
        let mut b = Vec::new();
        b.extend((-1i32).to_be_bytes());
        b.extend(version.to_be_bytes());
        b.extend(2i32.to_be_bytes());
        let n_floats = if version == 3 { 11 } else { 13 };
        for (id, name) in [(7i64, "Sirius|HIP 32349"), (9i64, "")] {
            for v in [PC_TO_U * 3.0, 0.0, PC_TO_U * 4.0] {
                b.extend(v.to_be_bytes());
            }
            for i in 0..n_floats {
                b.extend((i as f32).to_be_bytes());
            }
            b.extend(id.to_be_bytes());
            let units: Vec<u16> = name.encode_utf16().collect();
            b.extend((units.len() as i32).to_be_bytes());
            for u in units {
                b.extend(u.to_be_bytes());
            }
        }
        b
    }

    #[test]
    fn reads_version_3() {
        let cat = parse(&sample(3)).unwrap();
        assert_eq!(cat.version, 3);
        assert_eq!(cat.stars.len(), 2);
        let s = &cat.stars[0];
        assert_eq!(s.names, ["Sirius", "HIP 32349"]);
        assert_eq!(s.hip, 32349);
        assert!((s.distance_pc() - 5.0).abs() < 1e-9);
        assert_eq!((s.vx, s.app_mag, s.abs_mag, s.size, s.t_eff), (0.0, 6.0, 7.0, 9.0, 10.0));
        assert!(s.log_g.is_nan());
        let t = &cat.stars[1];
        assert_eq!(t.names, ["9"]);
        assert_eq!(t.hip, -1);
    }

    #[test]
    fn reads_version_4() {
        let cat = parse(&sample(4)).unwrap();
        assert_eq!((cat.stars[1].log_g, cat.stars[1].mh), (11.0, 12.0));
    }

    #[test]
    fn rejects_truncated_and_unknown() {
        let b = sample(3);
        assert!(matches!(parse(&b[..b.len() - 3]), Err(Error::Truncated(_))));
        assert!(matches!(parse(&sample(4)[..8]), Err(Error::Truncated(8))));
        let mut bad = sample(3);
        bad[7] = 9;
        assert!(matches!(parse(&bad), Err(Error::UnsupportedVersion(9))));
    }
}
