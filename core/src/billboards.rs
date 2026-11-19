//! Galaxies and nebulae that are drawn as a picture in space, from the catalog JSON files.

use std::path::{Path, PathBuf};

use glam::DVec3;
use serde_json::Value;

pub struct Billboard {
    pub name: String,
    /// Position in parsecs, in the internal equatorial axes.
    pub position_pc: DVec3,
    /// Side of the picture in parsecs.
    pub size_pc: f64,
    pub color: [f32; 4],
    /// Distances in parsecs between which the picture fades in, and the opacity at each end.
    pub fade_in: Option<([f64; 2], [f64; 2])>,
    pub texture: PathBuf,
}

/// The JSON files of the catalogs are written by hand and some have commas before a closing
/// bracket. This removes them, outside of strings.
fn strip_trailing_commas(text: &str) -> String {
    let mut out = String::with_capacity(text.len());
    let mut in_string = false;
    let mut escaped = false;
    for c in text.chars() {
        if in_string {
            out.push(c);
            if escaped {
                escaped = false;
            } else if c == '\\' {
                escaped = true;
            } else if c == '"' {
                in_string = false;
            }
            continue;
        }
        if c == '"' {
            in_string = true;
        } else if c == ']' || c == '}' {
            // Drop a comma that only has white space between it and this bracket.
            let trimmed = out.trim_end().len();
            if out[..trimmed].ends_with(',') {
                let tail = out.split_off(trimmed);
                out.pop();
                out.push_str(&tail);
            }
        }
        out.push(c);
    }
    out
}

fn pair(v: Option<&Value>) -> Option<[f64; 2]> {
    let a = v?.as_array()?;
    Some([a.first()?.as_f64()?, a.get(1)?.as_f64()?])
}

fn billboard(v: &Value, data: &Path) -> Option<Billboard> {
    let archetype = v.get("archetype")?.as_str()?;
    if archetype != "Billboard" && archetype != "BillboardGalaxy" {
        return None;
    }
    let name = match v.get("names").and_then(|n| n.get(0)) {
        Some(n) => n.as_str()?,
        None => v.get("name")?.as_str()?,
    };
    let eq = v.get("coordinates")?.get("equatorial")?.as_array()?;
    let (ra, dec, dist) = (eq.first()?.as_f64()?.to_radians(), eq.get(1)?.as_f64()?.to_radians(), eq.get(2)?.as_f64()?);
    let material = v.get("model")?.get("material")?;
    let spec = material.get("diffuse").or_else(|| material.get("base"))?.as_str()?;
    let texture = data.join(spec.strip_prefix("$data/")?);
    let color = v.get("color").and_then(Value::as_array);
    let channel = |i: usize| color.and_then(|c| c.get(i)).and_then(Value::as_f64).unwrap_or(1.0) as f32;
    let fade_in = pair(v.get("fadeIn")).map(|d| (d, pair(v.get("fadeInMap")).unwrap_or([0.0, 1.0])));
    Some(Billboard {
        name: name.to_owned(),
        position_pc: DVec3::new(dec.cos() * ra.sin(), dec.sin(), dec.cos() * ra.cos()) * dist,
        size_pc: v.get("sizePc")?.as_f64()?,
        color: [channel(0), channel(1), channel(2), channel(3)],
        fade_in,
        texture: texture.exists().then_some(texture)?,
    })
}

fn read(path: &Path, data: &Path, out: &mut Vec<Billboard>) -> Result<(), String> {
    let text = std::fs::read_to_string(path).map_err(|e| format!("{}: {e}", path.display()))?;
    let v: Value = serde_json::from_str(&strip_trailing_commas(&text)).map_err(|e| format!("{}: {e}", path.display()))?;
    let objects = v.get("objects").and_then(Value::as_array).ok_or(format!("{}: no \"objects\" array", path.display()))?;
    out.extend(objects.iter().filter_map(|o| billboard(o, data)));
    Ok(())
}

/// Reads the galaxies and nebulae of the catalogs of the map. `data` is the folder that `$data`
/// stands for in the texture paths.
pub fn load(catalogs: &[std::path::PathBuf], data: &Path) -> Result<Vec<Billboard>, String> {
    let mut out = Vec::new();
    for path in catalogs {
        read(path, data, &mut out)?;
    }
    Ok(out)
}

impl Billboard {
    /// Opacity seen from a distance in parsecs.
    pub fn alpha(&self, distance_pc: f64) -> f32 {
        let fade = match self.fade_in {
            Some(([near, far], [a, b])) if far > near => a + (b - a) * ((distance_pc - near) / (far - near)).clamp(0.0, 1.0),
            _ => 1.0,
        };
        self.color[3] * fade as f32
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn strips_commas() {
        assert_eq!(strip_trailing_commas(r#"{"a": [1, 2, ], "b": "x,]", }"#), r#"{"a": [1, 2 ], "b": "x,]" }"#);
        assert_eq!(strip_trailing_commas("[{\"a\": 1},\n  ]"), "[{\"a\": 1}\n  ]");
        assert_eq!(strip_trailing_commas(r#"["\",]"]"#), r#"["\",]"]"#);
    }

    #[test]
    fn reads_the_catalogs() {
        let Some(home) = std::env::var_os("HOME") else { return };
        let data = Path::new(&home).join(".gaiasky/data");
        // The data is not part of the repository. Without it there is nothing to check.
        if !data.join("catalog-nebulae/particles-nebulae.json").exists() || !data.join("catalog-nbg/particles-nbg.json").exists() {
            return;
        }
        let all = load(&[data.join("catalog-nbg/particles-nbg.json"), data.join("catalog-nebulae/particles-nebulae.json")], &data).unwrap();
        assert_eq!(all.len(), 55, "14 galaxies and 41 nebulae");
        let m33 = all.iter().find(|b| b.name == "Triangulum spiral").unwrap();
        assert!((m33.position_pc.length() - 840_000.0).abs() < 1.0);
        // Right ascension 23.46°, declination 30.66°.
        assert!((m33.position_pc.y / 840_000.0 - 30.6603f64.to_radians().sin()).abs() < 1e-6);
        assert!((m33.position_pc.x.atan2(m33.position_pc.z).to_degrees() - 23.4617).abs() < 1e-6);
        assert_eq!(m33.alpha(840_000.0), 1.0);
        assert_eq!(m33.alpha(0.0), 0.5);
        let crab = all.iter().find(|b| b.name == "Crab Nebula BB").unwrap();
        assert_eq!(crab.alpha(1.0), 0.0, "the picture gives way to the volume up close");
        assert_eq!(crab.alpha(200.0), 0.5);
    }
}
