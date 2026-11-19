//! A timed trip of the camera from where it is to a body.
//!
//! The speed profile is the one of `gaiasky.util.camera.JourneyProfile`: the distance covered
//! grows exponentially away from the departure and shrinks the same way into the arrival, so both
//! ends are seen at a sensible pace however far apart they are.

use glam::{DQuat, DVec3};

use crate::astro;

/// Share of the time spent easing in, and the same again easing out.
const RAMP: f64 = 0.1;
/// Seconds the camera takes to turn towards the destination.
const TURN_SECONDS: f64 = 4.0;

fn ease(alpha: f64) -> f64 {
    let a = alpha.clamp(0.0, 1.0);
    let norm = 2.0 * RAMP * (1.0 - RAMP);
    if a < RAMP {
        a * a / norm
    } else if a > 1.0 - RAMP {
        let b = 1.0 - a;
        1.0 - b * b / norm
    } else {
        (a - 0.5 * RAMP) / (1.0 - RAMP)
    }
}

fn grow(u: f64, ratio: f64) -> f64 {
    if ratio <= 1.0 + 1e-9 {
        u
    } else {
        (u * ratio.ln()).exp_m1() / (ratio - 1.0)
    }
}

fn ratio(half: f64, scale: f64) -> f64 {
    if scale > 0.0 && half > 0.0 {
        (half / scale).max(1.0)
    } else {
        1.0
    }
}

/// Distance covered at the fraction `alpha` of the time. The scales are the sizes of what is at
/// each end: the trip moves about that much per moment there.
pub fn covered(alpha: f64, distance: f64, departure_scale: f64, arrival_scale: f64) -> f64 {
    let e = ease(alpha);
    let half = 0.5 * distance;
    if e <= 0.5 {
        half * grow(2.0 * e, ratio(half, departure_scale))
    } else {
        distance - half * grow(2.0 * (1.0 - e), ratio(half, arrival_scale))
    }
}

pub struct Journey {
    /// The direction that stays up on the screen.
    pub up: DVec3,
    /// Its whole time, slowed part and all (`slow_within`).
    pub duration: f64,
    pub elapsed: f64,
    start: DVec3,
    start_orientation: DQuat,
    departure_scale: f64,
    arrival_scale: f64,
    /// The time of the profile (without the slowed part's extra).
    base: f64,
    /// Slowed within distances of the end (into a galaxy): from each (clock time, profile time)
    /// on, the profile goes at `rate` of the clock; before the first, as the clock.
    slow: Vec<(f64, f64, f64)>,
    /// ... or slowed smoothly by how much is left (`slow_by_left`): the clock's and the profile's
    /// times side by side, for both ways.
    warp: Vec<(f64, f64)>,
}

/// Where the camera is and where it looks at a moment of the trip.
pub struct Step {
    pub position: DVec3,
    pub orientation: DQuat,
    pub remaining_km: f64,
    pub done: bool,
}

impl Journey {
    pub fn new(duration: f64, start: DVec3, start_orientation: DQuat, departure_scale: f64, arrival_scale: f64) -> Journey {
        let duration = duration.max(1e-3);
        Journey { up: DVec3::Y, duration, elapsed: 0.0, start, start_orientation, departure_scale, arrival_scale, base: duration, slow: Vec::new(), warp: Vec::new() }
    }

    /// Goes the last `zone` (km) before `end` more slowly: the same way, `extra_s` seconds longer
    /// (into a galaxy: its stars near the ship are made, and seen, before they are passed).
    pub fn slow_within(&mut self, zone: f64, extra_s: f64, end: DVec3) {
        let Some(from) = self.profile_time_left(zone, end) else { return };
        let inside = self.base - from;
        if inside <= 0.0 || extra_s <= 0.0 {
            return;
        }
        self.slow_zones(&[(zone, inside + extra_s)], end);
    }

    /// The profile's time when what is left of the trip to `end` is `left` (km); None if it
    /// never is that far.
    fn profile_time_left(&self, left: f64, end: DVec3) -> Option<f64> {
        let distance = (end - self.start).length();
        if left <= 0.0 || distance <= left {
            return None;
        }
        let rest = |t: f64| distance - covered(t / self.base, distance, self.departure_scale, self.arrival_scale);
        let (mut lo, mut hi) = (0.0, self.base);
        for _ in 0..80 {
            let mid = 0.5 * (lo + hi);
            if rest(mid) > left {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        Some(hi)
    }

    /// Goes stretches before `end` at times of their own: each (km left at its start, seconds
    /// it takes), the farthest first, each to the next (the last to the end). The trip is the
    /// same way, longer or shorter by as much. Into a galaxy: slowed a lot where it comes into
    /// its stars, and quicker again inside.
    pub fn slow_zones(&mut self, zones: &[(f64, f64)], end: DVec3) {
        let starts: Vec<(f64, f64)> = zones.iter().filter_map(|&(left, secs)| self.profile_time_left(left, end).map(|p| (p, secs.max(1e-3)))).collect();
        if starts.is_empty() {
            return;
        }
        let mut clock = starts[0].0;
        self.slow.clear();
        for (i, &(p, secs)) in starts.iter().enumerate() {
            let next = starts.get(i + 1).map_or(self.base, |n| n.0);
            if next <= p {
                continue;
            }
            self.slow.push((clock, p, (next - p) / secs));
            clock += secs;
        }
        self.duration = clock;
    }

    /// Goes slower, smoothly, by how much is left of the trip to `end`: at each (km left, rate)
    /// the profile goes at `rate` of the clock, between them evenly in the logarithms of both, as
    /// the first before it and as the last after it (rates of 1 at both ends: as it was there).
    /// Into a galaxy: braking long before it comes into its stars, slow there, quicker inside,
    /// never a jump in speed.
    pub fn slow_by_left(&mut self, knots: &[(f64, f64)], end: DVec3) {
        if knots.is_empty() {
            return;
        }
        let distance = (end - self.start).length();
        let rate_at = |left: f64| -> f64 {
            let l = left.max(1e-9).ln();
            let first = knots[0];
            if l >= first.0.ln() {
                return first.1;
            }
            for w in knots.windows(2) {
                let (a, b) = (w[0], w[1]);
                if l >= b.0.ln() {
                    let t = (l - a.0.ln()) / (b.0.ln() - a.0.ln());
                    // (smoothed at the knots: no corner in the speed)
                    let t = t * t * (3.0 - 2.0 * t);
                    return (a.1.ln() + (b.1.ln() - a.1.ln()) * t).exp();
                }
            }
            knots[knots.len() - 1].1
        };
        const N: usize = 20_000;
        self.warp.clear();
        self.warp.push((0.0, 0.0));
        let mut clock = 0.0;
        let mut last_left = distance;
        for i in 1..=N {
            // finer near the end, where it all happens: the profile's time as 1 - (1 - u)^3
            let u = i as f64 / N as f64;
            let tau = self.base * (1.0 - (1.0 - u).powi(3));
            let prev = self.warp[i - 1].1;
            let left = distance - covered(tau / self.base, distance, self.departure_scale, self.arrival_scale);
            let rate = rate_at(0.5 * (left + last_left));
            last_left = left;
            clock += (tau - prev) / rate.max(1e-6);
            self.warp.push((clock, tau));
        }
        self.slow.clear();
        self.duration = clock;
    }

    /// The trip's time (seconds from its start) when what is left of it to `end` is `left` (km).
    pub fn time_left(&self, left: f64, end: DVec3) -> Option<f64> {
        let p = self.profile_time_left(left, end)?;
        if !self.warp.is_empty() {
            let i = self.warp.partition_point(|w| w.1 < p).clamp(1, self.warp.len() - 1);
            let (a, b) = (self.warp[i - 1], self.warp[i]);
            return Some(a.0 + (b.0 - a.0) * (p - a.1) / (b.1 - a.1).max(1e-12));
        }
        // the clock at a profile time: back through the stretches
        let mut t = p;
        for &(clock, from, rate) in &self.slow {
            if p >= from {
                t = clock + (p - from) / rate;
            }
        }
        Some(t)
    }

    /// The profile's time at a time of the trip (slowed within the zones).
    fn profile_time(&self, elapsed: f64) -> f64 {
        let t = elapsed.clamp(0.0, self.duration);
        if !self.warp.is_empty() {
            let i = self.warp.partition_point(|w| w.0 < t).clamp(1, self.warp.len() - 1);
            let (a, b) = (self.warp[i - 1], self.warp[i]);
            return (a.1 + (b.1 - a.1) * (t - a.0) / (b.0 - a.0).max(1e-12)).min(self.base);
        }
        match self.slow.iter().rev().find(|z| t >= z.0) {
            Some(&(clock, from, rate)) => (from + (t - clock) * rate).min(self.base),
            None => t.min(self.base),
        }
    }

    /// Where the trip is at a time of it (seconds from its start), not moving it: for what is
    /// ahead (the stars round where it will be, made before it gets there).
    pub fn position_at(&self, elapsed: f64, end: DVec3) -> DVec3 {
        let alpha = self.profile_time(elapsed) / self.base;
        let leg = end - self.start;
        let distance = leg.length();
        let gone = covered(alpha, distance, self.departure_scale, self.arrival_scale);
        if distance > 0.0 { self.start + leg * (gone / distance) } else { end }
    }

    /// Moves the trip on by `dt` seconds. `end` is where the trip stops and `centre` the point
    /// the camera looks at.
    pub fn advance(&mut self, dt: f64, end: DVec3, centre: DVec3) -> Step {
        self.elapsed = (self.elapsed + dt).min(self.duration);
        let alpha = self.profile_time(self.elapsed) / self.base;
        let leg = end - self.start;
        let distance = leg.length();
        let gone = covered(alpha, distance, self.departure_scale, self.arrival_scale);
        let position = if distance > 0.0 { self.start + leg * (gone / distance) } else { end };
        // Turn to the destination first, then keep it in the middle.
        let facing = astro::look_at(centre - position, self.up);
        let turn = (self.elapsed / TURN_SECONDS.min(self.duration * 0.5)).clamp(0.0, 1.0);
        let smooth = turn * turn * (3.0 - 2.0 * turn);
        Step { position, orientation: self.start_orientation.slerp(facing, smooth).normalize(), remaining_km: distance - gone, done: self.elapsed >= self.duration }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn covers_the_whole_distance_without_going_back() {
        let (d, a, b) = (5.9e9, 6371.0, 1190.0);
        assert_eq!(covered(0.0, d, a, b), 0.0);
        assert!((covered(1.0, d, a, b) - d).abs() < 1e-3);
        assert!((covered(0.5, d, a, b) - d / 2.0).abs() < 1.0);
        let mut last = 0.0;
        for i in 1..=1000 {
            let c = covered(i as f64 / 1000.0, d, a, b);
            assert!(c >= last, "went back at {i}");
            last = c;
        }
        // Slow next to the planets: after a tenth of the time the Earth is still close.
        assert!(covered(0.1, d, a, b) < 100.0 * a);
        assert!(d - covered(0.9, d, a, b) < 100.0 * b);
    }

    /// Slowed into a galaxy: the same way, longer; within the zone slower by as much, outside as
    /// before; it still arrives exactly; and where it will be is where it gets.
    #[test]
    fn slowed_near_the_end_it_still_arrives() {
        let km_ly = 9.4607e12;
        let (start, end) = (DVec3::ZERO, DVec3::new(2.5e6 * km_ly, 0.0, 0.0));
        let plain = Journey::new(60.0, start, DQuat::IDENTITY, 6371.0, 1.7e6);
        let mut slow = Journey::new(60.0, start, DQuat::IDENTITY, 6371.0, 1.7e6);
        slow.slow_within(30_000.0 * km_ly, 40.0, end);
        assert!((slow.duration - 100.0).abs() < 1e-9);
        let speed = |j: &Journey, t: f64| (j.position_at(t + 0.01, end) - j.position_at(t, end)).length() / 0.01 / km_ly;
        let (_, from, rate) = slow.slow[0];
        // before the zone the same; inside it slower by `rate`
        assert!((speed(&plain, from * 0.9) - speed(&slow, from * 0.9)).abs() < 1e-6 * speed(&plain, from * 0.9));
        let t_plain = from + 5.0;
        let t_slow = from + 5.0 / rate;
        assert!((speed(&slow, t_slow) / speed(&plain, t_plain) - rate).abs() < 0.01, "{} {}", speed(&slow, t_slow), speed(&plain, t_plain));
        let left = (end - slow.position_at(from, end)).length() / km_ly;
        eprintln!("into the zone at {from:.1} s ({left:.0} ly left); inside {:.1} s instead of {:.1}, {:.2} as fast", 100.0 - from, 60.0 - from, rate);
        assert!((left - 30_000.0).abs() < 1.0);
        let mut step = slow.advance(0.0, end, end);
        let ahead = slow.position_at(70.0, end);
        for _ in 0..(70 * 60) {
            step = slow.advance(1.0 / 60.0, end, end);
        }
        assert!((step.position - ahead).length() < 1.0 + 1e-9 * ahead.length(), "where it would be is where it is");
        for _ in 0..(31 * 60) {
            step = slow.advance(1.0 / 60.0, end, end);
        }
        assert!(step.done && step.position == end);
    }

    #[test]
    fn arrives_exactly_and_looks_at_the_destination() {
        let (start, end, centre) = (DVec3::new(1e8, 0.0, 0.0), DVec3::new(-4e9, 1e9, 2e9), DVec3::new(-4e9, 1e9, 2.000_005e9));
        let mut j = Journey::new(60.0, start, DQuat::IDENTITY, 6371.0, 1190.0);
        let mut step = j.advance(0.0, end, centre);
        assert_eq!(step.position, start);
        assert!(!step.done);
        for _ in 0..3601 {
            step = j.advance(1.0 / 60.0, end, centre);
        }
        assert!(step.done);
        assert_eq!(step.position, end);
        assert_eq!(step.remaining_km, 0.0);
        assert!(((step.orientation * DVec3::NEG_Z) - (centre - end).normalize()).length() < 1e-9);
    }

    /// Into a galaxy: slowed a lot where the ship comes into its stars, quicker again inside;
    /// it still arrives exactly, and the time it comes in is when it is there.
    #[test]
    fn slowed_where_it_comes_in_quicker_inside() {
        let km_ly = 9.4607e12;
        let (start, end) = (DVec3::ZERO, DVec3::new(2.5e6 * km_ly, 0.0, 0.0));
        let plain = Journey::new(60.0, start, DQuat::IDENTITY, 6371.0, 1.7e6);
        let mut j = Journey::new(60.0, start, DQuat::IDENTITY, 6371.0, 1.7e6);
        j.slow_zones(&[(40_000.0 * km_ly, 30.0), (8_000.0 * km_ly, 15.0)], end);
        let t_in = j.time_left(40_000.0 * km_ly, end).unwrap();
        let t_inside = j.time_left(8_000.0 * km_ly, end).unwrap();
        let left = |t: f64| (end - j.position_at(t, end)).length() / km_ly;
        assert!((left(t_in) - 40_000.0).abs() < 1.0 && (left(t_inside) - 8_000.0).abs() < 1.0);
        assert!((t_inside - t_in - 30.0).abs() < 1e-6 && (j.duration - t_inside - 15.0).abs() < 1e-6);
        let speed = |j: &Journey, t: f64| (j.position_at(t + 0.01, end) - j.position_at(t, end)).length() / 0.01 / km_ly;
        let p_in = plain.time_left(40_000.0 * km_ly, end).unwrap();
        let (coming_in, plain_in) = (speed(&j, t_in + 0.5), speed(&plain, p_in + 0.5 * (plain.time_left(8_000.0 * km_ly, end).unwrap() - p_in) / 30.0));
        // the inside: its speed just in, against just before the slowing gave way
        let (before, after) = (speed(&j, t_inside - 0.2), speed(&j, t_inside + 0.2));
        eprintln!("comes in at {t_in:.1} s at {coming_in:.0} ly/s (else {plain_in:.0}); {before:.0} ly/s 8000 ly out, then {after:.0}; arrives at {:.1} s", j.duration);
        assert!(coming_in < plain_in * 0.5, "slowed where it comes in");
        assert!(after > before, "quicker again inside");
        let mut step = j.advance(0.0, end, end);
        for _ in 0..(j.duration * 60.0) as usize + 60 {
            step = j.advance(1.0 / 60.0, end, end);
        }
        assert!(step.done && step.position == end);
    }

    /// Slowed smoothly into a galaxy: much slower where it comes in, quicker again inside, its
    /// speed never jumping, and it still arrives.
    #[test]
    fn slowed_smoothly_into_a_galaxy() {
        let km_ly = 9.4607e12;
        let (start, end) = (DVec3::ZERO, DVec3::new(2.5e6 * km_ly, 0.0, 0.0));
        let plain = Journey::new(60.0, start, DQuat::IDENTITY, 6371.0, 1.7e6);
        let mut j = Journey::new(60.0, start, DQuat::IDENTITY, 6371.0, 1.7e6);
        j.slow_by_left(&[(150_000.0 * km_ly, 1.0), (40_000.0 * km_ly, 0.05), (8_000.0 * km_ly, 0.05), (3_000.0 * km_ly, 0.3), (300.0 * km_ly, 1.0)], end);
        let speed = |j: &Journey, t: f64| (j.position_at(t + 0.01, end) - j.position_at(t, end)).length() / 0.01 / km_ly;
        let mut worst_jump: f64 = 1.0;
        let mut line = String::new();
        let mut t = 0.0;
        let mut last = speed(&j, 0.0);
        while t < j.duration - 0.05 {
            let v = speed(&j, t);
            if last > 1.0 && v > 1.0 {
                worst_jump = worst_jump.max((v / last).max(last / v));
            }
            last = v;
            t += 0.05;
        }
        for ly in [150_000.0, 40_000.0, 20_000.0, 8_000.0, 3_000.0, 1_000.0, 100.0] {
            let (tj, tp) = (j.time_left(ly * km_ly, end).unwrap(), plain.time_left(ly * km_ly, end).unwrap());
            line += &format!("{ly:.0} ly: {:.0} ly/s at {tj:.1} s (else {:.0} at {tp:.1} s); ", speed(&j, tj), speed(&plain, tp));
        }
        eprintln!("{line}arrives at {:.1} s; the speed changes by at most {worst_jump:.2} in 0.05 s", j.duration);
        let at = |ly: f64| speed(&j, j.time_left(ly * km_ly, end).unwrap());
        assert!(at(40_000.0) < 0.1 * speed(&plain, plain.time_left(40_000.0 * km_ly, end).unwrap()));
        assert!(at(3_000.0) > at(8_000.0) * 0.6, "quicker again (for how near the end it is)");
        assert!(worst_jump < 1.3, "no jump in speed: {worst_jump}");
        let mut step = j.advance(0.0, end, end);
        for _ in 0..(j.duration * 60.0) as usize + 60 {
            step = j.advance(1.0 / 60.0, end, end);
        }
        assert!(step.done && step.position == end);
    }
}
