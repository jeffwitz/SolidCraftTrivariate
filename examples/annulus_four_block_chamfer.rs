//! Edit the gallery annulus in its four original charts. The new axial knots
//! define three integration spans per volume. All nets are homogeneous in mm.
use solidcraft_trivariate::{bracket::BracketBlock, sections::section_blocks};

const HEIGHTS: [i64; 4] = [0, 1, 31, 32]; // quarter millimetres

fn control_net(source: &BracketBlock, chamfer: bool) -> Vec<[i64; 4]> {
    assert_eq!(source.degrees, [2, 1, 1]);
    let mut net = Vec::with_capacity(24);
    for angular in 0..3 {
        for radial in 0..2 {
            let [x, y, _, w] = source.controls[4 * angular + 2 * radial];
            for z in HEIGHTS {
                // Outer radius 4 -> 3.75 and inner radius 2 -> 2.25.
                let scale = if chamfer && (z == 0 || z == 32) {
                    [15, 18][radial]
                } else {
                    16
                };
                net.push([x * scale, y * scale, z * 4 * w, 16 * w]);
            }
        }
    }
    net
}

fn main() {
    let chamfer = match std::env::args().nth(1).as_deref() {
        Some("sharp") => false,
        Some("chamfer") | None => true,
        _ => panic!("expected sharp or chamfer"),
    };
    println!(
        r#"{{"degrees":[2,1,1],"knots":[[0,0,0,1,1,1],[0,0,1,1],[0,0,0.25,7.75,8,8]],"fastest_index":"axial","blocks":["#
    );
    let blocks = section_blocks("full_annulus_v1").unwrap();
    for (index, block) in blocks.iter().enumerate() {
        let comma = if index + 1 == blocks.len() { "" } else { "," };
        let net = control_net(block, chamfer);
        println!(r#"{{"id":"{}","xyzw":{net:?}}}{comma}"#, block.id);
    }
    println!("]}}");
}

#[cfg(test)]
mod tests {
    use super::*;
    use solidcraft_trivariate::spatial_certificate::{certify_spatial, CertificateStatus};

    #[test]
    fn three_closed_spans_per_original_volume_are_regular_and_match_their_neighbors() {
        let sources = section_blocks("full_annulus_v1").unwrap();
        assert_eq!(sources.len(), 4);
        for chamfer in [false, true] {
            let nets: Vec<_> = sources.iter().map(|s| control_net(s, chamfer)).collect();
            let mut lower = f64::INFINITY;
            for (i, net) in nets.iter().enumerate() {
                assert_eq!(net[16..24], nets[(i + 1) % 4][0..8]);
                for layer in 0..3 {
                    let span: Vec<_> = net
                        .chunks_exact(4)
                        .flat_map(|row| row[layer..=layer + 1].iter().copied())
                        .collect();
                    let result = certify_spatial([2, 1, 1], &span, Default::default());
                    assert_eq!(result.status, CertificateStatus::Valid);
                    lower = lower.min(result.determinant_lower_bound_mm3.unwrap());
                }
                if !chamfer {
                    for (old, new) in sources[i].controls.chunks_exact(2).zip(net.chunks_exact(4)) {
                        assert_eq!(new[0], old[0].map(|x| 16 * x));
                        assert_eq!(new[3], old[1].map(|x| 16 * x));
                    }
                }
            }
            println!("chamfer={chamfer}, closed span determinant >= {lower:.12} mm³");
        }
    }

    #[test]
    fn axial_stations_integrate_to_the_reference_step_volume() {
        let sources = section_blocks("full_annulus_v1").unwrap();
        for (chamfer, reference) in [
            (false, 96.0 * std::f64::consts::PI),
            (true, 299.2367002544223),
        ] {
            let net = control_net(&sources[0], chamfer);
            let areas = HEIGHTS.map(|h| {
                let k = HEIGHTS.iter().position(|&z| z == h).unwrap();
                let [outer, inner] =
                    [net[k], net[4 + k]].map(|[x, y, _, w]| (x as f64).hypot(y as f64) / w as f64);
                assert!((outer + inner - 6.).abs() < 1e-13);
                std::f64::consts::PI * (outer - inner) * (outer + inner)
            });
            // The sum of the radii is constant, so the area is affine on each span.
            let volume: f64 = (0..3)
                .map(|k| (areas[k] + areas[k + 1]) * (HEIGHTS[k + 1] - HEIGHTS[k]) as f64 / 8.)
                .sum();
            assert!((volume - reference).abs() < 1e-10);
        }
    }
}
