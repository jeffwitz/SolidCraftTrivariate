//! Three edge states on the four charts of full_annulus_v1.
//! The fillet loses rank at the four circular rails on the end faces.
use crate::{bracket::BracketBlock, sections::section_blocks};

use super::EdgeProfile as Profile;

fn net(source: &BracketBlock, state: Profile) -> Vec<[i64; 4]> {
    assert_eq!(source.degrees, [2, 1, 1]);
    let mut controls = Vec::with_capacity(63);
    for angle in 0..3 {
        let [x, y, _, w] = source.controls[4 * angle]; // outer angular row, nominal radius 4
        for radial in 0..3_i64 {
            for layer in 0..7 {
                let stage = layer.min(6 - layer);
                let weight = [1, 1, 2, 2][stage];
                // Homogeneous setback in eighths of a millimetre.
                let setback = match state {
                    Profile::Sharp => 0,
                    Profile::Chamfer => [2, 1, 0, 0][stage],
                    Profile::Fillet => [2, 0, 0, 0][stage],
                };
                let radius = 8 * (4 - radial) * weight + (radial - 1) * setback;
                let mut height = if state == Profile::Fillet {
                    // The interior radial pole lifts the cap parametrization away
                    // from the edge; its coefficient vanishes on both walls.
                    [0, i64::from(radial == 1), 4, 64][stage]
                } else {
                    [0, 1, 4, 64][stage]
                };
                if layer > 3 {
                    height = 64 * weight - height;
                }
                controls.push([x * radius, y * radius, 4 * w * height, 32 * w * weight]);
            }
        }
    }
    controls
}

use super::{numbers, Json};
pub(super) fn document(state: Profile) -> Result<Json, &'static str> {
    let blocks = section_blocks("full_annulus_v1")?;
    Ok(super::original_document(
        [2, 2, 2],
        [
            numbers([0, 0, 0, 1, 1, 1]),
            numbers([0, 0, 0, 1, 1, 1]),
            super::axial_knots(2),
        ],
        blocks.iter().map(|b| (&b.id, net(b, state))),
    ))
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn all_three_states_preserve_the_same_weights_and_all_four_exact_interfaces() {
        let sectors = section_blocks("full_annulus_v1").unwrap();
        assert_eq!(sectors.len(), 4);
        for state in [Profile::Sharp, Profile::Chamfer, Profile::Fillet] {
            let nets: Vec<_> = sectors.iter().map(|s| net(s, state)).collect();
            for (i, poles) in nets.iter().enumerate() {
                assert_eq!(poles.len(), 63);
                assert_eq!(poles[42..], nets[(i + 1) % 4][..21]);
                assert!(poles
                    .iter()
                    .zip(net(&sectors[i], Profile::Sharp))
                    .all(|(a, b)| a[3] == b[3]));
            }
        }
    }
}
