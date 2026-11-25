//! Ordered regularized circular cuts and additions in an axially layered plate.
use crate::bracket::BracketBlock;
use std::collections::{BTreeSet, VecDeque};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CylinderOp {
    Add,
    Cut,
}
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CylinderFeature {
    pub id: String,
    pub op: CylinderOp,
    pub radius_mm: i32,
    pub z_mm: [i32; 2],
}
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CoaxialProgram {
    pub base_height_mm: i32,
    pub features: Vec<CylinderFeature>,
}

/// The most recent feature that actually changed the occupancy of this region.
/// Zero denotes the initial plate. Positive indices start at one in the history.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Region {
    pub material: bool,
    pub changed_by: usize,
}

#[derive(Clone, Debug)]
pub struct LayeredRegions {
    pub cuts_z: Vec<i32>,
    pub radii: Vec<i32>,
    /// Each row is an axial slab; column zero contains the axis and the final
    /// column reaches the outer square. Each boundary circle is shared.
    pub slabs: Vec<Vec<Region>>,
}

impl CoaxialProgram {
    pub fn regions(&self) -> Result<LayeredRegions, &'static str> {
        if !(1..=16).contains(&self.base_height_mm) || self.features.len() > 8 {
            return Err("coaxial_program_range");
        }
        let mut names = BTreeSet::new();
        let mut heights = BTreeSet::from([0, self.base_height_mm]);
        let mut radii = BTreeSet::new();
        for feature in &self.features {
            let id = &feature.id;
            if id.is_empty()
                || id.len() > 48
                || id == "base"
                || !id
                    .bytes()
                    .all(|c| c.is_ascii_alphanumeric() || c == b'_' || c == b'-')
                || !names.insert(id)
            {
                return Err("coaxial_feature_id");
            }
            let [lo, hi] = feature.z_mm;
            if !(1..=8).contains(&feature.radius_mm) || lo < 0 || hi > 24 || lo >= hi {
                return Err("coaxial_feature_range");
            }
            heights.extend([lo, hi]);
            radii.insert(feature.radius_mm);
        }
        let cuts_z: Vec<_> = heights.into_iter().collect();
        let radii: Vec<_> = radii.into_iter().collect();
        let slabs = cuts_z
            .windows(2)
            .map(|z| {
                (0..=radii.len())
                    .map(|band| {
                        let mut state = Region {
                            material: z[1] <= self.base_height_mm,
                            changed_by: 0,
                        };
                        for (index, feature) in self.features.iter().enumerate() {
                            if z[0] >= feature.z_mm[0]
                                && z[1] <= feature.z_mm[1]
                                && radii.get(band).is_some_and(|&r| r <= feature.radius_mm)
                            {
                                let next = feature.op == CylinderOp::Add;
                                if state.material != next {
                                    state = Region {
                                        material: next,
                                        changed_by: index + 1,
                                    };
                                }
                            }
                        }
                        state
                    })
                    .collect()
            })
            .collect();
        let grid = LayeredRegions {
            cuts_z,
            radii,
            slabs,
        };
        grid.check_material()?;
        Ok(grid)
    }
}

impl LayeredRegions {
    fn check_material(&self) -> Result<(), &'static str> {
        let nz = self.slabs.len();
        let nr = self.radii.len() + 1;
        for pair in self.slabs.windows(2) {
            for r in 0..nr - 1 {
                let a = pair[0][r].material;
                let b = pair[0][r + 1].material;
                if a != b && a == pair[1][r + 1].material && b == pair[1][r].material {
                    return Err("coaxial_nonmanifold_edge_contact");
                }
            }
        }
        let total = self.slabs.iter().flatten().filter(|r| r.material).count();
        let (z, r) = self
            .slabs
            .iter()
            .enumerate()
            .find_map(|(z, row)| row.iter().position(|s| s.material).map(|r| (z, r)))
            .ok_or("coaxial_empty_material")?;
        let mut seen = vec![vec![false; nr]; nz];
        seen[z][r] = true;
        let mut queue = VecDeque::from([(z, r)]);
        let mut count = 0;
        while let Some((z, r)) = queue.pop_front() {
            count += 1;
            for (a, b) in [
                (z.wrapping_sub(1), r),
                (z + 1, r),
                (z, r.wrapping_sub(1)),
                (z, r + 1),
            ] {
                if a < nz && b < nr && self.slabs[a][b].material && !seen[a][b] {
                    seen[a][b] = true;
                    queue.push_back((a, b));
                }
            }
        }
        if count != total {
            return Err("coaxial_disconnected_material");
        }
        Ok(())
    }
}

type Curve = [[i64; 3]; 3];

/// Angular charts alternate their quarter-circle gauge at cardinal endpoints.
/// This makes the radial seam's weights on both sheets identical, including
/// the transition between constant-weight straight edges and round profiles.
fn circle_curve(radius: i32, sector: usize) -> Curve {
    let q = sector / 2;
    let mut half = sector % 2;
    if q % 2 == 1 {
        half = 1 - half;
    }
    let mut net = match half {
        0 => [[4, 0, 4], [4, 2, 4], [3, 4, 5]],
        _ => [[3, 4, 5], [2, 6, 6], [0, 8, 8]],
    };
    if q % 2 == 1 {
        net.reverse();
    }
    for p in &mut net {
        let [x, y, w] = *p;
        let (x, y) = match q {
            0 => (x, y),
            1 => (-x, y),
            2 => (-x, -y),
            _ => (x, -y),
        };
        *p = [4 * i64::from(radius) * x, 4 * i64::from(radius) * y, 4 * w];
    }
    net
}

fn square_curve(quarter_mm: i64, sector: usize) -> Curve {
    let points = [
        [1, 0],
        [1, 1],
        [0, 1],
        [-1, 1],
        [-1, 0],
        [-1, -1],
        [0, -1],
        [1, -1],
    ];
    let [x, y] = points[sector];
    let [xx, yy] = points[(sector + 1) % 8];
    [
        [2 * quarter_mm * x, 2 * quarter_mm * y, 8],
        [quarter_mm * (x + xx), quarter_mm * (y + yy), 8],
        [2 * quarter_mm * xx, 2 * quarter_mm * yy, 8],
    ]
}
fn rectangle(id: String, x: [i64; 2], y: [i64; 2], z: [i32; 2]) -> BracketBlock {
    let mut controls = Vec::new();
    for x in x {
        for y in y {
            for z in z {
                controls.push([x, y, 4 * i64::from(z), 4]);
            }
        }
    }
    BracketBlock {
        id,
        degrees: [1, 1, 1],
        controls,
    }
}
fn sector_block(id: String, outer: Curve, inner: Curve, z: [i32; 2]) -> BracketBlock {
    let mut controls = Vec::new();
    for i in 0..3 {
        for [x, y, w] in [outer[i], inner[i]] {
            for z in z {
                controls.push([x, y, i64::from(z) * w, w]);
            }
        }
    }
    BracketBlock {
        id,
        degrees: [2, 1, 1],
        controls,
    }
}

/// Build the exact integer solids of an admissible ordered coaxial history.
pub fn build_coaxial_program(program: &CoaxialProgram) -> Result<Vec<BracketBlock>, &'static str> {
    let grid = program.regions()?;
    let mut blocks = Vec::new();
    for (layer, z) in grid.cuts_z.windows(2).enumerate() {
        let z = [z[0], z[1]];
        if grid.radii.is_empty() {
            blocks.push(rectangle(
                format!("plate_{}_{}", z[0], z[1]),
                [-40, 40],
                [-40, 40],
                z,
            ));
            continue;
        }
        let core = i64::from(grid.radii[0]);
        for (band, region) in grid.slabs[layer].iter().enumerate() {
            if !region.material {
                continue;
            }
            if band == 0 {
                for (i, x) in [[-core, 0], [0, core]].into_iter().enumerate() {
                    for (j, y) in [[-core, 0], [0, core]].into_iter().enumerate() {
                        blocks.push(rectangle(
                            format!("core_{}_{}_{}_{}", z[0], z[1], i, j),
                            x,
                            y,
                            z,
                        ));
                    }
                }
            }
            for sector in 0..8 {
                let inner = if band == 0 {
                    square_curve(core, sector)
                } else {
                    circle_curve(grid.radii[band - 1], sector)
                };
                let outer = if band == grid.radii.len() {
                    square_curve(40, sector)
                } else {
                    circle_curve(grid.radii[band], sector)
                };
                blocks.push(sector_block(
                    format!("band_{}_{}_{}_{}", z[0], z[1], band, sector),
                    outer,
                    inner,
                    z,
                ));
            }
        }
    }
    Ok(blocks)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::incidence::enumerate_faces;
    use crate::planar::certify_planar_refined;

    fn feature(id: &str, op: CylinderOp, radius_mm: i32, z_mm: [i32; 2]) -> CylinderFeature {
        CylinderFeature {
            id: id.into(),
            op,
            radius_mm,
            z_mm,
        }
    }

    #[test]
    fn sleeve_and_pierced_boss_have_distinct_ordered_histories() {
        let boss = feature("boss", CylinderOp::Add, 6, [8, 12]);
        let bore = feature("bore", CylinderOp::Cut, 4, [0, 12]);
        let a = CoaxialProgram {
            base_height_mm: 8,
            features: vec![boss.clone(), bore.clone()],
        };
        let b = CoaxialProgram {
            base_height_mm: 8,
            features: vec![bore, boss],
        };
        let first = a.regions().unwrap();
        let second = b.regions().unwrap();
        assert_eq!(
            first.slabs[1][0],
            Region {
                material: false,
                changed_by: 2
            }
        );
        assert_eq!(
            second.slabs[1][0],
            Region {
                material: true,
                changed_by: 2
            }
        );
        let blocks = build_coaxial_program(&a).unwrap();
        let topology = enumerate_faces(&blocks).unwrap();
        assert_eq!(
            2 * topology.interfaces.len() + topology.exterior.len(),
            6 * blocks.len()
        );
    }

    #[test]
    fn all_integer_core_and_outer_rings_have_continuously_positive_maps() {
        for radius in 1..=8 {
            let program = CoaxialProgram {
                base_height_mm: 8,
                features: vec![feature("pocket", CylinderOp::Cut, radius, [4, 8])],
            };
            let blocks = build_coaxial_program(&program).unwrap();
            let topo = enumerate_faces(&blocks).unwrap();
            assert_eq!(
                2 * topo.interfaces.len() + topo.exterior.len(),
                6 * blocks.len()
            );
            for b in blocks {
                let net = b
                    .controls
                    .iter()
                    .step_by(2)
                    .map(|p| [p[0], p[1], p[3]])
                    .collect::<Vec<_>>();
                assert!(
                    certify_planar_refined([b.degrees[0], b.degrees[1]], &net).unwrap() > 0.
                );
            }
        }
    }

    #[test]
    fn rejects_a_contact_ring_and_a_detached_plug() {
        let pocket = feature("pocket", CylinderOp::Cut, 5, [4, 8]);
        for (r, lo, reason) in [
            (5, 8, "coaxial_nonmanifold_edge_contact"),
            (3, 6, "coaxial_disconnected_material"),
        ] {
            let p = CoaxialProgram {
                base_height_mm: 8,
                features: vec![
                    pocket.clone(),
                    feature("boss", CylinderOp::Add, r, [lo, 10]),
                ],
            };
            assert_eq!(p.regions().unwrap_err(), reason);
        }
    }
}
