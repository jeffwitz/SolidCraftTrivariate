//! Four-sector bore hosts with exact rectangular infill.
//!
//! Parallel bores first use full-section strips separated at exact integer
//! mid-gaps. A single bore can use its largest centered square host plus box
//! infill, which keeps distant stock ends outside the curved charts. Unsupported
//! overlaps are left to the established planner fallback.

use crate::{
    bracket::BracketBlock,
    orthogonal_sketch::{self, Bore, OrthogonalSketch, SweepPartition},
    sketch_layers::{area2, edges},
};
use std::collections::{BTreeMap, BTreeSet};

pub(crate) type Bounds = [[i32; 2]; 3];

pub(crate) struct HoledBox {
    pub blocks: Vec<BracketBlock>,
    pub physical: BTreeMap<(String, usize, bool), String>,
    pub guards: Vec<(String, Bounds)>,
}

fn bbox(j: &OrthogonalSketch) -> Bounds {
    std::array::from_fn(|k| {
        if k == 2 {
            [0, j.height]
        } else {
            j.outer.vertices.iter().fold([i32::MAX, i32::MIN], |r, p| {
                [r[0].min(p[k]), r[1].max(p[k])]
            })
        }
    })
}

fn host(bounds: Bounds, bore: &Bore) -> Option<Bounds> {
    if bore.axis > 2 || bore.radius <= 0 {
        return None;
    }
    if bore.center[bore.axis] <= bounds[bore.axis][0]
        || bore.center[bore.axis] >= bounds[bore.axis][1]
    {
        return None;
    }
    let a = (bore.axis + 1) % 3;
    let b = (bore.axis + 2) % 3;
    // The largest centered square contained in the body gives a scale-aware
    // local host. Its circle remains strictly interior by the input contract.
    let half = [a, b]
        .into_iter()
        .map(|k| (bore.center[k] - bounds[k][0]).min(bounds[k][1] - bore.center[k]))
        .min()?;
    if half <= bore.radius {
        return None;
    }
    let mut h = bounds;
    for k in [a, b] {
        h[k] = [bore.center[k] - half, bore.center[k] + half];
    }
    Some(h)
}

fn hosts_for_group<'a>(
    bounds: Bounds,
    axis: usize,
    bores: &[&'a Bore],
    full_single: bool,
) -> Option<Vec<(&'a Bore, Bounds)>> {
    if bores.len() == 1 {
        return Some(vec![(
            bores[0],
            if full_single {
                bounds
            } else {
                host(bounds, bores[0])?
            },
        )]);
    }
    if bores.len() > 1 {
        for split_axis in [(axis + 1) % 3, (axis + 2) % 3] {
            let mut ordered = bores.to_vec();
            ordered.sort_by_key(|b| b.center[split_axis]);
            if ordered.windows(2).any(|w| {
                let mid = w[0].center[split_axis]
                    + (w[1].center[split_axis] - w[0].center[split_axis]) / 2;
                w[0].center[split_axis] + w[0].radius >= mid
                    || w[1].center[split_axis] - w[1].radius <= mid
            }) {
                continue;
            }
            let other = if split_axis == (axis + 1) % 3 {
                (axis + 2) % 3
            } else {
                (axis + 1) % 3
            };
            if ordered.iter().any(|b| {
                b.center[other] - b.radius <= bounds[other][0]
                    || b.center[other] + b.radius >= bounds[other][1]
            }) {
                continue;
            }
            let mut cuts = Vec::with_capacity(ordered.len() + 1);
            cuts.push(bounds[split_axis][0]);
            for pair in ordered.windows(2) {
                cuts.push(
                    pair[0].center[split_axis]
                        + (pair[1].center[split_axis] - pair[0].center[split_axis]) / 2,
                );
            }
            cuts.push(bounds[split_axis][1]);
            let strips = ordered
                .iter()
                .enumerate()
                .map(|(i, bore)| {
                    let mut h = bounds;
                    h[split_axis] = [cuts[i], cuts[i + 1]];
                    (*bore, h)
                })
                .collect::<Vec<_>>();
            if strips.iter().all(|(b, h)| {
                h[split_axis][0] < b.center[split_axis] - b.radius
                    && b.center[split_axis] + b.radius < h[split_axis][1]
            }) {
                return Some(strips);
            }
        }
    }
    let initial = bores
        .iter()
        .map(|b| host(bounds, b).map(|h| (*b, h)))
        .collect::<Option<Vec<_>>>()?;
    let overlap = |hs: &[(&Bore, Bounds)]| {
        (0..hs.len()).any(|i| {
            (0..i).any(|j| {
                (0..3).all(|k| hs[i].1[k][0].max(hs[j].1[k][0]) < hs[i].1[k][1].min(hs[j].1[k][1]))
            })
        })
    };
    if !overlap(&initial) {
        return Some(initial);
    }
    None
}

pub(crate) fn plain_box(id: String, c: Bounds) -> BracketBlock {
    let mut controls = Vec::with_capacity(8);
    for x in c[0] {
        for y in c[1] {
            for z in c[2] {
                controls.push([i64::from(x), i64::from(y), i64::from(z), 1]);
            }
        }
    }
    BracketBlock {
        id,
        degrees: [1, 1, 1],
        controls,
    }
}

/// Construct one or more parallel, disjoint bore hosts in a rectangular prism.
/// Returns `None` when the host squares overlap or a bore is not strictly
/// interior; certification failures are errors so the caller can retain its
/// established fallback.
pub(crate) fn holed_box(
    axis: usize,
    bounds: Bounds,
    bores: &[&Bore],
    prefix: &str,
    full_single: bool,
) -> Result<Option<HoledBox>, &'static str> {
    if bores.is_empty() || axis > 2 || bores.iter().any(|b| b.axis != axis) {
        return Ok(None);
    }
    if bores.iter().any(|b| {
        (0..3).any(|k| {
            if k == axis {
                b.center[k] <= bounds[k][0] || b.center[k] >= bounds[k][1]
            } else {
                b.center[k] - b.radius <= bounds[k][0] || b.center[k] + b.radius >= bounds[k][1]
            }
        })
    }) {
        return Ok(None);
    }
    let Some(hosts) = hosts_for_group(bounds, axis, bores, full_single) else {
        return Ok(None);
    };
    let mut levels: [BTreeSet<i32>; 3] = std::array::from_fn(|k| BTreeSet::from(bounds[k]));
    for (_, h) in &hosts {
        for k in 0..3 {
            levels[k].extend(h[k]);
        }
    }
    let levels = levels.map(|v| v.into_iter().collect::<Vec<_>>());
    let mut blocks = Vec::new();
    let mut physical = BTreeMap::new();
    for (b, h) in &hosts {
        let charts = crate::compact_wing::build(
            axis,
            *h,
            b.center,
            b.radius,
            &format!("{prefix}_{}", b.id),
        )?;
        if charts.len() != 4 {
            return Err("compact_host_sector_count");
        }
        for chart in charts {
            orthogonal_sketch::certify_circular_side(&chart, axis, b.center, b.radius)?;
            physical.insert((chart.id.clone(), 1, true), format!("bore/{}/wall", b.id));
            blocks.push(chart);
        }
    }
    for ix in levels[0].windows(2) {
        for iy in levels[1].windows(2) {
            for iz in levels[2].windows(2) {
                let c: Bounds = [[ix[0], ix[1]], [iy[0], iy[1]], [iz[0], iz[1]]];
                if hosts
                    .iter()
                    .any(|(_, h)| (0..3).all(|k| h[k][0] <= c[k][0] && c[k][1] <= h[k][1]))
                {
                    continue;
                }
                let id = format!(
                    "{prefix}_fill_{}_{}_{}_{}_{}_{}",
                    c[0][0], c[0][1], c[1][0], c[1][1], c[2][0], c[2][1]
                );
                blocks.push(plain_box(id, c));
            }
        }
    }
    let guards = hosts.into_iter().map(|(b, h)| (b.id.clone(), h)).collect();
    Ok(Some(HoledBox {
        blocks,
        physical,
        guards,
    }))
}

/// Compact local-host path for a rectangular plate with one or more parallel bores.
pub(crate) fn try_plate(j: &OrthogonalSketch) -> Result<Option<SweepPartition>, &'static str> {
    if !j.fillets.is_empty()
        || j.outer.vertices.len() != 4
        || j.bores.is_empty()
        || j.height <= 0
        || edges(&j.outer).any(|(a, b)| a[0] != b[0] && a[1] != b[1])
    {
        return Ok(None);
    }
    let bounds = bbox(j);
    let axis = j.bores[0].axis;
    let Some(HoledBox {
        blocks,
        physical,
        guards,
    }) = holed_box(
        axis,
        bounds,
        &j.bores.iter().collect::<Vec<_>>(),
        "plate",
        false,
    )?
    else {
        return Ok(None);
    };
    let mut planes: [BTreeSet<i32>; 3] = std::array::from_fn(|k| BTreeSet::from(bounds[k]));
    for (_, h) in &guards {
        for k in 0..3 {
            planes[k].extend(h[k]);
        }
    }
    let planes = planes.map(|p| p.into_iter().collect());
    let four_pi = j.bores.iter().try_fold(0_i128, |s, b| {
        s.checked_sub(
            4 * i128::from(b.radius).pow(2) * i128::from(bounds[b.axis][1] - bounds[b.axis][0]),
        )
        .ok_or("compact_host_volume_overflow")
    })?;
    let p = orthogonal_sketch::finish_partition(
        j,
        blocks,
        physical,
        i128::from(area2(&j.outer)) * i128::from(j.height),
        four_pi,
        planes,
        guards,
    )?;
    Ok(Some(p))
}
