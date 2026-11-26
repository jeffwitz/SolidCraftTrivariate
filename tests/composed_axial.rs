#[path = "../examples/composed_axial.rs"]
mod composed_axial;

fn map_bounds(block: &solidcraft_trivariate::bracket::BracketBlock, axis: usize) -> [f64; 2] {
    let mut values = block
        .controls
        .iter()
        .map(|point| point[axis] as f64 / point[3] as f64);
    let first = values.next().unwrap();
    values.fold([first, first], |[lo, hi], value| {
        [lo.min(value), hi.max(value)]
    })
}

fn carried_coverage(result: &composed_axial::ComposedResult, point: [f64; 3]) -> usize {
    result
        .carried_maps
        .iter()
        .filter(|map| {
            (0..3).all(|axis| {
                let [lo, hi] = map_bounds(&map.block, axis);
                lo < point[axis] && point[axis] < hi
            })
        })
        .count()
}

fn carried_volume_mm3(result: &composed_axial::ComposedResult) -> f64 {
    result
        .carried_maps
        .iter()
        .map(|map| {
            (0..3)
                .map(|axis| {
                    let [lo, hi] = map_bounds(&map.block, axis);
                    hi - lo
                })
                .product::<f64>()
        })
        .sum()
}

#[derive(Debug)]
struct TensionComponent {
    range: [i32; 2],
    tip_displacement_mm: f64,
    strain_energy_n_mm: f64,
}

fn solve_dense(mut matrix: Vec<Vec<f64>>, mut rhs: Vec<f64>) -> Vec<f64> {
    let n = rhs.len();
    for pivot in 0..n {
        let best = (pivot..n)
            .max_by(|&a, &b| matrix[a][pivot].abs().total_cmp(&matrix[b][pivot].abs()))
            .unwrap();
        assert!(matrix[best][pivot].abs() > 1e-14);
        matrix.swap(pivot, best);
        rhs.swap(pivot, best);
        let diagonal = matrix[pivot][pivot];
        for entry in matrix[pivot].iter_mut().skip(pivot) {
            *entry /= diagonal;
        }
        rhs[pivot] /= diagonal;
        let pivot_row = matrix[pivot].clone();
        for (row_index, row) in matrix.iter_mut().enumerate() {
            if row_index == pivot {
                continue;
            }
            let factor = row[pivot];
            for (j, entry) in row.iter_mut().enumerate().skip(pivot) {
                *entry -= factor * pivot_row[j];
            }
            rhs[row_index] -= factor * rhs[pivot];
        }
    }
    rhs
}

fn solve_carried_tension(
    result: &composed_axial::ComposedResult,
    young_modulus_n_mm2: f64,
    area_mm2: f64,
    force_n: f64,
) -> Vec<TensionComponent> {
    use std::collections::BTreeSet;

    let mut elements: Vec<_> = result
        .carried_maps
        .iter()
        .map(|map| map.result_interval_mm)
        .collect();
    elements.sort_unstable();
    assert!(!elements.is_empty());
    assert!(elements.windows(2).all(|w| w[0] != w[1]));

    let mut ranges = Vec::<[i32; 2]>::new();
    for [lo, hi] in &elements {
        assert!(lo < hi);
        if let Some(last) = ranges.last_mut() {
            if *lo < last[1] {
                panic!("carried axial elements overlap in their interiors");
            }
            if *lo == last[1] {
                last[1] = *hi;
                continue;
            }
        }
        ranges.push([*lo, *hi]);
    }

    ranges
        .into_iter()
        .map(|range| {
            let local_elements: Vec<_> = elements
                .iter()
                .copied()
                .filter(|e| range[0] <= e[0] && e[1] <= range[1])
                .collect();
            let nodes: Vec<_> = local_elements
                .iter()
                .flat_map(|e| [e[0], e[1]])
                .collect::<BTreeSet<_>>()
                .into_iter()
                .collect();
            assert_eq!(nodes.first().copied(), Some(range[0]));
            assert_eq!(nodes.last().copied(), Some(range[1]));
            let index: std::collections::BTreeMap<_, _> =
                nodes.iter().enumerate().map(|(i, z)| (*z, i)).collect();
            let mut stiffness = vec![vec![0.0; nodes.len()]; nodes.len()];
            for [lo, hi] in &local_elements {
                let length = f64::from(hi - lo);
                let k = young_modulus_n_mm2 * area_mm2 / length;
                let i = index[lo];
                let j = index[hi];
                stiffness[i][i] += k;
                stiffness[i][j] -= k;
                stiffness[j][i] -= k;
                stiffness[j][j] += k;
            }

            // Dirichlet clamp at the component's lower endpoint, traction at its upper endpoint.
            let free_count = nodes.len() - 1;
            let reduced = stiffness[1..].iter().map(|row| row[1..].to_vec()).collect();
            let mut rhs = vec![0.0; free_count];
            *rhs.last_mut().unwrap() = force_n;
            let free = solve_dense(reduced, rhs);
            let mut displacement = vec![0.0; nodes.len()];
            displacement[1..].copy_from_slice(&free);
            let tip_displacement_mm = displacement[nodes.len() - 1];
            let strain_energy_n_mm = 0.5
                * displacement
                    .iter()
                    .enumerate()
                    .map(|(i, ui)| {
                        displacement
                            .iter()
                            .enumerate()
                            .map(|(j, uj)| ui * stiffness[i][j] * uj)
                            .sum::<f64>()
                    })
                    .sum::<f64>();
            assert!((strain_energy_n_mm - 0.5 * force_n * tip_displacement_mm).abs() < 1e-13);
            let analytic_tip =
                force_n * f64::from(range[1] - range[0]) / (young_modulus_n_mm2 * area_mm2);
            let analytic_energy = force_n * force_n * f64::from(range[1] - range[0])
                / (2.0 * young_modulus_n_mm2 * area_mm2);
            assert!((tip_displacement_mm - analytic_tip).abs() <= 1e-11 * analytic_tip.abs());
            assert!((strain_energy_n_mm - analytic_energy).abs() <= 1e-11 * analytic_energy.abs());
            for (z, u) in nodes.iter().zip(displacement) {
                let exact = force_n * f64::from(*z - range[0]) / (young_modulus_n_mm2 * area_mm2);
                assert!((u - exact).abs() <= 1e-11 * analytic_tip.abs());
            }
            TensionComponent {
                range,
                tip_displacement_mm,
                strain_energy_n_mm,
            }
        })
        .collect()
}

#[test]
fn composed_union_difference_preserves_volume_containment_and_lineage() {
    let cases = [
        ([0, 8], [4, 12], [5, 9]),
        ([-2, 5], [2, 10], [-1, 1]),
        ([0, 6], [6, 12], [4, 8]),
        ([-8, 2], [-4, 6], [1, 4]),
    ];
    for (a, b, c) in cases {
        let result = composed_axial::compose(a, b, c).unwrap();
        let first_length = a[1].max(b[1]) - a[0].min(b[0]);
        assert!((result.first.volume_mm3 - 400.0 * f64::from(first_length)).abs() < 1e-10);
        let union_lo = a[0].min(b[0]);
        let union_hi = a[1].max(b[1]);
        assert_eq!(result.composite_operand.interval_mm, [union_lo, union_hi]);

        let overlap = (union_hi.min(c[1]) - union_lo.max(c[0])).max(0);
        let expected_volume = 400.0 * f64::from((union_hi - union_lo) - overlap);
        assert!((result.second.volume_mm3 - expected_volume).abs() < 1e-10);
        let carried_volume = carried_volume_mm3(&result);
        assert!((carried_volume - expected_volume).abs() < 1e-10);
        let tension = solve_carried_tension(&result, 210_000.0, 400.0, 1_000.0);
        assert_eq!(tension.len(), 2);
        assert!(tension.iter().all(|component| {
            component.range[0] < component.range[1]
                && component.tip_displacement_mm > 0.0
                && component.strain_energy_n_mm > 0.0
        }));

        // Independent point classification over interior section points and axial cells.
        for z4 in (4 * (union_lo - 1) + 1..4 * (union_hi + 1)).step_by(2) {
            let z = f64::from(z4) / 4.0;
            let expected = (f64::from(a[0]) < z && z < f64::from(a[1])
                || f64::from(b[0]) < z && z < f64::from(b[1]))
                && !(f64::from(c[0]) < z && z < f64::from(c[1]));
            for x in [-9.0, 0.0, 9.0] {
                for y in [-9.0, 0.0, 9.0] {
                    let coverage = carried_coverage(&result, [x, y, z]);
                    assert_eq!(
                        coverage,
                        usize::from(expected),
                        "carried map coverage at ({x},{y},{z}) for {a:?}, {b:?}, {c:?}"
                    );
                }
            }
        }
        let lineage: Vec<_> = result
            .lineage
            .iter()
            .flat_map(|(_, labels)| labels)
            .copied()
            .collect();
        assert!(lineage.contains(&"C"));
        assert!(lineage.iter().any(|&label| label == "A" || label == "B"));
        if a == [0, 8] && b == [4, 12] && c == [5, 9] {
            assert!(lineage.contains(&"A"));
            assert!(lineage.contains(&"B"));
        }
        assert!(result.carried_maps.iter().any(|map| {
            map.source_interval_mm != map.result_interval_mm
                && map.source_block_id.starts_with("axial_boolean_cell_")
                && !map.original_operands.is_empty()
        }));
    }
}

#[test]
fn disconnected_intermediate_is_not_silently_reduced_to_its_bounding_interval() {
    // This union is connected only if its operands touch or overlap; it cannot be
    // represented by the single-interval AxialOperand consumed by stage two.
    assert!(composed_axial::compose([0, 2], [5, 8], [1, 6]).is_err());
}
