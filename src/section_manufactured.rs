//! Manufactured transverse loading on an inherited multi-support displacement space.
use crate::{
    analysis_space::face_axes,
    bracket_json::decode_blocks,
    constrained_elasticity::solve_constrained,
    elasticity::gauss,
    extracted_space::assemble_sampled,
    incidence::enumerate_faces,
    inherited_space::InheritedSpace,
    json::{object, Json},
};
use std::collections::BTreeMap;
type V = [f64; 3];
type M = [V; 3];
fn n(v: impl ToString) -> Json {
    Json::Number(v.to_string())
}
fn dot(a: V, b: V) -> f64 {
    a.into_iter().zip(b).map(|(x, y)| x * y).sum()
}
fn cross(a: V, b: V) -> V {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
enum Profile {
    Constant,
    Power {
        degree: i32,
        origin: f64,
        width: f64,
    },
    Gaussian {
        center: V,
        lengths: V,
    },
}
impl Profile {
    fn jets(&self, x: V) -> (f64, V, M) {
        match self {
            Self::Constant => (1., [0.; 3], [[0.; 3]; 3]),
            Self::Power {
                degree,
                origin,
                width,
            } => {
                let t = (x[0] - origin) / width;
                let mut h = [[0.; 3]; 3];
                h[0][0] = f64::from(degree * (degree - 1)) * t.powi(degree - 2) / width.powi(2);
                (
                    t.powi(*degree),
                    [f64::from(*degree) * t.powi(degree - 1) / width, 0., 0.],
                    h,
                )
            }
            Self::Gaussian { center, lengths } => {
                let t: V = std::array::from_fn(|i| (x[i] - center[i]) / lengths[i]);
                let f = (-0.5 * dot(t, t)).exp();
                let a: V = std::array::from_fn(|i| -t[i] / lengths[i]);
                let d = a.map(|v| v * f);
                let h = std::array::from_fn(|i| {
                    std::array::from_fn(|j| {
                        f * (a[i] * a[j] - if i == j { 1. / lengths[i].powi(2) } else { 0. })
                    })
                });
                (f, d, h)
            }
        }
    }
}
struct Exact {
    profile: Profile,
    axial_sine: bool,
    z0: f64,
    height: f64,
    lambda: f64,
    mu: f64,
}
impl Exact {
    fn axial(&self, z: f64) -> (f64, f64, f64) {
        if self.axial_sine {
            let wave = std::f64::consts::PI / self.height;
            let phase = wave * (z - self.z0);
            (phase.sin(), wave * phase.cos(), -wave * wave * phase.sin())
        } else {
            ((z - self.z0) / self.height, 1. / self.height, 0.)
        }
    }
    fn field(&self, x: V) -> (V, M) {
        let (f, df, _) = self.profile.jets(x);
        let (h, dh, _) = self.axial(x[2]);
        let mut gradient = df.map(|a| a * h);
        gradient[2] += f * dh;
        ([0., 0., f * h], [[0.; 3], [0.; 3], gradient])
    }
    fn stress(&self, g: M) -> M {
        let div = g[0][0] + g[1][1] + g[2][2];
        std::array::from_fn(|i| {
            std::array::from_fn(|j| {
                self.mu * (g[i][j] + g[j][i]) + if i == j { self.lambda * div } else { 0. }
            })
        })
    }
    fn body(&self, x: V) -> V {
        let (f, df, d2f) = self.profile.jets(x);
        let (h, dh, d2h) = self.axial(x[2]);
        let laplace = h * (d2f[0][0] + d2f[1][1] + d2f[2][2]) + 2. * df[2] * dh + f * d2h;
        std::array::from_fn(|i| {
            let mixed = h * d2f[i][2] + dh * df[i] + if i == 2 { df[2] * dh + f * d2h } else { 0. };
            -(self.lambda + self.mu) * mixed - if i == 2 { self.mu * laplace } else { 0. }
        })
    }
}

/// E = 2.6 MPa, nu = 0.3. Prescribe zero displacement on the lowest exterior
/// face(s), use -div(sigma) in the volume and sigma*n over the full remaining
/// boundary, including pocket floors. The auxiliary bottom traction is also
/// integrated for the independent unrestricted virtual-work check.
pub fn solve(model: &Json, space: &InheritedSpace, job: &Json) -> Result<Json, &'static str> {
    let z0 = space
        .cells
        .iter()
        .map(|c| c.z_bounds[0])
        .fold(f64::INFINITY, f64::min);
    let z1 = space
        .cells
        .iter()
        .map(|c| c.z_bounds[1])
        .fold(f64::NEG_INFINITY, f64::max);
    let faces = enumerate_faces(&decode_blocks(model)?)?.exterior;
    let cells: BTreeMap<_, _> = space.cells.iter().map(|c| (c.id.as_str(), c)).collect();
    let bottom: Vec<_> = faces
        .iter()
        .filter(|s| s.axis == 2 && !s.upper && cells[s.block.as_str()].z_bounds[0] == z0)
        .cloned()
        .collect();
    let essential = space.zero_face_lift(&bottom)?;
    let indices: BTreeMap<_, _> = space
        .cells
        .iter()
        .enumerate()
        .map(|(i, c)| (c.id.as_str(), i))
        .collect();
    let boundary: Vec<_> = faces
        .iter()
        .map(|s| (indices[s.block.as_str()], s.axis, s.upper))
        .collect();
    let domain = IntegrationDomain {
        scalar_dofs: space.ndof,
        cell_count: space.cells.len(),
        z_bounds: [z0, z1],
        boundary: &boundary,
        essential: &essential,
    };
    let mut result = solve_sampled(&domain, job, |i, p| space.cells[i].sample(p))?;
    if let Json::Object(fields) = &mut result {
        for (key, value) in [
            ("inherited_joins", space.common_parent_faces),
            ("physical_c1_joins", space.transverse_faces),
            ("remaining_c0_joins", space.retained_faces),
        ] {
            fields.insert(key.into(), n(value));
        }
    }
    Ok(result)
}

pub struct IntegrationDomain<'a> {
    pub scalar_dofs: usize,
    pub cell_count: usize,
    pub z_bounds: [f64; 2],
    /// (leaf index, parametric axis, high side), excluding internal faces.
    pub boundary: &'a [(usize, usize, bool)],
    pub essential: &'a crate::constraint_system::Lift,
}

pub fn solve_sampled<F>(
    domain: &IntegrationDomain<'_>,
    job: &Json,
    sample: F,
) -> Result<Json, &'static str>
where
    F: Fn(
        usize,
        V,
    ) -> Result<
        (
            crate::extracted_space::ElementGeometrySample,
            crate::extracted_space::ElementBasisSample,
        ),
        &'static str,
    >,
{
    let [z0, z1] = domain.z_bounds;
    if !z0.is_finite()
        || !z1.is_finite()
        || z0 >= z1
        || domain.boundary.is_empty()
        || domain
            .boundary
            .iter()
            .any(|&(i, a, _)| i >= domain.cell_count || a >= 3)
    {
        return Err("invalid manufactured integration domain");
    }
    let schema = job.field("schema")?.text()?;
    let profile = match schema {
        "section_axial_sine_mms_v1" => {
            job.keys(&["schema", "quadrature"])?;
            Profile::Constant
        }
        "section_quadratic_mms_v1" | "section_cubic_mms_v1" | "section_quintic_mms_v1" => {
            job.keys(&["schema", "quadrature", "x_origin_mm", "x_scale_mm"])?;
            let origin = job.field("x_origin_mm")?.number::<f64>()?;
            let width = job.field("x_scale_mm")?.number::<f64>()?;
            if !origin.is_finite() || !width.is_finite() || width <= 0. {
                return Err("invalid manufactured polynomial scale");
            }
            Profile::Power {
                degree: match schema {
                    "section_cubic_mms_v1" => 3,
                    "section_quintic_mms_v1" => 5,
                    _ => 2,
                },
                origin,
                width,
            }
        }
        "local_gaussian_mms_v1" => {
            job.keys(&["schema", "quadrature", "center_mm", "lengths_mm"])?;
            let vector = |key: &str| -> Result<V, &'static str> {
                let v = job.field(key)?.items()?;
                if v.len() != 3 {
                    return Err("Gaussian load expects three components");
                }
                Ok([v[0].number()?, v[1].number()?, v[2].number()?])
            };
            let center = vector("center_mm")?;
            let lengths = vector("lengths_mm")?;
            if center.iter().any(|x| !x.is_finite())
                || lengths.iter().any(|x| !x.is_finite() || *x <= 0.)
            {
                return Err("invalid Gaussian load parameters");
            }
            Profile::Gaussian { center, lengths }
        }
        _ => return Err("unknown section verification job"),
    };
    let order = job.field("quadrature")?.number()?;
    let exact = Exact {
        profile,
        axial_sine: schema == "section_axial_sine_mms_v1",
        z0,
        height: z1 - z0,
        lambda: 1.5,
        mu: 1.,
    };
    let start = std::time::Instant::now();
    let mut assembly = assemble_sampled(
        domain.scalar_dofs,
        domain.cell_count,
        exact.lambda,
        exact.mu,
        order,
        &sample,
        |x| exact.body(x),
    )?;
    let q = gauss(order)?;
    for &(leaf, axis, upper) in domain.boundary {
        let [a, b] = face_axes(axis);
        for &(u, wu) in &q {
            for &(v, wv) in &q {
                let mut pt = [0.; 3];
                pt[a] = u;
                pt[b] = v;
                pt[axis] = usize::from(upper) as f64;
                let (g, basis) = sample(leaf, pt)?;
                let j = g.jacobian_columns;
                let area = cross(j[(axis + 1) % 3], j[(axis + 2) % 3]);
                let sigma = exact.stress(exact.field(g.position).1);
                let factor = wu * wv * if upper { 1. } else { -1. };
                for (&id, &r) in basis.dofs.iter().zip(&basis.values) {
                    for (i, &row) in sigma.iter().enumerate() {
                        assembly.force[3 * id + i] += r * dot(row, area) * factor;
                    }
                }
            }
        }
    }
    let assembly_ms = start.elapsed().as_secs_f64() * 1000.;
    let start = std::time::Instant::now();
    let solution = solve_constrained(
        &assembly.stiffness,
        [domain.essential; 3],
        &assembly.force,
        1e-11,
        12000,
    )?;
    let solve_ms = start.elapsed().as_secs_f64() * 1000.;
    let mut weak = vec![0.; 3 * domain.scalar_dofs];
    let mut integrated_uz = 0.;
    let mut totals = [0f64; 6];
    for leaf in 0..domain.cell_count {
        for &(u, wu) in &q {
            for &(v, wv) in &q {
                for &(w, ww) in &q {
                    let (g, b) = sample(leaf, [u, v, w])?;
                    let (x, d) = exact.field(g.position);
                    let mut uh = [0.; 3];
                    let mut du = [[0.; 3]; 3];
                    let sig = exact.stress(d);
                    let weight = g.determinant * wu * wv * ww;
                    for (k, &id) in b.dofs.iter().enumerate() {
                        for i in 0..3 {
                            let coeff = solution.displacement_mm[3 * id + i];
                            uh[i] += b.values[k] * coeff;
                            for (j, entry) in du[i].iter_mut().enumerate() {
                                *entry += b.physical_gradients[k][j] * coeff;
                            }
                            weak[3 * id + i] += weight * dot(sig[i], b.physical_gradients[k]);
                        }
                    }
                    integrated_uz += weight * uh[2];
                    let e: V = std::array::from_fn(|i| uh[i] - x[i]);
                    let de: M =
                        std::array::from_fn(|i| std::array::from_fn(|j| du[i][j] - d[i][j]));
                    let ds = exact.stress(de);
                    totals[0] += weight * dot(e, e);
                    totals[1] += weight * dot(x, x);
                    for k in 0..3 {
                        totals[2] += weight * dot(ds[k], de[k]);
                        totals[3] += weight * dot(sig[k], d[k]);
                        totals[4] += weight * dot(ds[k], ds[k]);
                        totals[5] += weight * dot(sig[k], sig[k]);
                    }
                }
            }
        }
    }
    let weak_defect = (weak
        .iter()
        .zip(&assembly.force)
        .map(|(a, b)| (a - b).powi(2))
        .sum::<f64>()
        / weak.iter().map(|x| x * x).sum::<f64>())
    .sqrt();
    if totals.iter().any(|v| !v.is_finite())
        || !weak_defect.is_finite()
        || !integrated_uz.is_finite()
        || [1, 3, 5].iter().any(|&i| totals[i] <= 0.)
    {
        return Err("indeterminate: manufactured physical norm is nonfinite or has zero reference");
    }
    let ratio = |i: usize| (totals[i] / totals[i + 1]).max(0.).sqrt();
    Ok(object([
        (
            "schema",
            Json::String("section_manufactured_result_v1".into()),
        ),
        ("manufactured_field", Json::String(schema.to_string())),
        ("quadrature", n(order)),
        ("scalar_dofs", n(domain.scalar_dofs)),
        ("free_vector_dofs", n(3 * domain.essential.free_roots.len())),
        (
            "mean_displacement_z_mm",
            n(integrated_uz / assembly.stiffness.integrated_volume),
        ),
        ("volume_mm3", n(assembly.stiffness.integrated_volume)),
        ("quadrature_points", n(assembly.quadrature_points)),
        ("block_nnz", n(assembly.stiffness.columns.len())),
        ("assembly_ms", n(assembly_ms)),
        ("solve_ms", n(solve_ms)),
        ("virtual_work_relative_defect", n(weak_defect)),
        ("relative_residual", n(solution.relative_residual)),
        ("energy_N_mm", n(solution.strain_energy_n_mm)),
        ("relative_displacement_L2", n(ratio(0))),
        ("relative_energy_norm", n(ratio(2))),
        ("relative_stress_L2", n(ratio(4))),
    ]))
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn gaussian_body_matches_independent_stress_divergence() {
        let exact = Exact {
            profile: Profile::Gaussian {
                center: [1.7, 2., 2.5],
                lengths: [0.4, 0.35, 0.55],
            },
            axial_sine: false,
            z0: 0.,
            height: 5.,
            lambda: 1.5,
            mu: 1.,
        };
        let x = [1.9, 2.1, 2.4];
        check_body(&exact, x);
        let cube = Exact {
            profile: Profile::Gaussian {
                center: [0.31, 0.57, 0.68],
                lengths: [0.30, 0.41, 0.36],
            },
            height: 1.,
            ..exact
        };
        check_body(&cube, [0.24, 0.49, 0.76]);
        let pocket = Exact {
            profile: Profile::Power {
                degree: 5,
                origin: -10.,
                width: 70.,
            },
            height: 20.,
            ..cube
        };
        check_body(&pocket, [17., 12., 7.]);
    }

    fn check_body(exact: &Exact, x: V) {
        let body = exact.body(x);
        let step = 1e-5;
        let mut numerical = [0.; 3];
        for j in 0..3 {
            let mut a = x;
            let mut b = x;
            a[j] -= step;
            b[j] += step;
            let stress_a = exact.stress(exact.field(a).1);
            let stress_b = exact.stress(exact.field(b).1);
            for i in 0..3 {
                numerical[i] -= (stress_b[i][j] - stress_a[i][j]) / (2. * step);
            }
        }
        for (expected, found) in body.into_iter().zip(numerical) {
            assert!((expected - found).abs() < 1e-7);
        }
    }

    #[test]
    fn axial_sine_body_matches_stress_divergence_and_bottom_trace_is_zero() {
        let exact = Exact {
            profile: Profile::Constant,
            axial_sine: true,
            z0: 0.,
            height: 1.,
            lambda: 1.5,
            mu: 1.,
        };
        let x = [0.37, 0.61, 0.42];
        check_body(&exact, x);
        assert_eq!(exact.field([0.5, 0.5, 0.]).0, [0.; 3]);
        assert!((exact.field([0.5, 0.5, 1.]).0[2]).abs() < 1e-15);
        assert!((exact.field([0.5, 0.5, 1.]).1[2][2] + std::f64::consts::PI).abs() < 1e-14);
    }
}
