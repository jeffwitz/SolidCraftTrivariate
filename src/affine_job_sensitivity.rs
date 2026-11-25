//! Fixed-regime material derivative of a displacement job under a global affine velocity.
use crate::{
    analysis_space::AnalysisSpace,
    constraint_system::Lift,
    elasticity::Vector,
    elasticity_sensitivity::{assemble_sensitivity, MaterialJet},
    json::{object, Json},
    mechanical_sensitivity::{differentiate, LiftRate, MechanicalDerivative},
};

pub(crate) struct Velocity {
    linear: [Vector; 3],
    offset: Vector,
    endpoint: Option<EndpointFlow>,
    pub(crate) exact: Option<crate::prescribed_motion::Motion>,
}

/// Piecewise-affine vertical motion of one physical cap and its sheared axis.
struct EndpointFlow {
    lower: f64,
    moving: f64,
    upper: f64,
    axis: Vector,
}

impl EndpointFlow {
    fn amplitude(&self, z: f64) -> f64 {
        if z < self.moving {
            ((z - self.lower) / (self.moving - self.lower)).clamp(0., 1.)
        } else {
            ((self.upper - z) / (self.upper - self.moving)).clamp(0., 1.)
        }
    }

    fn slope(&self, z: f64) -> f64 {
        if z > self.lower && z < self.moving {
            1. / (self.moving - self.lower)
        } else if z >= self.moving && z < self.upper {
            -1. / (self.upper - self.moving)
        } else {
            0.
        }
    }
}

impl Velocity {
    pub(crate) fn parse(model: &Json, spec: &Json) -> Result<Self, &'static str> {
        if spec.field("schema")?.text()? == "solidcraft_interval_endpoint_velocity_v1" {
            spec.keys(&["schema", "bore", "endpoint"])?;
            if model.field("schema")?.text()? != crate::interval_bores_json::FAMILY {
                return Err("endpoint direction requires a certified interval bore model");
            }
            let recipe = crate::interval_bores_json::parse(model)?;
            let id = spec.field("bore")?.text()?;
            let end = match spec.field("endpoint")?.text()? {
                "lower" => 0,
                "upper" => 1,
                _ => return Err("unknown interval endpoint"),
            };
            let bore = recipe
                .bores
                .iter()
                .find(|b| b.bore.id == id)
                .ok_or("unknown interval bore")?;
            let z = bore.z_interval_mm[end];
            let h = recipe.height_mm;
            if z <= 0 || z >= h {
                return Err("external endpoint requires event analysis");
            }
            let cuts: Vec<_> = recipe.bores.iter().flat_map(|b| b.z_interval_mm).collect();
            if cuts.iter().filter(|&&v| v == z).count() != 1 {
                return Err("coincident endpoint requires event analysis");
            }
            let lower = cuts.iter().copied().filter(|&v| v < z).max().unwrap_or(0);
            let upper = cuts.iter().copied().filter(|&v| v > z).min().unwrap_or(h);
            let axis = [
                f64::from(recipe.top_shift_mm[0]) / f64::from(h),
                f64::from(recipe.top_shift_mm[1]) / f64::from(h),
                1.,
            ];
            return Ok(Self {
                linear: [[0.; 3]; 3],
                offset: [0.; 3],
                exact: None,
                endpoint: Some(EndpointFlow {
                    lower: lower.into(),
                    moving: z.into(),
                    upper: upper.into(),
                    axis,
                }),
            });
        }
        if spec.field("schema")?.text()? == "solidcraft_affine_velocity_v2" {
            spec.keys(&["schema", "affine_mm"])?;
            let map = crate::prescribed_motion::Motion::parse(spec.field("affine_mm")?)?;
            let coefficients = map.as_f64();
            return Ok(Self {
                linear: coefficients.map(|r| [r[1], r[2], r[3]]),
                offset: coefficients.map(|r| r[0]),
                exact: Some(map),
                endpoint: None,
            });
        }
        spec.keys(&["schema", "linear", "offset_mm"])?;
        if spec.field("schema")?.text()? != "solidcraft_affine_velocity_v1" {
            return Err("unknown mechanical velocity schema");
        }
        let [a, b, c] = spec.field("linear")?.items()? else {
            return Err("expected three velocity rows");
        };
        let row = crate::mechanical_job::vector;
        Ok(Self {
            linear: [row(a)?, row(b)?, row(c)?],
            offset: row(spec.field("offset_mm")?)?,
            exact: None,
            endpoint: None,
        })
    }

    pub(crate) fn divergence(&self, x: Vector) -> f64 {
        if let Some(cap) = &self.endpoint {
            return cap.slope(x[2]);
        }
        self.linear[0][0] + self.linear[1][1] + self.linear[2][2]
    }

    pub(crate) fn at(&self, x: Vector) -> Vector {
        if let Some(cap) = &self.endpoint {
            return cap.axis.map(|v| v * cap.amplitude(x[2]));
        }
        let linear = self.linear_velocity(x);
        std::array::from_fn(|i| linear[i] + self.offset[i])
    }

    /// The oriented area vector transforms by the cofactor of the affine map.
    pub(crate) fn oriented_area_rate(&self, x: Vector, n: Vector) -> Vector {
        let trace = self.divergence(x);
        if let Some(cap) = &self.endpoint {
            let dot: f64 = cap.axis.into_iter().zip(n).map(|(a, b)| a * b).sum();
            return [trace * n[0], trace * n[1], trace * (n[2] - dot)];
        }
        std::array::from_fn(|i| {
            trace * n[i] - (0..3).map(|k| self.linear[k][i] * n[k]).sum::<f64>()
        })
    }

    fn linear_velocity(&self, x: Vector) -> Vector {
        self.linear
            .map(|r| r.into_iter().zip(x).map(|(a, b)| a * b).sum())
    }
}

pub(crate) struct SensitivityJob<'a> {
    pub space: &'a AnalysisSpace,
    pub force: &'a [f64],
    pub force_rate: &'a [f64],
    pub includes_body_load: bool,
    pub lifts: [&'a Lift; 3],
    pub support_rates: Option<[Vec<LiftRate>; 3]>,
    pub material: MaterialJet,
    pub order: usize,
    pub relative_solve_tolerance: f64,
}

pub(crate) fn solve(velocity: &Velocity, data: SensitivityJob<'_>) -> Result<Json, &'static str> {
    let space = data.space;
    let mut point_velocities = vec![[0.; 3]; space.ndof];
    let mut jets = Vec::with_capacity(space.blocks.len());
    for block in &space.blocks {
        let mut local = Vec::with_capacity(block.controls.len());
        for (p, &i) in block.controls.iter().zip(&block.dofs) {
            let weight = p[3] as f64;
            let physical = std::array::from_fn(|c| p[c] as f64 / weight);
            point_velocities[i] = velocity.at(physical);
            // Drop only the global translation from the stiffness jet.
            let [x, y, z] = if velocity.endpoint.is_some() {
                point_velocities[i].map(|v| v * weight)
            } else {
                velocity.linear_velocity(std::array::from_fn(|c| p[c] as f64))
            };
            local.push([x, y, z, 0.]);
        }
        jets.push(local);
    }
    let k = assemble_sensitivity(space, &jets, data.material, data.order)?;
    let zero = vec![LiftRate::default(); space.ndof];
    let rates = data
        .support_rates
        .as_ref()
        .map_or([zero.as_slice(); 3], |rates| {
            rates.each_ref().map(Vec::as_slice)
        });
    let result = differentiate(
        MechanicalDerivative {
            stiffness: &k.stiffness,
            stiffness_rate: &k.derivative,
            force: data.force,
            force_rate: data.force_rate,
            lifts: data.lifts,
            lift_rates: rates,
            objective_gradient: data.force,
            objective_partial: 0.,
            fixed_regime: true,
        },
        data.relative_solve_tolerance,
        20000,
    )?;
    let dot = |a: &[f64], b: &[f64]| a.iter().zip(b).map(|(x, y)| x * y).sum::<f64>();
    let explicit = dot(data.force_rate, &result.primal.displacement_mm);
    let compliance = dot(data.force, &result.primal.displacement_mm);
    let num = |v: f64| Json::Number(v.to_string());
    let numbers = |v: &[f64]| Json::Array(v.iter().copied().map(num).collect());
    let resultant =
        std::array::from_fn::<_, 3, _>(|c| data.force_rate.iter().skip(c).step_by(3).sum());
    let mut report = object([
        (
            "schema",
            Json::String("solidcraft_affine_job_derivative_v1".into()),
        ),
        (
            "parameter",
            Json::String("dimensionless s in x(s)=x+s(Bx+b)".into()),
        ),
        ("fixed_regime", Json::Bool(true)),
        ("volume_mm3", num(k.stiffness.integrated_volume)),
        ("volume_rate_mm3", num(k.integrated_volume_derivative_mm3)),
        ("compliance_n_mm", num(compliance)),
        (
            "compliance_rate_tangent_n_mm",
            num(result.objective_tangent + explicit),
        ),
        (
            "compliance_rate_adjoint_n_mm",
            num(result.objective_adjoint + explicit),
        ),
        (
            if data.includes_body_load {
                "applied_load_explicit_work_rate_n_mm"
            } else {
                "traction_work_explicit_rate_n_mm"
            },
            num(explicit),
        ),
        (
            "physical_displacement_rate_mm_xyz",
            numbers(&result.physical_tangent),
        ),
        (
            "geometry_velocity_mm_xyz",
            Json::Array(point_velocities.iter().map(|v| numbers(v)).collect()),
        ),
        ("nodal_force_rate_n_xyz", numbers(data.force_rate)),
        ("applied_resultant_rate_n", numbers(&resultant)),
        (
            "relative_residuals",
            numbers(&[
                result.primal.relative_residual,
                result.tangent_relative_residual,
                result.adjoint_relative_residual,
            ]),
        ),
    ]);
    if let (Some(cap), Json::Object(fields)) = (&velocity.endpoint, &mut report) {
        fields.insert(
            "schema".into(),
            Json::String("solidcraft_interval_endpoint_derivative_v1".into()),
        );
        fields.insert(
            "parameter".into(),
            Json::String("positive endpoint height change in mm".into()),
        );
        fields.insert(
            "active_planes_mm".into(),
            numbers(&[cap.lower, cap.moving, cap.upper]),
        );
        fields.insert("material_axis_per_mm".into(), numbers(&cap.axis));
    }
    if let (Some(_), Json::Object(fields)) = (velocity.exact, &mut report) {
        let u = &result.primal.displacement_mm;
        let du = &result.physical_tangent;
        let ku = k.stiffness.apply(u)?;
        let dku = k.derivative.apply(u)?;
        let kdu = k.stiffness.apply(du)?;
        let reaction: Vec<_> = ku.iter().zip(data.force).map(|(a, b)| a - b).collect();
        let reaction_rate: Vec<_> = dku
            .iter()
            .zip(kdu)
            .zip(data.force_rate)
            .map(|((a, b), c)| a + b - c)
            .collect();
        let energy_shape = 0.5 * dot(u, &dku);
        let reaction_lift_work: f64 = reaction
            .as_chunks::<3>()
            .0
            .iter()
            .enumerate()
            .map(|(i, r)| (0..3).map(|c| r[c] * rates[c][i].offset).sum::<f64>())
            .sum();
        let strain_rate = energy_shape + dot(&ku, du);
        let adjoint_strain_rate = energy_shape + result.objective_adjoint + reaction_lift_work;
        let support_rate = dot(&reaction_rate, u) + dot(&reaction, du);
        let applied_rate = result.objective_tangent + explicit;
        fields.insert(
            "schema".into(),
            Json::String("solidcraft_affine_job_derivative_v2".into()),
        );
        fields.remove("compliance_n_mm");
        fields.remove("compliance_rate_tangent_n_mm");
        fields.remove("compliance_rate_adjoint_n_mm");
        for (key, value) in [
            ("applied_work_n_mm", compliance),
            ("applied_work_rate_tangent_n_mm", applied_rate),
            (
                "applied_work_rate_adjoint_n_mm",
                result.objective_adjoint + explicit,
            ),
            ("strain_energy_n_mm", result.primal.strain_energy_n_mm),
            ("strain_energy_rate_tangent_n_mm", strain_rate),
            ("strain_energy_rate_adjoint_n_mm", adjoint_strain_rate),
            ("support_work_n_mm", dot(&reaction, u)),
            ("support_work_rate_n_mm", support_rate),
            (
                "work_rate_closure_n_mm",
                applied_rate + support_rate - 2. * strain_rate,
            ),
        ] {
            if !value.is_finite() {
                return Err("nonfinite material work derivative");
            }
            fields.insert(key.into(), num(value));
        }
    }
    Ok(report)
}
