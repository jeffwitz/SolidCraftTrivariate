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
    pub(crate) exact: Option<crate::prescribed_motion::Motion>,
}

impl Velocity {
    pub(crate) fn parse(spec: &Json) -> Result<Self, &'static str> {
        if spec.field("schema")?.text()? == "solidcraft_affine_velocity_v2" {
            spec.keys(&["schema", "affine_mm"])?;
            let map = crate::prescribed_motion::Motion::parse(spec.field("affine_mm")?)?;
            let coefficients = map.as_f64();
            return Ok(Self {
                linear: coefficients.map(|r| [r[1], r[2], r[3]]),
                offset: coefficients.map(|r| r[0]),
                exact: Some(map),
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
        })
    }

    pub(crate) fn divergence(&self) -> f64 {
        self.linear[0][0] + self.linear[1][1] + self.linear[2][2]
    }

    pub(crate) fn at(&self, x: Vector) -> Vector {
        let linear = self.linear_velocity(x);
        std::array::from_fn(|i| linear[i] + self.offset[i])
    }

    /// The oriented area vector transforms by the cofactor of the affine map.
    pub(crate) fn oriented_area_rate(&self, n: Vector) -> Vector {
        let trace = self.divergence();
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
}

pub(crate) fn solve(velocity: &Velocity, data: SensitivityJob<'_>) -> Result<Json, &'static str> {
    let space = data.space;
    let mut point_velocities = vec![[0.; 3]; space.ndof];
    let mut jets = Vec::with_capacity(space.blocks.len());
    for block in &space.blocks {
        let mut local = Vec::with_capacity(block.controls.len());
        for (p, &i) in block.controls.iter().zip(&block.dofs) {
            // A global translation cancels identically from all gradients and
            // integration measures. Omit it in K' and restore it in the velocity field.
            let [x, y, z] = velocity.linear_velocity(std::array::from_fn(|k| p[k] as f64));
            local.push([x, y, z, 0.]);
            point_velocities[i] =
                std::array::from_fn(|k| [x, y, z][k] / p[3] as f64 + velocity.offset[k]);
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
        1e-11,
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
