//! Integration of a spatial affine force density on the audited material volume.
use crate::{
    affine_job_sensitivity::Velocity,
    analysis_space::AnalysisSpace,
    elasticity::{basis_geometry, gauss, Vector},
    json::Json,
    mechanical_job::vector,
};

pub(crate) struct Density {
    constant: Vector,
    gradient: [Vector; 3],
}

impl Density {
    pub(crate) fn parse(spec: &Json) -> Result<Self, &'static str> {
        spec.keys(&["constant_n_per_mm3", "gradient_n_per_mm4"])?;
        let [a, b, c] = spec.field("gradient_n_per_mm4")?.items()? else {
            return Err("expected three body force gradient rows");
        };
        Ok(Self {
            constant: vector(spec.field("constant_n_per_mm3")?)?,
            gradient: [vector(a)?, vector(b)?, vector(c)?],
        })
    }

    fn variation(&self, x: Vector) -> Vector {
        self.gradient
            .map(|r| r.into_iter().zip(x).map(|(a, b)| a * b).sum())
    }

    pub(crate) fn integrate(
        &self,
        space: &AnalysisSpace,
        order: usize,
        velocity: Option<&Velocity>,
    ) -> Result<VolumeLoad, &'static str> {
        let q = gauss(order)?;
        let size = space
            .ndof
            .checked_mul(3)
            .ok_or("body load dimension overflow")?;
        let mut out = VolumeLoad {
            force: vec![0.; size],
            force_rate: velocity.map(|_| vec![0.; size]),
            volume: 0.,
            scale: 0.,
            moment_scale: 0.,
            resultant: [0.; 3],
            moment: [0.; 3],
        };
        for block in &space.blocks {
            if block.dofs.len() != 27 || block.dofs.iter().any(|&id| id >= space.ndof) {
                return Err("invalid body load block DOFs");
            }
            for &(u, a) in &q {
                for &(v, b) in &q {
                    for &(w, c) in &q {
                        let sample = basis_geometry(block, [u, v, w])?;
                        let measure = sample.determinant * a * b * c;
                        let x = sample.position;
                        let linear = self.variation(x);
                        let density: Vector = std::array::from_fn(|i| self.constant[i] + linear[i]);
                        let df = density.map(|f| f * measure);
                        out.volume += measure;
                        out.scale += df[0].hypot(df[1]).hypot(df[2]);
                        let moment: Vector = std::array::from_fn(|i| {
                            let j = (i + 1) % 3;
                            let k = (i + 2) % 3;
                            x[j] * df[k] - x[k] * df[j]
                        });
                        out.moment_scale += moment[0].hypot(moment[1]).hypot(moment[2]);
                        for i in 0..3 {
                            out.resultant[i] += df[i];
                            out.moment[i] += moment[i];
                        }
                        let distribute = |nodal: &mut [f64], load: Vector| {
                            for (&id, shape) in block.dofs.iter().zip(sample.values) {
                                for (i, f) in load.iter().enumerate() {
                                    nodal[3 * id + i] += shape * f;
                                }
                            }
                        };
                        distribute(&mut out.force, df);
                        if let (Some(flow), Some(rate)) = (velocity, &mut out.force_rate) {
                            let gradient_rate = self.variation(flow.at(x));
                            let div = flow.divergence(x);
                            distribute(
                                rate,
                                std::array::from_fn(|i| {
                                    measure * (gradient_rate[i] + div * density[i])
                                }),
                            );
                        }
                    }
                }
            }
        }
        if !out.scale.is_finite()
            || !out.moment_scale.is_finite()
            || !out.volume.is_finite()
            || out
                .force
                .iter()
                .chain(out.force_rate.iter().flatten())
                .chain(&out.resultant)
                .chain(&out.moment)
                .any(|x| !x.is_finite())
        {
            return Err("nonfinite integrated body load");
        }
        Ok(out)
    }
}

pub(crate) struct VolumeLoad {
    pub force: Vec<f64>,
    pub force_rate: Option<Vec<f64>>,
    pub volume: f64,
    pub scale: f64,
    moment_scale: f64,
    pub resultant: Vector,
    pub moment: Vector,
}

/// Adjacent global tensor Gauss rules for the volume right-hand side.
pub(crate) struct BodyQuadrature {
    min_order: usize,
    max_order: usize,
    relative_tolerance: f64,
    max_points: usize,
}

impl BodyQuadrature {
    pub(crate) fn parse(spec: &Json) -> Result<Self, &'static str> {
        spec.keys(&["min_order", "max_order", "relative_tolerance", "max_points"])?;
        let config = Self {
            min_order: spec.field("min_order")?.number()?,
            max_order: spec.field("max_order")?.number()?,
            relative_tolerance: spec.field("relative_tolerance")?.number()?,
            max_points: spec.field("max_points")?.number()?,
        };
        if !(2..=16).contains(&config.min_order)
            || config.max_order > 16
            || config.min_order >= config.max_order
            || !config.relative_tolerance.is_finite()
            || config.relative_tolerance <= 0.
            || config.relative_tolerance > 1.
            || config.max_points == 0
        {
            return Err("invalid body quadrature configuration");
        }
        Ok(config)
    }

    pub(crate) fn integrate(
        &self,
        density: &Density,
        space: &AnalysisSpace,
        velocity: Option<&Velocity>,
    ) -> Result<(VolumeLoad, Json), &'static str> {
        use crate::json::object;
        let n = |v: usize| Json::Number(v.to_string());
        let mut points = 0usize;
        let mut levels = Vec::new();
        let mut previous: Option<VolumeLoad> = None;
        let mut order = self.min_order;
        loop {
            let sample_count = space
                .blocks
                .len()
                .checked_mul(order.pow(3))
                .ok_or("body quadrature count overflow")?;
            points = points
                .checked_add(sample_count)
                .ok_or("body quadrature count overflow")?;
            if points > self.max_points {
                return Err("body quadrature point budget exceeded");
            }
            let integrated = density.integrate(space, order, velocity)?;
            let estimate = previous
                .as_ref()
                .map(|old| body_discrepancy(old, &integrated))
                .transpose()?;
            levels.push(object([
                ("order", n(order)),
                ("points", n(sample_count)),
                (
                    "relative_discrepancies",
                    estimate.map_or(Json::Null, |v| {
                        Json::Array(v.map(|x| Json::Number(x.to_string())).to_vec())
                    }),
                ),
            ]));
            if estimate.is_some_and(|values| values.iter().all(|v| *v <= self.relative_tolerance)) {
                return Ok((
                    integrated,
                    object([
                        ("selected_order", n(order)),
                        ("points", n(points)),
                        (
                            "metrics",
                            Json::Array(
                                ["volume", "nodal_force", "resultant", "moment"]
                                    .map(|v| Json::String(v.into()))
                                    .to_vec(),
                            ),
                        ),
                        ("history", Json::Array(levels)),
                    ]),
                ));
            }
            if order == self.max_order {
                return Err("body quadrature tolerance not met within order budget");
            }
            previous = Some(integrated);
            order = (order + 2).min(self.max_order);
        }
    }
}

fn body_discrepancy(old: &VolumeLoad, new: &VolumeLoad) -> Result<[f64; 4], &'static str> {
    fn relative(a: &[f64], b: &[f64], scale: f64) -> f64 {
        let delta = a.iter().zip(b).fold(0_f64, |s, (x, y)| s.hypot(x - y));
        let norm = b.iter().fold(0_f64, |s, x| s.hypot(*x));
        let denominator = if norm > 1e-12 * scale { norm } else { scale };
        if denominator == 0. && delta == 0. {
            0.
        } else {
            delta / denominator
        }
    }
    let values = [
        (old.volume - new.volume).abs() / new.volume,
        relative(&old.force, &new.force, new.scale),
        relative(&old.resultant, &new.resultant, new.scale),
        relative(&old.moment, &new.moment, new.moment_scale),
    ];
    if values.iter().any(|x| !x.is_finite()) {
        return Err("invalid body quadrature estimate");
    }
    Ok(values)
}
