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
                        for i in 0..3 {
                            let j = (i + 1) % 3;
                            let k = (i + 2) % 3;
                            out.resultant[i] += df[i];
                            out.moment[i] += x[j] * df[k] - x[k] * df[j];
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
                            let div = flow.divergence();
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
    pub resultant: Vector,
    pub moment: Vector,
}
