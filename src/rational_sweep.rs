//! Parametric rational sweep with a stable physical wedge and adaptive charts.
//!
//! The physical domain is the extruded annular wedge described by
//! `theta=2*atan(t*u)`, radius `r=2-v`, and `z=h*w`. Its physical signature is
//! independent of whether the current representation uses one chart or two.
//! Controls and their analytic jets are materialized in binary64; this module
//! does not provide a rounding enclosure or a geometry-adoption certificate.

/// Why a rational sweep request was refused.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SweepError {
    NonFiniteParameter,
    ParameterOutOfDomain,
    NumericalRange,
}

/// Current chart count, which may change while the physical wedge stays smooth.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SweepRepresentation {
    SingleBlock,
    TwoBlocks,
}

/// Materialization status of the binary64 control net.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MaterializationStatus {
    /// Controls are rounded binary64 values; no error enclosure is certified.
    Binary64Uncertified,
}

/// Stable physical-domain signature for an extruded annular wedge.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SweepDomainSignature {
    pub family: &'static str,
    pub start_angle_rad: f64,
    pub end_angle_rad: f64,
    pub inner_radius_mm: f64,
    pub outer_radius_mm: f64,
    pub height_mm: f64,
}

/// One local tensor-product chart, restricted from the global master parameter.
#[derive(Debug, Clone, PartialEq)]
pub struct SweepBlock {
    pub id: &'static str,
    pub parameter_interval: [f64; 2],
    pub degrees: [usize; 3],
    /// Rounded homogeneous `(wx,wy,wz,w)` controls; `w` index varies fastest.
    pub controls: Vec<[f64; 4]>,
    /// Analytic control jet at fixed local parameters; absent at the chart event.
    pub dcontrols_dt: Option<Vec<[f64; 4]>>,
    /// Analytic control jet for the height parameter.
    pub dcontrols_dh: Vec<[f64; 4]>,
}

/// A materialized sweep, its representation signature, and physical volume jet.
#[derive(Debug, Clone, PartialEq)]
pub struct RationalSweep {
    pub t: f64,
    pub h: f64,
    pub domain: SweepDomainSignature,
    pub representation: SweepRepresentation,
    pub representation_event: bool,
    pub materialization: MaterializationStatus,
    pub blocks: Vec<SweepBlock>,
    pub volume_mm3: f64,
    pub dvolume_dt_mm3: f64,
    pub dvolume_dh_mm2: f64,
}

/// Physical map and analytic derivatives with respect to master parameters.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SweepEvaluation {
    pub position_mm: [f64; 3],
    /// Derivative columns with respect to master `(u,v,w)`.
    pub parameter_derivatives: [[f64; 3]; 3],
    pub dposition_dt_mm: [f64; 3],
    pub dposition_dh_mm: [f64; 3],
}

const T_MIN: f64 = 0.25;
const T_MAX: f64 = 4.0;
const H_MIN: f64 = 0.5;
const H_MAX: f64 = 32.0;
const DEGREES: [usize; 3] = [2, 1, 1];

fn validate_parameters(t: f64, h: f64) -> Result<(), SweepError> {
    if !t.is_finite() || !h.is_finite() {
        return Err(SweepError::NonFiniteParameter);
    }
    if !(T_MIN..=T_MAX).contains(&t) || !(H_MIN..=H_MAX).contains(&h) {
        return Err(SweepError::ParameterOutOfDomain);
    }
    Ok(())
}

type Net = Vec<[f64; 4]>;

fn base_controls(t: f64, h: f64) -> (Net, Net, Net) {
    // Bernstein controls for (1-(tu)^2, 2tu, 1+(tu)^2), with derivatives in t.
    let arc = [
        ([1.0, 0.0, 1.0], [0.0, 0.0, 0.0]),
        ([1.0, t, 1.0], [0.0, 1.0, 0.0]),
        (
            [1.0 - t * t, 2.0 * t, 1.0 + t * t],
            [-2.0 * t, 2.0, 2.0 * t],
        ),
    ];
    let mut controls = Vec::with_capacity(12);
    let mut dt = Vec::with_capacity(12);
    let mut dh = Vec::with_capacity(12);
    for (x, dx) in arc {
        for radius in [2.0, 1.0] {
            for z_fraction in [0.0, 1.0] {
                controls.push([radius * x[0], radius * x[1], h * z_fraction * x[2], x[2]]);
                dt.push([
                    radius * dx[0],
                    radius * dx[1],
                    h * z_fraction * dx[2],
                    dx[2],
                ]);
                dh.push([0.0, 0.0, z_fraction * x[2], 0.0]);
            }
        }
    }
    (controls, dt, dh)
}

fn bernstein2(u: f64) -> [f64; 3] {
    let a = 1.0 - u;
    [a * a, 2.0 * a * u, u * u]
}

fn bernstein2_derivative(u: f64) -> [f64; 3] {
    [-2.0 * (1.0 - u), 2.0 - 4.0 * u, 2.0 * u]
}

fn restrict_quadratic(
    source: &[[f64; 4]],
    interval: [f64; 2],
) -> Result<Vec<[f64; 4]>, SweepError> {
    if source.len() != 12 {
        return Err(SweepError::NumericalRange);
    }
    let [a, b] = interval;
    let delta = b - a;
    let left = bernstein2(a);
    let left_derivative = bernstein2_derivative(a);
    let right = bernstein2(b);
    let mut result = Vec::with_capacity(12);
    let coefficients = [
        left,
        std::array::from_fn(|i| left[i] + 0.5 * delta * left_derivative[i]),
        right,
    ];
    for c in coefficients {
        for radial in 0..2 {
            for height in 0..2 {
                let mut point = [0.0; 4];
                for (source_angle, factor) in c.iter().enumerate() {
                    let source_index = (source_angle * 2 + radial) * 2 + height;
                    for component in 0..4 {
                        point[component] += factor * source[source_index][component];
                    }
                }
                result.push(point);
            }
        }
    }
    if result.iter().flatten().any(|value| !value.is_finite()) {
        return Err(SweepError::NumericalRange);
    }
    Ok(result)
}

fn restricted_block(
    id: &'static str,
    interval: [f64; 2],
    base: &[[f64; 4]],
    base_dt: &[[f64; 4]],
    base_dh: &[[f64; 4]],
    event: bool,
) -> Result<SweepBlock, SweepError> {
    Ok(SweepBlock {
        id,
        parameter_interval: interval,
        degrees: DEGREES,
        controls: restrict_quadratic(base, interval)?,
        dcontrols_dt: if event {
            None
        } else {
            Some(restrict_quadratic(base_dt, interval)?)
        },
        dcontrols_dh: restrict_quadratic(base_dh, interval)?,
    })
}

/// Construct the rational annular sweep and its analytic control and volume jets.
pub fn build_rational_sweep(t: f64, h: f64) -> Result<RationalSweep, SweepError> {
    validate_parameters(t, h)?;
    let threshold = 1.0 + 2.0_f64.sqrt();
    let event = t == threshold;
    let split = t > threshold;
    let (base, base_dt, base_dh) = base_controls(t, h);
    let (representation, blocks) = if split {
        (
            SweepRepresentation::TwoBlocks,
            vec![
                restricted_block("first", [0.0, 0.5], &base, &base_dt, &base_dh, event)?,
                restricted_block("second", [0.5, 1.0], &base, &base_dt, &base_dh, event)?,
            ],
        )
    } else {
        (
            SweepRepresentation::SingleBlock,
            vec![restricted_block(
                "whole",
                [0.0, 1.0],
                &base,
                &base_dt,
                &base_dh,
                event,
            )?],
        )
    };

    let volume = 3.0 * h * t.atan();
    let dvolume_dt = 3.0 * h / (1.0 + t * t);
    let dvolume_dh = 3.0 * t.atan();
    if !volume.is_finite() || !dvolume_dt.is_finite() || !dvolume_dh.is_finite() {
        return Err(SweepError::NumericalRange);
    }
    Ok(RationalSweep {
        t,
        h,
        domain: SweepDomainSignature {
            family: "extruded_annular_wedge_v1",
            start_angle_rad: 0.0,
            end_angle_rad: 2.0 * t.atan(),
            inner_radius_mm: 1.0,
            outer_radius_mm: 2.0,
            height_mm: h,
        },
        representation,
        representation_event: event,
        materialization: MaterializationStatus::Binary64Uncertified,
        blocks,
        volume_mm3: volume,
        dvolume_dt_mm3: dvolume_dt,
        dvolume_dh_mm2: dvolume_dh,
    })
}

/// Evaluate the physical map and derivatives at global master parameters.
pub fn evaluate_sweep(t: f64, h: f64, [u, v, w]: [f64; 3]) -> Result<SweepEvaluation, SweepError> {
    validate_parameters(t, h)?;
    if [u, v, w].iter().any(|value| !value.is_finite()) {
        return Err(SweepError::NonFiniteParameter);
    }
    if [u, v, w].iter().any(|value| !(0.0..=1.0).contains(value)) {
        return Err(SweepError::ParameterOutOfDomain);
    }
    let q = t * u;
    let q2 = q * q;
    let denominator = 1.0 + q2;
    let denominator2 = denominator * denominator;
    let radius = 2.0 - v;
    let cos = (1.0 - q2) / denominator;
    let sin = 2.0 * q / denominator;
    let evaluation = SweepEvaluation {
        position_mm: [radius * cos, radius * sin, h * w],
        parameter_derivatives: [
            [
                -4.0 * radius * q * t / denominator2,
                2.0 * radius * (1.0 - q2) * t / denominator2,
                0.0,
            ],
            [-cos, -sin, 0.0],
            [0.0, 0.0, h],
        ],
        dposition_dt_mm: [
            -4.0 * radius * q * u / denominator2,
            2.0 * radius * (1.0 - q2) * u / denominator2,
            0.0,
        ],
        dposition_dh_mm: [0.0, 0.0, w],
    };
    if evaluation
        .position_mm
        .iter()
        .chain(evaluation.parameter_derivatives.iter().flatten())
        .chain(evaluation.dposition_dt_mm.iter())
        .chain(evaluation.dposition_dh_mm.iter())
        .any(|value| !value.is_finite())
    {
        return Err(SweepError::NumericalRange);
    }
    Ok(evaluation)
}
