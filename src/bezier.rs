//! Rational tensor-product Bézier maps, with the third control index varying fastest.

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MapError {
    UnsupportedDegrees,
    WrongControlCount,
    InvalidControl,
    OutsideDomain,
    NumericalRange,
}

impl std::fmt::Display for MapError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            Self::UnsupportedDegrees => "unsupported degree family",
            Self::WrongControlCount => "incorrect number of homogeneous controls",
            Self::InvalidControl => "controls must be finite with positive weights",
            Self::OutsideDomain => "parameters must lie in the closed unit cube",
            Self::NumericalRange => "result outside the supported numerical range",
        })
    }
}

impl std::error::Error for MapError {}

#[derive(Debug, Clone)]
pub struct RationalBezierMap {
    degrees: [usize; 3],
    controls: Vec<[f64; 4]>,
}

#[derive(Debug, Clone, Copy)]
pub struct Evaluation {
    pub position: [f64; 3],
    /// Physical derivatives, in (du, dv, dw) order.
    pub jacobian_columns: [[f64; 3]; 3],
}

impl RationalBezierMap {
    /// Each control is (weight*x, weight*y, weight*z, weight).
    pub fn new(degrees: [usize; 3], controls: Vec<[f64; 4]>) -> Result<Self, MapError> {
        let count = degrees
            .iter()
            .try_fold(1usize, |n, d| n.checked_mul(d.checked_add(1)?));
        if count != Some(controls.len()) {
            return Err(MapError::WrongControlCount);
        }
        if controls
            .iter()
            .any(|p| p[3] <= 0.0 || p.iter().any(|c| !c.is_finite()))
        {
            return Err(MapError::InvalidControl);
        }
        if degrees.iter().any(|d| !(1..=3).contains(d)) {
            return Err(MapError::UnsupportedDegrees);
        }
        Ok(Self { degrees, controls })
    }

    pub fn evaluate(&self, [u, v, w]: [f64; 3]) -> Result<Evaluation, MapError> {
        if [u, v, w].iter().any(|x| !(0.0..=1.0).contains(x)) {
            return Err(MapError::OutsideDomain);
        }
        let (bu, du) = basis(self.degrees[0], u);
        let (bv, dv) = basis(self.degrees[1], v);
        let (bw, dw) = basis(self.degrees[2], w);
        let w_count = self.degrees[2] + 1;
        let v_count = self.degrees[1] + 1;
        let mut h = [[0.0; 4]; 4];
        for (index, p) in self.controls.iter().enumerate() {
            let (i, j, k) = (
                index / (w_count * v_count),
                (index / w_count) % v_count,
                index % w_count,
            );
            let factors = [
                bu[i] * bv[j] * bw[k],
                du[i] * bv[j] * bw[k],
                bu[i] * dv[j] * bw[k],
                bu[i] * bv[j] * dw[k],
            ];
            for (row, factor) in h.iter_mut().zip(factors) {
                for (sum, coord) in row.iter_mut().zip(p) {
                    *sum += factor * coord;
                }
            }
        }
        let weight = h[0][3];
        if weight <= 0.0 || !h.iter().flatten().all(|x| x.is_finite()) {
            return Err(MapError::NumericalRange);
        }
        let position = std::array::from_fn(|c| h[0][c] / weight);
        let jacobian_columns: [[f64; 3]; 3] = std::array::from_fn(|axis| {
            let row = h[axis + 1];
            std::array::from_fn(|c| (row[c] - position[c] * row[3]) / weight)
        });
        if !position
            .iter()
            .chain(jacobian_columns.iter().flatten())
            .all(|x| x.is_finite())
        {
            return Err(MapError::NumericalRange);
        }
        Ok(Evaluation {
            position,
            jacobian_columns,
        })
    }
}

pub(crate) fn basis(degree: usize, t: f64) -> ([f64; 4], [f64; 4]) {
    let a = 1. - t;
    match degree {
        1 => ([a, t, 0., 0.], [-1., 1., 0., 0.]),
        2 => (
            [a * a, 2. * a * t, t * t, 0.],
            [-2. * a, 2. * (a - t), 2. * t, 0.],
        ),
        3 => (
            [a * a * a, 3. * a * a * t, 3. * a * t * t, t * t * t],
            [
                -3. * a * a,
                3. * a * (a - 2. * t),
                3. * t * (2. * a - t),
                3. * t * t,
            ],
        ),
        _ => unreachable!("degree validated in RationalBezierMap::new"),
    }
}
