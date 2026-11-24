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
        if !matches!(degrees, [1, 1, 1] | [2, 1, 1]) {
            return Err(MapError::UnsupportedDegrees);
        }
        Ok(Self { degrees, controls })
    }

    pub fn evaluate(&self, [u, v, w]: [f64; 3]) -> Result<Evaluation, MapError> {
        if [u, v, w].iter().any(|x| !(0.0..=1.0).contains(x)) {
            return Err(MapError::OutsideDomain);
        }
        let a = 1.0 - u;
        let (bu, du) = if self.degrees[0] == 1 {
            ([a, u, 0.0], [-1.0, 1.0, 0.0])
        } else {
            (
                [a * a, 2.0 * a * u, u * u],
                [-2.0 * a, 2.0 * (a - u), 2.0 * u],
            )
        };
        let mut h = [[0.0; 4]; 4];
        for (index, p) in self.controls.iter().enumerate() {
            let (i, j, k) = (index / 4, (index / 2) % 2, index % 2);
            let bv = [1.0 - v, v][j];
            let bw = [1.0 - w, w][k];
            let factors = [
                bu[i] * bv * bw,
                du[i] * bv * bw,
                bu[i] * [-1.0, 1.0][j] * bw,
                bu[i] * bv * [-1.0, 1.0][k],
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
