//! Ordered additive and subtractive cartesian boxes with exact integer planes.
//! Material connectivity and local octant links are checked before charts are returned.
use crate::bracket::BracketBlock;
use std::collections::{BTreeSet, VecDeque};

pub type Bounds = [[i32; 2]; 3];
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BoxOp {
    Add,
    Cut,
}
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct BoxFeature {
    pub id: String,
    pub operation: BoxOp,
    pub bounds_mm: Bounds,
}
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct OrthogonalProgram {
    pub base_mm: Bounds,
    pub features: Vec<BoxFeature>,
}
fn valid_bounds(bounds: &Bounds) -> bool {
    bounds
        .iter()
        .all(|&[lo, hi]| lo >= -64 && lo < hi && hi <= 64)
}
fn encloses(a: &Bounds, b: &Bounds) -> bool {
    (0..3).all(|j| a[j][0] <= b[j][0] && b[j][1] <= a[j][1])
}
impl OrthogonalProgram {
    pub fn validate_parameters(&self) -> Result<(), &'static str> {
        if !valid_bounds(&self.base_mm) || self.features.len() > 8 {
            return Err("orthogonal_program_range");
        }
        let mut ids = BTreeSet::new();
        for f in &self.features {
            if f.id.is_empty()
                || f.id.len() > 48
                || f.id == "base"
                || !f.id.bytes().all(|x| x.is_ascii_alphanumeric() || x == b'_')
                || !ids.insert(&f.id)
            {
                return Err("orthogonal_feature_id");
            }
            if !valid_bounds(&f.bounds_mm) {
                return Err("orthogonal_feature_bounds");
            }
        }
        Ok(())
    }
}

struct Grid {
    axes: [Vec<i32>; 3],
    dims: [usize; 3],
    solid: Vec<bool>,
}
impl Grid {
    fn new(program: &OrthogonalProgram) -> Result<Self, &'static str> {
        program.validate_parameters()?;
        let axes: [Vec<i32>; 3] = std::array::from_fn(|j| {
            let mut x = BTreeSet::from(program.base_mm[j]);
            for f in &program.features {
                x.extend(f.bounds_mm[j]);
            }
            x.into_iter().collect()
        });
        let dims = axes.each_ref().map(|a| a.len() - 1);
        let mut grid = Self {
            axes,
            dims,
            solid: Vec::with_capacity(dims.iter().product()),
        };
        for i in 0..dims[0] {
            for j in 0..dims[1] {
                for k in 0..dims[2] {
                    let b = grid.bounds([i, j, k]);
                    let mut inside = encloses(&program.base_mm, &b);
                    for f in &program.features {
                        if encloses(&f.bounds_mm, &b) {
                            inside = f.operation == BoxOp::Add;
                        }
                    }
                    grid.solid.push(inside);
                }
            }
        }
        Ok(grid)
    }
    fn bounds(&self, id: [usize; 3]) -> Bounds {
        std::array::from_fn(|k| [self.axes[k][id[k]], self.axes[k][id[k] + 1]])
    }
    fn index(&self, p: [usize; 3]) -> usize {
        (p[0] * self.dims[1] + p[1]) * self.dims[2] + p[2]
    }
    fn coordinates(&self, n: usize) -> [usize; 3] {
        [
            n / (self.dims[1] * self.dims[2]),
            n / self.dims[2] % self.dims[1],
            n % self.dims[2],
        ]
    }
    fn occupied(&self, p: [isize; 3]) -> bool {
        if (0..3).any(|a| p[a] < 0 || p[a] >= self.dims[a] as isize) {
            return false;
        }
        self.solid[self.index(p.map(|x| x as usize))]
    }
    fn check_connectivity(&self) -> Result<(), &'static str> {
        let seed = self
            .solid
            .iter()
            .position(|x| *x)
            .ok_or("orthogonal_empty_material")?;
        let mut seen = vec![false; self.solid.len()];
        seen[seed] = true;
        let mut found = 0;
        let mut queue = VecDeque::from([seed]);
        while let Some(i) = queue.pop_front() {
            found += 1;
            let p = self.coordinates(i);
            for a in 0..3 {
                for d in [-1, 1] {
                    let mut q = p.map(|x| x as isize);
                    q[a] += d;
                    if self.occupied(q) {
                        let j = self.index(q.map(|x| x as usize));
                        if !seen[j] {
                            seen[j] = true;
                            queue.push_back(j);
                        }
                    }
                }
            }
        }
        if found != self.solid.iter().filter(|x| **x).count() {
            return Err("orthogonal_disconnected_material");
        }
        Ok(())
    }
    fn check_links(&self) -> Result<(), &'static str> {
        for i in 0..=self.dims[0] {
            for j in 0..=self.dims[1] {
                for k in 0..=self.dims[2] {
                    let vertex = [i, j, k];
                    let mask = (0..8).fold(0u8, |mask, octant| {
                        let point = std::array::from_fn(|a| {
                            vertex[a] as isize - 1 + ((octant >> a) & 1) as isize
                        });
                        mask | (u8::from(self.occupied(point)) << octant)
                    });
                    if mask == 0 || mask == 255 {
                        continue;
                    }
                    // Each face of the octant cube is the cyclic transverse link of a grid edge.
                    for a in 0..3 {
                        let b = (a + 1) % 3;
                        let c = (a + 2) % 3;
                        for s in 0..2 {
                            let base = s << a;
                            let values = [
                                base,
                                base | (1 << b),
                                base | (1 << b) | (1 << c),
                                base | (1 << c),
                            ]
                            .map(|v| (mask >> v) & 1);
                            if values == [1, 0, 1, 0] || values == [0, 1, 0, 1] {
                                return Err("orthogonal_nonmanifold_edge");
                            }
                        }
                    }
                    if !octants_connected(mask) || !octants_connected(!mask) {
                        return Err("orthogonal_nonmanifold_vertex");
                    }
                }
            }
        }
        Ok(())
    }
}
fn octants_connected(mask: u8) -> bool {
    if mask == 0 {
        return true;
    }
    let mut seen = 1u8 << mask.trailing_zeros();
    loop {
        let mut next = seen;
        for i in 0..8 {
            if seen & (1 << i) != 0 {
                for a in 0..3 {
                    next |= (1 << (i ^ (1 << a))) & mask;
                }
            }
        }
        if next == seen {
            return seen == mask;
        }
        seen = next;
    }
}
pub fn build_orthogonal_program(
    program: &OrthogonalProgram,
) -> Result<Vec<BracketBlock>, &'static str> {
    let grid = Grid::new(program)?;
    grid.check_links()?;
    grid.check_connectivity()?;
    let mut blocks = Vec::new();
    for (i, &solid) in grid.solid.iter().enumerate() {
        if solid {
            let [x, y, z] = grid.bounds(grid.coordinates(i));
            let id = format!("box_{}_{}_{}_{}_{}_{}", x[0], x[1], y[0], y[1], z[0], z[1]);
            let mut controls = Vec::with_capacity(8);
            for a in x {
                for b in y {
                    for c in z {
                        controls.push([i64::from(a), i64::from(b), i64::from(c), 1]);
                    }
                }
            }
            blocks.push(BracketBlock {
                id,
                degrees: [1, 1, 1],
                controls,
            });
        }
    }
    Ok(blocks)
}

#[cfg(test)]
mod tests {
    use super::*;
    fn composed(bounds_mm: Bounds) -> OrthogonalProgram {
        OrthogonalProgram {
            base_mm: [[0, 4]; 3],
            features: vec![BoxFeature {
                id: "addition".into(),
                operation: BoxOp::Add,
                bounds_mm,
            }],
        }
    }
    #[test]
    fn contacts_have_different_link_witnesses() {
        assert_eq!(
            build_orthogonal_program(&composed([[4, 6], [4, 6], [1, 3]])),
            Err("orthogonal_nonmanifold_edge")
        );
        assert_eq!(
            build_orthogonal_program(&composed([[4, 6]; 3])),
            Err("orthogonal_nonmanifold_vertex")
        );
        assert_eq!(
            build_orthogonal_program(&composed([[5, 6]; 3])),
            Err("orthogonal_disconnected_material")
        );
        assert!(build_orthogonal_program(&composed([[4, 6], [1, 3], [1, 3]])).is_ok());
    }
}
