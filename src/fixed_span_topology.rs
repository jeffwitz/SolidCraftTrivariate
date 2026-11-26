//! Verify the original full-face neighbors after axial knot insertion within charts.
use crate::{
    bracket::BracketBlock,
    incidence::{enumerate_faces, BlockSide},
};
use std::collections::BTreeSet;

fn parent_side(side: &BlockSide) -> Result<(BlockSide, usize), &'static str> {
    let (owner, index) = side.block.rsplit_once('#').ok_or("missing span index")?;
    let mut result = side.clone();
    result.block = owner.to_string();
    Ok((result, index.parse().map_err(|_| "invalid span index")?))
}

pub(crate) fn verify(parent: &[BracketBlock], cells: &[BracketBlock]) -> Result<(), &'static str> {
    let source = enumerate_faces(parent)?;
    let edited = enumerate_faces(cells)?;
    let mut pairs = BTreeSet::new();
    let mut knots = BTreeSet::new();
    for e in edited.interfaces {
        let (a, i) = parent_side(&e.first)?;
        let (b, j) = parent_side(&e.second)?;
        if a.block == b.block {
            if i.abs_diff(j) != 1
                || a.axis != b.axis
                || a.upper != (i < j)
                || b.upper != (j < i)
                || e.orientation.swap
                || e.orientation.flip_first
                || e.orientation.flip_second
                || !knots.insert((a.block, i.min(j)))
            {
                return Err("incompatible consecutive knot spans");
            }
        } else {
            pairs.insert((
                a,
                b,
                e.orientation.swap,
                e.orientation.flip_first,
                e.orientation.flip_second,
            ));
        }
    }
    let expected = source
        .interfaces
        .into_iter()
        .map(|e| {
            (
                e.first,
                e.second,
                e.orientation.swap,
                e.orientation.flip_first,
                e.orientation.flip_second,
            )
        })
        .collect();
    let exterior = edited
        .exterior
        .iter()
        .map(|s| parent_side(s).map(|(p, _)| p))
        .collect::<Result<BTreeSet<_>, _>>()?;
    if parent.len().checked_add(knots.len()) != Some(cells.len())
        || pairs != expected
        || exterior != source.exterior.into_iter().collect()
    {
        return Err("edge treatment changed logical face incidence");
    }
    Ok(())
}
