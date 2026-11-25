//! Proper signed-permutation frames for equal-radius orthogonal cylinder pairs.
use crate::{bracket::BracketBlock, steinmetz_core::steinmetz_box};

#[derive(Clone, Copy)]
struct Frame {
    world_axes: [usize; 3],
}

fn frame(axes: [usize; 2]) -> Result<Frame, &'static str> {
    if axes[0] >= 3 || axes[1] >= 3 || axes[0] == axes[1] {
        return Err("steinmetz_axes_must_be_distinct");
    }
    let middle = (0..3)
        .find(|axis| *axis != axes[0] && *axis != axes[1])
        .ok_or("steinmetz_invalid_axis_frame")?;
    // Use only proper cyclic permutations. The canonical X/Z cylinders are
    // interchangeable, so the requested order does not require a reflection.
    let world_axes = match middle {
        0 => [2, 0, 1],
        1 => [0, 1, 2],
        2 => [1, 2, 0],
        _ => return Err("steinmetz_invalid_axis_frame"),
    };
    Ok(Frame { world_axes })
}

fn to_local_bounds(bounds: [[i32; 2]; 3], frame: Frame) -> Result<[[i32; 2]; 3], &'static str> {
    let mut local = [[0_i32; 2]; 3];
    for axis in 0..3 {
        local[axis] = bounds[frame.world_axes[axis]];
    }
    Ok(local)
}

fn to_local_center(center: [i32; 3], frame: Frame) -> Result<[i32; 3], &'static str> {
    Ok(std::array::from_fn(|axis| center[frame.world_axes[axis]]))
}

fn to_world(block: BracketBlock, frame: Frame) -> Result<BracketBlock, &'static str> {
    let controls = block
        .controls
        .into_iter()
        .map(|p| {
            let mut world = [0_i64; 4];
            for local_axis in 0..3 {
                world[frame.world_axes[local_axis]] = p[local_axis];
            }
            world[3] = p[3];
            Ok(world)
        })
        .collect::<Result<Vec<_>, &'static str>>()?;
    Ok(BracketBlock {
        id: block.id,
        degrees: block.degrees,
        controls,
    })
}

/// Build equal-radius through-cylinder intersections along any two world axes.
/// The local core uses canonical X/Z cylinders; a proper cyclic permutation
/// maps its axes to the requested pair and their common normal.
pub fn intersecting_cylinders(
    bounds: [[i32; 2]; 3],
    center: [i32; 3],
    radius: i32,
    guard: i32,
    axes: [usize; 2],
    id: &str,
) -> Result<Vec<BracketBlock>, &'static str> {
    if id.is_empty() || radius <= 0 || guard <= radius {
        return Err("steinmetz_invalid_parameters");
    }
    for axis_bounds in bounds {
        if axis_bounds[0] >= axis_bounds[1] {
            return Err("steinmetz_invalid_box_bounds");
        }
    }
    let frame = frame(axes)?;
    let local_bounds = to_local_bounds(bounds, frame)?;
    let local_center = to_local_center(center, frame)?;
    steinmetz_box(local_bounds, local_center, radius, guard, id)?
        .into_iter()
        .map(|block| to_world(block, frame))
        .collect()
}
