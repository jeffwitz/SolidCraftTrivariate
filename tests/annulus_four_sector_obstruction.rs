//! Rank obstruction at a free edge shared by two smooth physical faces.
//!
//! For a regular 3D block map, let the free edge be the intersection of local
//! faces u=0 and v=0. If both physical faces have the same tangent plane along
//! that edge, then F_u, F_v, and the edge tangent F_w all lie in that plane.
//! The closed-cube Jacobian therefore has rank at most two. Knot insertion,
//! orientation changes, and nonsingular changes of parameter or world basis
//! preserve this rank loss. The example below is an explicit embedded map
//! whose Jacobian is strictly positive in the open cube but degenerates on the
//! free edge, so interior sampling alone cannot expose the obstruction.

type Vec3 = [f64; 3];
type Mat3 = [[f64; 3]; 3]; // row-major

fn dot(a: Vec3, b: Vec3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross(a: Vec3, b: Vec3) -> Vec3 {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn determinant(columns: [Vec3; 3]) -> f64 {
    dot(columns[0], cross(columns[1], columns[2]))
}

fn multiply(a: Mat3, b: Mat3) -> Mat3 {
    std::array::from_fn(|i| {
        std::array::from_fn(|j| (0..3).map(|k| a[i][k] * b[k][j]).sum())
    })
}

fn columns(a: Mat3) -> [Vec3; 3] {
    std::array::from_fn(|j| std::array::from_fn(|i| a[i][j]))
}

fn map_jacobian(u: f64, v: f64) -> Mat3 {
    // F(u,v,w) = (u-v, u^2+u*v, w). Its two boundary faces map to a line and a
    // parabola. They meet with a common tangent plane, as at a plane/fillet
    // junction, but have different second-order geometry.
    [[1.0, -1.0, 0.0], [2.0 * u + v, u, 0.0], [0.0, 0.0, 1.0]]
}

fn assert_close_zero(value: f64, scale: f64) {
    assert!(value.abs() <= 2e-14 * scale.max(1.0), "{value:.17e}");
}

#[test]
fn smooth_free_edge_forces_closed_cube_rank_loss_despite_positive_interior() {
    // The map is one-to-one: from x=u-v and y=u(u+v), u is the unique
    // nonnegative root of 2u^2-x*u-y=0, then v=u-x. Its exact determinant is
    // 3u+v, strictly positive in the open square and zero at the corner.
    for u in [0.05, 0.2, 0.5, 0.9] {
        for v in [0.05, 0.3, 0.8] {
            let j = map_jacobian(u, v);
            assert!((determinant(columns(j)) - (3.0 * u + v)).abs() < 1e-14);
            assert!(determinant(columns(j)) > 0.0);
        }
    }

    let j = map_jacobian(0.0, 0.0);
    let c = columns(j);
    assert_close_zero(determinant(c), 1.0);
    assert!(dot(cross(c[0], c[2]), cross(c[0], c[2])) > 0.5);
    assert_close_zero(dot(c[1], cross(c[0], c[2])), 1.0);

    // The two exterior face tangent bases coincide at the seam: each is
    // spanned by the same transverse vector and the same free-edge tangent.
    let face_u0_transverse = [1.0, 0.0, 0.0];
    let face_v0_transverse = [1.0, 0.0, 0.0];
    let edge_tangent = [0.0, 0.0, 1.0];
    assert_close_zero(determinant([
        face_u0_transverse,
        face_v0_transverse,
        edge_tangent,
    ]), 1.0);
    assert!(dot(cross(face_u0_transverse, edge_tangent), cross(face_v0_transverse, edge_tangent)) > 0.9);
}

#[test]
fn arbitrary_interior_knot_refinement_preserves_the_local_rank_obstruction() {
    // Knot insertion/refinement does not change the physical map. For any
    // locally regular reparameterization, its first span multiplies the two
    // transverse columns by nonzero one-sided derivatives; the edge column
    // remains in the common tangent plane. Exercise strongly nonuniform
    // interior knot locations and realistic positive local span scales.
    let cases: [[[f64; 3]; 3]; 5] = [
        [[0.01, 0.42, 0.99], [0.13, 0.6, 0.9], [0.7, 0.8, 0.95]],
        [[0.2, 0.2001, 0.8], [0.04, 0.5, 0.9], [0.17, 0.45, 0.7]],
        [[0.37, 0.61, 0.93], [0.25, 0.51, 0.88], [0.1, 0.3, 0.9]],
        [[0.00001, 0.5, 0.50001], [0.8, 0.9, 0.99], [0.3, 0.4, 0.5]],
        [[0.11, 0.13, 0.17], [0.9, 0.91, 0.92], [0.2, 0.4, 0.6]],
    ];
    let base = map_jacobian(0.0, 0.0);
    for knot_vectors in cases {
        for knots in knot_vectors {
            assert!(knots[0] > 0.0 && knots[0] < knots[1] && knots[1] < knots[2] && knots[2] < 1.0);
        }
        // Normalize each first knot span to [0,1]; these scales come from
        // the actual nonuniform knots rather than a test-only type marker.
        let first_span_scales = knot_vectors.map(|knots| knots[0]);
        let reparameterization = [
            [first_span_scales[0], 0.0, 0.0],
            [0.0, first_span_scales[1], 0.0],
            [0.0, 0.0, first_span_scales[2]],
        ];
        let transformed = columns(multiply(base, reparameterization));
        assert_close_zero(determinant(transformed), first_span_scales.iter().product());
        assert!(dot(cross(transformed[0], transformed[2]), cross(transformed[0], transformed[2])) > 0.0);
    }
}

#[test]
fn realistic_parameter_and_world_basis_changes_cannot_restore_rank_three() {
    let base = map_jacobian(0.0, 0.0);
    let q = std::f64::consts::FRAC_1_SQRT_2;
    let changes: [Mat3; 7] = [
        [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], // swap u,v
        [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], // reversed u
        [[1.0, 0.25, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], // shear
        [[2.0, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]], // rescale
        [[q, -q, 0.0], [q, q, 0.0], [0.0, 0.0, 1.0]], // rotate u,v
        [[1.0, 0.0, 0.2], [0.0, 1.0, -0.3], [0.0, 0.0, 1.0]], // mix edge axis
        [[q, 0.0, q], [0.0, 1.0, 0.0], [-q, 0.0, q]], // rotate about v
    ];
    for change in changes {
        let changed = columns(multiply(base, change));
        assert_close_zero(determinant(changed), 1.0);
        let minors = [
            cross(changed[0], changed[1]),
            cross(changed[0], changed[2]),
            cross(changed[1], changed[2]),
        ];
        assert!(minors.iter().map(|m| dot(*m, *m)).sum::<f64>() > 1e-6);

        // Rigid world-frame rotation and a nonsingular shear also preserve
        // rank. The world shear is intentionally not orthogonal.
        let world_changes = [
            [[q, -q, 0.0], [q, q, 0.0], [0.0, 0.0, 1.0]],
            [[1.0, 0.4, 0.0], [0.0, 1.0, 0.2], [0.0, 0.0, 1.0]],
        ];
        for world in world_changes {
            assert_close_zero(determinant(columns(multiply(world, multiply(base, change)))), 1.0);
        }
    }
}
