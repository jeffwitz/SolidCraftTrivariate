//! Straight-tool volumes and curved wall controls provide independent checks
//! of material selection, embedding and boundary export.
use solidcraft_trivariate::{
    face_profile_extrusion::{FaceProfileExtrusion, build},
    interchange, step_boundary,
};

fn rectangle(operation: &str) -> FaceProfileExtrusion {
    let corners = [[3., 2., 1.], [7., 2., 1.], [7., 6., 1.], [3., 6., 1.]];
    FaceProfileExtrusion {
        spans: (0..4)
            .map(|i| vec![corners[i], corners[(i + 1) % 4]])
            .collect(),
        bounds: [[0., 10.], [0., 8.]],
        thickness: 6.,
        depth: 2.,
        operation: operation.into(),
        frame: [[0., 0., 6.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
        tolerance: 1e-6,
    }
}

#[test]
fn complete_box_and_tool_have_closed_exportable_shells() {
    for operation in ["tool", "cut", "add"] {
        let recipe = rectangle(operation);
        let model = build(&recipe).unwrap_or_else(|e| panic!("{operation}: {e}"));
        assert_eq!(
            interchange::audit(&model)
                .unwrap()
                .field("status")
                .unwrap()
                .text(),
            Ok("valid")
        );
        let step =
            step_boundary::export_step(&model).unwrap_or_else(|e| panic!("{operation}: {e}"));
        assert_eq!(step.matches("MANIFOLD_SOLID_BREP(").count(), 1);
    }
}
#[test]
fn open_profile_is_not_admitted() {
    let mut recipe = rectangle("cut");
    recipe.spans.pop();
    assert!(build(&recipe).is_err());
}

#[test]
fn quintic_profile_preserves_degree_in_the_exported_boundary() {
    let mut r = rectangle("tool");
    r.spans[0] = vec![
        [3., 2., 1.],
        [3.5, 1.9, 1.],
        [4.5, 1.7, 1.],
        [5.5, 1.65, 1.],
        [6.5, 1.95, 1.],
        [7., 2., 1.],
    ];
    let model = build(&r).unwrap();
    let step = step_boundary::export_step(&model).unwrap();
    assert!(step.contains("B_SPLINE_CURVE(5,"));
    assert!(step.contains("B_SPLINE_SURFACE(5,1,") || step.contains("B_SPLINE_SURFACE(1,5,"));
}
