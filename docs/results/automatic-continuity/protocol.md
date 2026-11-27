# Single-solve acceptance protocol

Baseline: b0b505319ed42f5674c2f7be7382453efe0d1fc3. The native field compiler and mechanical dispatcher select constraints from model data before stiffness assembly. Manual plans remain replayable. Compare against explicit choices in regression tests, outside the user solve. An edited geometry must lose any witness it no longer satisfies. Reject unsupported interface physics. Test classification, field dimension, the physical load and atomic CLI output.

The physical control uses the published curved blind pocket, one homogeneous material with E = 10 MPa and zero Poisson ratio, zero body force and a fixed lower face. Apply unit axial tension to the top and pocket floor. The exact solution is u = (0,0,z/10) mm. Source boundary moments give V = 48,976 mm³ and projected base area 2,500 mm². Expect applied force (0,0,2500) N, reaction (0,0,-2500) N and strain energy 2,448.8 N mm. Relative tolerances: 1e-8 on integrated quantities and 1e-9 on the reported residual. No mechanical solve may be used to select a field. Rebuild the English manual and check its diagrams and links.

## Full-field comparison

For the offline comparison use the same model, the same quadratic axial field with C¹ axial inheritance, cubic section fields, and a single section span. The reference plan is `{"schema":"inherited_field_plan_v1","plane_degree":3,"axial_degree":2,"axial_continuity":1,"restore_strip":false,"transverse_interfaces":[]}`. The automatic plan is `{"schema":"inherited_field_auto_plan_v1","axial_degree":2}`. Export each with `export-inherited-space`. Solve each saved package with `solve-inherited-space` and `{"schema":"section_affine_power_mms_v1","quadrature":10,"coefficients":[1,2,0],"offset":10,"scale":150,"degree":5}`.

| Lateral field | Free vector DOFs | Relative energy-norm error |
| --- | ---: | ---: |
| C⁰ on every lateral interface | 456 | 0.00523539912 |
| Automatic radial C¹ | 276 | 0.00821331628 |

The raw tension job uses the second space directly. The measured energy is 2,448.8000000 N mm and the reaction along z is −2,500.0000000 N. A successful physical patch test establishes reproduction on this load. The quintic comparison measures the loss of approximation capacity for a smoother field at this coarse resolution. Neither comparison is part of production field selection.
