# SC7: import extrusions on origin XZ and YZ planes

## Starting point and files

- Starting revision: `8aaf73c8a0c9cb58141bed1436c2c38c341c5a99` (SC5/SC6 adapters
  are the current uncommitted worktree inputs).
- Native mapping: use the existing SC5 polygon/G17 or SC6 circular-polygon/G20
  recipe as the source model, then represent the world placement by a G19
  `affine_image_v1` transform with determinant +1, denominator 1 and zero
  translation. Do not rotate coefficients by an independent code path.
- New code is isolated to a SolidCraft plane adapter, its tests and reference
  corpus. CLI dispatch is coordinated with the root integrator.

## Frozen host subset

Accept only `solidcraft/66` version `0.48.1` timelines already admitted by SC5
or SC6, with a sketch on the origin XZ or YZ plane and no plane flip, offset,
attachment, parameterized placement or additional feature. Keep the original
SC5/SC6 requirements for geometry, profile selection, extrusion extent and
integral dimensions. XY imports continue to use their existing routes.

The SolidCraft command catalogue states that the sketch basis is `(+X,-Z)` on
XZ and `(+Y,+Z)` on YZ. Host-generated asymmetric STEP geometry confirmed those
bases and positive extrusion directions: XZ occupies x `[10,22]`, y `[0,6]`,
z `[-27,-20]`; YZ occupies x `[0,6]`, y `[10,22]`, z `[20,27]` for a local
profile bounding box `[10,22]×[20,27]` extruded by 6 mm. Both contain one valid
solid of volume 360 mm³ for the initial probe. The frozen qualification profile
is widened asymmetrically to leave room for a certified bore guard. Positive
extrusion maps are therefore

| Plane | Local sketch coordinates and extrusion direction | G19 matrix `A` |
| --- | --- | --- |
| XZ | `(u,v,w) -> (u,w,-v)`; positive distance along world `+Y` | `[[1,0,0],[0,0,1],[0,-1,0]]` |
| YZ | `(u,v,w) -> (w,u,v)`; positive distance along world `+X` | `[[0,0,1],[1,0,0],[0,1,0]]` |

Both matrices have determinant +1, preserving the certified source
orientation. A host result inconsistent with either basis or direction blocks
the mapping until the convention is corrected and documented.

## Qualification and refusal rules

Freeze real host `.sc` documents and extracted histories for an asymmetric
polygon-only extrusion and a polygon with one circular hole on each plane.
Use unequal side lengths and non-centred bore locations so axis permutations
and sign errors are observable. Require the G20 certificate for bore cases.
The independent FreeCAD STEP oracle checks one valid solid, bounding box,
analytic transformed volume, surface area and centroid, symmetric-difference
volume below `1e-8 mm^3`, volume error below `1e-8 mm^3`, area error below
`1e-8 mm^2`, centroid error below `1e-9 mm`, and matching material/bore/exterior
point classifications away from boundaries.

Include two atomic refusals made from real host histories: (1) an offset custom
XZ plane; (2) negative distance. The host accepts negative distance and records
it as `direction:"positive", distance:"-6"`; the importer refuses it because
SC5/SC6 certify only positive extrusion distances. The refusal must leave an
existing native output byte-for-byte unchanged. Preserve the exact refusal
reason and history. Re-run SC2 and SC5 dispatch/oracles to ensure the added
plane routing does not change existing XY imports.

## Evidence and limits

The host omits a zero constraint counter (`k`) from some sketch histories. The
adapter adds `k:0` only to its private XY-normalized clone before passing it to
the existing SC5/SC6 decoder; this is a lossless spelling normalization.
Host STEP placement establishes the plane basis and positive extrusion sense;
G19 supplies the exact coefficient transform and source-certificate
revalidation. FreeCAD comparisons are differential evidence on the frozen
cases. They do not extend SC5/SC6 beyond their recognized XY timeline grammar
or establish a general attached-plane implementation.
