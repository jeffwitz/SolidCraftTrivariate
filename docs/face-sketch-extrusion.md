# Face sketch extrusion

Baseline: main 608d88be. The compact ribbon previously accepted an origin-XY polygon and created a separate body. The new entry resolves the sketch, its attached planar face and its owning body before editing the document. The original sketch and supporting feature must survive undo and save/reopen.

Acceptance uses real SolidCraft face attachments and both fit-point and control-point spline fixtures. Run addition outward and removal inward from identical inputs. Check the number of closed bodies and compare material change against an independent SolidCraft sketch extrusion. Rejected operations must leave the source intact. Geometry and topology, stock recognition and host boundary composition are reported separately.

## Geometry contract

The ribbon reads the authored contours in sketch coordinates, skipping construction curves and linked face outlines. The supporting face and target body come from the recorded attachment. A planar face normal from the host resolves the material side, including when the sketch axes point inward. A positive depth always means outward for addition and inward for removal.

The importer extracts homogeneous Bézier spans from control B-splines by knot insertion. Fit-point cubics use chord-length interpolation with the periodic end conditions for a closed fit. A weighted spline retains its weights. Input controls from the host use binary64; the native model stores integer homogeneous controls reconstructed within an explicit positional budget of one nanometre (0.000001 mm). Half-span subdivision in the charts uses integer de Casteljau arithmetic.

The present partition admits one closed loop. Its endpoint polygon and curved span control hulls must admit the radial partition checked by the constructor. A scaled inner polygon separates its central quadrilaterals from the curved walls. For rectangular stock, the builder adds an outer ring to the stock boundary and selects material by depth on matching layers. Every emitted volume undergoes a Bernstein determinant test. Higher-degree swept geometry has its own admission path; mechanical analysis of degree-five fields is outside this qualification.

* **Full result.** When the outer partition of a six-face rectangular stock aligned with the sketch axes passes certification, the selected pocket or boss is exported from the native volume model. The host replaces the stock at the end of its timeline, preserving the earlier supporting feature and sketch. The JSON filename ends in `.result.model.json`.
* **Extrusion tool.** On other admitted planar faces, or when the stock partition fails but the tool partition passes, the kernel builds the complete curved tool. SolidCraft applies its boundary to the specific owner with a join or cut. The JSON filename ends in `.tool.model.json`. This represents only the tool; the final complex stock has a host B-rep. The plugin returns the stock partition diagnostic when applicable. The result must be a closed body with a material change in the requested direction.

The audit reconstructs the volume maps from the stored source recipe. One plugin transaction covers all document edits and companion-file writes. To compare another depth or edit the contour, undo the extrusion before rerunning it. The imported STEP boundary records the shape at execution time.
