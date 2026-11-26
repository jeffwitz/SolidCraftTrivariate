# Decision rule for the planar lookahead experiment

Evaluate the current automatic planner and the opt-in experimental strategy on identical source recipes. The frozen population is 45 public polygonal sketches and 12 previously published supplemental sketches. The external 11-sketch review set is a separate population; its original recipes were received before the prototype freeze. Only mixed_sizes was used for development. The other ten recipes were kept from the planner developer until the freeze.

Every admitted candidate needs the existing continuous local-regularity and exact section checks, conformity, equal physical volume and a successful native audit. Compare block counts before quality. At equal or lower block count, measure maximum in-plane spectral condition number, minimum corner sine, thickness-direction aspect and full 3D condition number separately on the published grids. Retain timeouts and failures in denominators. Timings are machine-dependent.

The experimental strategy is evaluated as an opt-in study. Keep the established default while waiting for independent out-of-corpus replay and a reviewed trade-off between in-plane quality, 3D quality and runtime. Report neutral and adverse results. The experiment makes no claim to solve concave or multi-level partitioning.
