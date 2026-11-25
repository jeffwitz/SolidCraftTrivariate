# J4d local fillet cost run closure

The J4d corpus has eight prescribed radius/path slots: four radii times
`from_h1` and `from_h2`. The frozen baseline summary revision is
`3d40313e3c376bb7ff2ecb3187b9959c16ea17b7`; the original 76 baseline trials
remain included. The detailed measurements and raw-output references are in
[`local-fillet-cost-20261125.json`](local-fillet-cost-20261125.json).

All four `from_h1` ladders completed through five states. Each first qualified
at step 3, 624 leaves; q8 measurement times were 3.53–3.65 s and cumulative
search times were 17.92–18.16 s. The final 736-leaf state also qualified for
all four radii. Per-trial GNU time records show 25.71–26.16 s and peak RSS
116,952–117,324 KiB for those runner processes.

For `from_h2`, the retained historical run report records timeouts at 1800.10 s
for `[40,80]`, 1800.12 s for `[80,40]`, and 1800.10 s for `[50,74]`. The
`[50,74]` raw JSON contains only step 0 (2,304 leaves, 953.81 s, q8 solve
791.31 s, 9,135 iterations); that partial state did not qualify. The GNU time
resource files for these attempts are empty, so no peak RSS value is claimed.

The `[51,74] from_h2` slot must not be classified as unstarted: PID 2117636 was
directly observed running `local_fillet_cost --path 51 74 from_h2` at about
20m27s elapsed. No final output or resource line remains. The report therefore
records `launched_termination_unconfirmed`; final completion, termination
time, and final RSS are unknown. A final process check found no active J4d
process, but cannot recover the missing terminal record. Thus the campaign
has four complete h1 ladders, three measured h2 timeout attempts (one with a
partial first state), and one launched h2 run with unconfirmed termination.

The report was reconstructed from retained JSON and resource files without
rerunning the expensive paths. Rebuild commands are recorded in the JSON.
