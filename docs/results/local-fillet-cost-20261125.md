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

The `[51,74] from_h2` slot was directly observed running as PID 2117636 at
about 20m27s elapsed. Its later retained wrapper record reports
`timeout_1800s` at 1800.10 s, resolving the earlier in-progress observation.
Its GNU time resource file is empty, so final RSS is unavailable. Thus all
four h2 paths have timeout records; the `[50,74]` timeout retains one
nonqualified first state, and the other three retain no final JSON output.
A final process check found no active J4d process.

The report was reconstructed from retained JSON and resource files without
rerunning the expensive paths. Rebuild commands are recorded in the JSON.

The terminal record for the final `[51,74] from_h2` timeout is preserved in
[`j4d_terminal_records.json`](../../tests/references/j4d_terminal_records.json).
It ties the run to the 8d197a8 J4d manifest, 51/74 sixteenth-mm radii, the
`from_h2` depth-two / three-round path, target roots and chart coordinates,
and the 1,800-second Python wrapper. The runner source matches revision
d7d764a. The release binary was built at 07:28 local, before that commit at
07:50; its review-time SHA-256 is recorded, but no build metadata identifies
the exact source revision, so the report leaves that field unknown. The
preserved timeout diagnostic and elapsed time prove termination by the
wrapper; the empty GNU time file supplies no RSS result.
