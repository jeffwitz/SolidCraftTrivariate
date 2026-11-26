# Exact cut propagation on filleted brackets

The frozen native kernel built ten of ten public controls. Six of ten passed with the reference executable. All four previously refused designs now have a native volume partition and an independently checked closed STEP solid. The general polygon benchmark uses a separate corpus. [Protocol and inputs](propagated-bores-protocol.md) · [Both executable fingerprints and per-case data](propagated-bores.json)

| Public control | Reference blocks | Candidate blocks | Reference max κ | Candidate max κ | Mapped closure rounds |
|---|---:|---:|---:|---:|---:|
| Original bracket | 10 | 10 | 13.33 | 13.33 | 0 |
| Two coplanar bores in one wing | 14 | 14 | 13.33 | 13.33 | 0 |
| East wing stack | 112 | 17 | 128.00 | 8.10 | 1 |
| North wing stack | 112 | 17 | 128.00 | 8.83 | 1 |
| Offset stack | Refused | 17 | — | 10.77 | 1 |
| Diagonal pair | 142 | 14 | 128.00 | 13.33 | 0 |
| Three bore elevations | Refused | 24 | — | 11.49 | 1 |
| A stack on each wing | Refused | 20 | — | 8.10 | 1 |
| Axial and transverse bores | Refused | 39 | — | 12.81 | 1 |
| Stack on unequal-width bracket | 112 | 17 | 192.00 | 7.41 | 1 |

Across the ten outputs, the independent oracle checked 253 exact oriented interfaces and 2,096 material probes. Every exported STEP file had one valid closed solid and one shell. The largest relative difference from analytic volume was 5.65e-16. Each map had a positive closed-domain determinant certificate. κ is the sampled spectral condition number on 125 points per map. There is no established block-count minimum for these 3D complexes.

Each east or north stack induces two requests on the hanging interface. Closing them gives three distinct cuts, including both bend maps and the opposite wing. On the diagonal pair the separator between the two tool centers stays within its wing and the final host faces already match. The mixed-axis case needs 18 cuts at closure. Removing every hanging junction from this placement therefore costs more blocks than on the parallel stacks.

This study covers orthogonal bores and axis-aligned rectangular hanging traces of the local hosts. The checks include all six admitted baseline models replayed with the frozen candidate. Source provenance and the exact cut parameter are available in each native model under `partition_plan.cut_sources`.
