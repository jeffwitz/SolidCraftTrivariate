# SC8 protocol: blind bores through G21

## Revision and scope

Starting revision: `25676c258521462814353fc72e9e10228cb1f634`.
This adapter extends the direct SolidCraft command-history importer with one
explicitly qualified blind cylindrical cut, represented by G21's exact axial
interval. Existing SC2, SC3, SC5, SC6 and SC7 inputs must continue to route to
their existing adapters byte-for-byte.

## Frozen host operation

Use the real SolidCraft CLI (`0.48.1`) to execute a two-command history:
`solid.box` creates an axis-aligned rectangular base and `solid.hole` cuts a
simple, flat-bottomed cylindrical hole from its top face. The first witness is
a 40×40×12 mm box at `[-20,-20,0]`, with an 8 mm diameter hole centered at
`[0,0,12]` and depth 4 mm. Its exact G21 interpretation is a CCW rectangle,
height 12, zero top shift, center `[0,0]`, radius 4, a disjoint guard radius 8,
and material-removal interval `[8,12]`. The through-all default, drilled tip,
countersink, counterbore, threaded/modelled variants, and any non-axis-aligned
or non-lattice input are outside this SC8 subset and must be refused without
changing the prior model.

The host must save and reopen the result as `.sc` and STEP. Independent OCCT
checks compare one connected solid, bounds, volume, center of mass, and
point-classification witnesses. The blind floor is the plane `z=8`; points
under the bore remain material below that plane, are void immediately above
the floor, and the circular wall extends only from `z=8` to `z=12`. The analytic
volume is `40*40*12 - pi*4^2*4 = 19200 - 64*pi mm^3`. G21's independent audit
must report the upper cap at z=8 as `bore/<id>/upper`, the outer top at z=12,
and no bore-wall selection below z=8.

## Additional coverage and refusals

Before changing production code, preserve the host-generated history and
reopenable `.sc`/STEP artifacts. Add one lower-face blind-hole witness if the
host emits it as a supported direct history; map it to `[0,depth]` only after
checking the host's direction semantics. Add a decimal-zero case only if a
host-admitted saved history actually contains a `0.0` coordinate; otherwise
record that the parser accepts exact zero spelling but do not claim SC8
coverage for it. Reject malformed interval/depth, nonintegral coordinates,
unsupported hole types, failed host replay, and cap/tag or volume mismatches
atomically. Report the exact host version, native/G21 audit results, STEP and
OCCT measurements, classification witnesses, old-subset replay, and limits.
