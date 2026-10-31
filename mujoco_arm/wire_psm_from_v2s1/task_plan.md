# V2-S1 steel/TPU muscle interfaces

Authorized task: determine steel paths and anchors on the frozen V2-S1 bones, preserving wire+TPU patterned muscles (not bare tendons as final actuators).

User confirmed: allow a TPU belly plus short passive steel tendon. A is fixed belly end; C is mobile belly end/active-wire termination, not a bone mounting point. User further confirmed: extending bellies is optional; do not force longer compartments or enlarge bones/brackets merely to hit length targets. Keep actual lengths and report coverage.

- [x] Reconcile V2 route units, bone-surface roots, muscle ends, wire terminations and guide centers.
- [x] Diagnose original routes on all accepted skeletal motion paths.
- [ ] Optimize mounting candidates and physical guide routes without changing frozen centers/chain.
- [x] Screen TPU muscle occupied volume and endpoint separation/stroke; retain uncalibrated limits as explicit scenarios.
- [x] Export candidate CAD-coordinate tables, viewable route model and verified/failed motion ranges. Native reference part128 point sketches saved/read back. Not a machining freeze.

Current selected review package: exports/wire_psm_from_v2s1/PSM-W1. Zero wire/bone, finite-guide/bone, wire-pair and belly/bone-envelope checks pass. Dynamic routing is NOT fully accepted; optimization remains open. Do not mark the overall anchor task complete or freeze holes based on this candidate.
