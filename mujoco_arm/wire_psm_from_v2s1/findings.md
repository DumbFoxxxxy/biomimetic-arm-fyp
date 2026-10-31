# Findings

V2-S1 geometry baseline: 7 coordinates, 807 accepted forward poses. Original 16 wire channels are seeds; muscle count is not wire count. PSM fulltext confirms internal wire runs between muscle endpoints and TPU deforms with the path; curved/wide muscle shapes are allowed.

Critical metadata issue: route_metadata_v2.json labels local coordinates local_mm, but values are metres (matching XML), while SolidWorks_hardpoints_v2.csv is global design millimetres. Read XML as truth and explicitly convert.

TPU reference coupon is unfinished (no outer pattern or calibrated operating range). 25% cell-shortening rod assay is a screening case; no physical rating. Interface planning must separate fixed bone roots, offset muscle fixture eyes, internal wire path, guide centers and external motor feed.

Previous automation task interrupted by user: typed official C# ISldWorks probe succeeded, RevisionNumber34.0.0, zero open documents. No registry changes were made. Dynamic PowerShell dispatch still failed; typed bridge is the established working route.

New architecture selected after user confirmation: motor wire -> A -> TPU belly -> C termination; a separate fixed-length passive tendon joins C to guide and B. Initial free spherical guide caused bypass/nonmonotonic C coupling; rejected in favor of a finite-radius swivelling sheave with fixed incoming tangent G. Exact relation Lp=(AG-s)+D(G,B,q), so s=AG+D-Lp follows from physical reeving, not a general assumption equating arbitrary cable stroke to muscle shortening.

Implementation checks: preserve spatial tendon geom references when removing physical proxy geoms. Clear MuJoCo site_sameframe/geom_sameframe fast paths before mutating coordinates born at zero; otherwise updates are silently ignored. Sparse ten_J must not be reshaped as dense; use finite differences. These issues were corrected before accepting candidates. Old spherical-interface solver artifacts are not selected.

Swivelling sheave independent check:199 random samples, circle/tangent residual3.63e-15mm; <=2deg arc-polyline length error0.001224mm. Long-compartment candidate generated; user steering now makes length a bounded preference. Baseline V2 audit rerun with retained cylinder wraps:36 wire/bone failure groups across807 poses. No anchors frozen yet.
