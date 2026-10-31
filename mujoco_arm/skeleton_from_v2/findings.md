# Findings

V2 has seven hinge coordinates (SC2, AC1, GH3, elbow1). The DESIGN_BRIEF describes a mechanically different SC2/AC-ball background and a constrained three-motion shoulder model; current implementation must explicitly retain V2's chain rather than silently switching architectures.

Baseline exact visual meshes overlap at zero: scapula/thorax 18559 mm3; clavicle/scapula 3.4; scapula/humerus 748; humerus/ulna 201. Joint overlap requires actual solid inspection, not blanket exclusions.

RhinoCode list currently reports no running Rhino session. No running SolidWorks process was found. Native CAD backwrite is not yet verified.

Implementation experiment: smooth plate stand-off must first subdivide coarse planar faces. Without subdivision, displacement changes apparent thickness/volume substantially. Refined 4mm triangles, 8mm stand-off, retained AC/GH zones, R10 AC ball and local zero relief remove zero-pose intersections. This candidate is NOT acceptable: first collision K1~15.44deg, K2~3.13deg, K3~0.69deg, K4~3.13deg, K5~3deg. It demonstrates that zero-pose relief is insufficient, not final ROM.

Scapular fixed-normal swing/twist measure is used explicitly; a pure 50deg rotation cannot simply be assigned to one coordinate. Minimum-coordinate-norm solve reaches 50deg twist at [53.3514,-14.2484,-15.6676]deg in retained SC2/AC1 chain. Its geometry remains to be audited. This metric is a declared engineering measure, not a full clinical scapular-angle convention.

First combined swept-union dilation experiment was stopped for excessive compute cost; replaced with individually dilated cutters and incremental differences. SolidWorks installed at D:/SolidWorks2026/SOLIDWORKS/SLDWORKS.exe and a hidden session was authorized/started; initial COM version read failed with TYPE_E_ELEMENTNOTFOUND. Use typed interop rather than claiming CAD success.

Full shoulder chest subtraction rejected: ~89% scapular loss and three positive components. Retained-shape joint relief costs ~9% scapular, ~7% humeral material; larger disconnected machining chip 9.19mm3 explicitly removed, main joint carrier retained. SC fork corner-only trial rejected as a route to full shoulder travel because required cut enters protected bearing carriers. Better SC mount slot removes ~424mm3 from the chest, no intersection with protected bearing regions, and one connected chest solid.

Important corrected implementation error: unioning signed decomposed shells fills internal cavities. Preserve the original solid and subtract only positive isolated chips. Rebuilt joint cache with this correction; subsequent SC cache retains humeral screw cavity and hollow ulna. All old ROM reports preceding corrected cache are exploratory only.

Analytic chest proxy/vertex search is insufficient: exact pairs detect upper scapular contact even when proxy has +0.15mm margin. Added a growing analytic margin (up to2.15mm) and reran path search; found knots through28deg, next local optimization at29 failed. Not a global infeasibility result. Full mesh verification is authoritative.

OpenCascade7.9.3.1.1 and VTK9.6.2 installed into isolated tmp/skeleton_cad_vendor310 after sandbox network retry. Intended deliverable is explicitly faceted editable BRep STEP, not restoration of SolidWorks source feature history. SolidWorks active original document has unsaved changes; do not save, replace or close it. COM/ROT connection remains unavailable.

PowerShell default decoding garbles Chinese UTF-8 documents; use -Encoding UTF8. Root AGENTS.md absent. A guessed refine_geometry.py path was absent; actual helper is refine_gimbal_geometry.py (reference only).
