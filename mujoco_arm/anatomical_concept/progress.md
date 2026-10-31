# Progress

2026-10-30: Last user-facing reply only explained benefits of supplied anatomy;
it did not change state. Revalidated by inspecting actual muscle directory.
Found 14 existing labeled anatomical meshes; built independent raw-coordinate
19-mesh reference with no joint chain. Added inventory, source hashes, license,
edge topology diagnostics and explicit repair copies of triangular gaps.
MuJoCo reload passed (19meshes,19geoms,0joints). Rendered and inspected anterior/
posterior overlay. Continuing attachment/fiber and new joint-design work; full
goal remains active. Existing copies are references, not TPU/actuation proof.

Roundtrip audit caught an input-cleanup interpretation error: 35 triangular
boundaries were isolated zero-thickness faces, not holes. Updated cleanup to
remove these islands rather than add duplicate caps; preserve source vertices
and files. Initial fit audit correctly rejected the capped-copy roundtrip.

2026-10-30: Generated actual volume-clipped BCC structures for all14 source
muscles, then compared a surface-following triangular cage. All source islands
and lattice fragments retained, not silently removed. Surface-cage component
counts match14/14 source shapes. Material fraction20%-56% flags overly dense
thin regions for later anisotropic/surface-pattern refinement. Channel, endpoint,
joint and contraction tasks remain incomplete; no old chain is inherited.
Exported full-precision OBJ and STL independently; corrected a false-negative
legacy STL loader diagnostic by preserving oriented facets during roundtrip.

Verified all14 cage OBJ/STL roundtrips; maximum relative volume discrepancy
1.11e-7, maximum material outside original volume3.15e-13mm3 (numerical zero).
MuJoCo19-mesh/19-geom static model reload and full-anatomy renders passed.
Rendered and inspected FJ1486 source/BCC/cage detail comparison. Surface cage
restores connectivity but is dense/irregular, not a finalized compliant design.
Source/cage component counts14/14 match, including three source micro-islands.

User-requested consolidation: copied13 verified right-side supplemental OBJ
files from D:/XJTLU/isa_BP3D_4.0_obj_99 into muscle/. All13 copies match source
SHA256, original14 unchanged. muscle/ now contains27 OBJ files. Provenance in
muscle/supplement_sources.json. Prior inventory and lattice outputs describe
the original14 only; rebuild expanded inputs before making27-model claims.

Expanded inventory and reference scene to27 muscles+5 bones. All27 cleaned
copies reimport as valid volumes; added a separate rule for a detached4-face
open scrap in FJ1521(span~0.00672mm), preserving closed micro-volumes. Source
hashes unchanged;135 muscle/bone input pairs checked,45 intersecting pairs
reported without accepting them as legal engineered contact.

Built a sourced anatomical structure atlas for all27 muscles, exact nearest-
triangle centroid probes,46 NPZ candidate proximity regions and support-missing
records(72 total structure-role entries). Candidate regions are not confirmed
landmarks, drilling coordinates or fiber directions. FJ1512 insertion has no
2mm radius proximity (minimum36.446mm): common tendon/model segmentation or
an engineered tendon needs review. Trapezius whole-muscle stop sites are not
automatically imposed on every subdivision. Located source reference IDs for
sternum, ribs, costal cartilage, cervical/thoracic spine and occiput; no robotic
fixtures invented for missing soft-tissue supports.

Annotation viewer initially hid sitegroup3; enabled it and added ghosted anatomy
so green/orange region samples are actually visible. Static reference reload
passes32 geoms/32meshes. Joint constraints, steel channels and end interfaces
remain next work; goal not complete.

User confirms remaining bones in original BodyParts3D folder. Imported39
catalogued support references into 骨骼/支撑参考 (32bone refs,7costal cartilage),
no independent alignment, all source/copy SHA256 checks pass. Inventory/reference
now71 meshes:5limb+27muscle+39support; MuJoCo reload passed. Rerendered and inspected
whole overlay with skull/neck/thorax context. Updated atlas proximity extraction
to query imported support categories, leaving only soft-tissue supports undefined.

Independently fitted new joint candidates from current source and manubrium,
not old V2/A1 points. Clavicle/end-surface sample gaps SC2.475mm andAC4.468mm;
GH sphere radius21.482mm, radialRMS1.323mm. Elbow centre calibration shifts3.188mm
from initial distal-band estimate; 0..150 degree 1-degree single-direction
bone-only probe has no >.001mm3 intersection. This is NOT a whole-arm or hardware
acceptance. Built SC3/AC3/GH3/elbow1 mathematical joint candidate, nq13/nv10;
retention, full motion, controllability and compliant actuation remain open.

Continued joint validation:32 fixed bone envelopes expose original neutral
clavicle/manubrium overlap319.292mm3; original chain is blocked at zero. Built
independent S1 local SC ball/head+socket base+two split caps+eight nominal M3
fastener/washer solids, preserving source files and connected bone material.
Twelve positive components, no neutral pair penetration. Added actual hardware
to full-body moving/fixed audit, with42 fixed shapes and5 moving bones.

CAD failure investigation: manifold unions retained microscopic ancestry-boundary
faces and zero-volume shells, invalidating OpenCascade sewing. Resetting ancestry
with as_original(), controlled0.01mm simplification and zero-volume-shell cleanup
resolved STEP validity. Added proper washer counterseats after small neutral
washer/base intersections appeared. All12 STEP parts plus12-solid combined STEP
pass roundtrip. OBJ/STL checks and actual compiled-zero placement of16 geometries
pass; max mesh placement discrepancy~1.97e-5mm. No real force/strength acceptance.

Corrected global-axis GH labeling: old160-degree globalX probe is not anatomical
flexion. Reprobed about axes perpendicular to fitted humerus, used true rotation
norm for step≤1degree, not maximum vector component. Continuous nominal results:
flex154.5625, abd119.3125, elbow150, combined humerus164.049degrees. SC elevation
27.5degree isolated probe limited by cap, protraction1.6875 by scapula/rib,
retraction8.625 by clavicle/rib. Need coupled SC/AC paths rather than pair exclusion.
2026 forward/return sampled poses across13 probe paths; all limitations retained.

Six straight rigid-ball pullout probes block within0.2–0.8mm with caps held fixed;
not an arbitrary-path escape proof or a printed retention load test. Clavicle
ROI impact check found five prior region sets remain within1mm of modified bone.
Added audited-prefix kinematic animation and viewer keys1–5; static muscles hide
on motion so images do not imply actual TPU contraction. Goal remains active.

AC physical candidates: preserved SC baseline and generated large and compact
AC ball/socket/split-cap/fastener designs with unchanged anatomical centres.
All12 parts each are connected; neutral all-pair and humerus overlap checks clear.
Classified introduced abduction collision into added socket material, not source
scapular surface. Compact ball6/socketouter9 restores isolated119.3125deg and
combined humerus164.049deg vs large70.9375deg isolated; AC upward tradeoff
47.125->35.5deg documented rather than silently reducing a target.
Full13 forward/return paths include15 moving and42 fixed solids; AC accessories
correctly follow scapula. Rigid same-owner pair checks may reuse results by
transform invariance, never skip collisions. Scene registration26 geometries
max discrepancy1.97e-5mm; actual quaternion sample increment1deg within rounding.
Both12 individual STEP roundtrips pass. Limited six-direction pullout probes
block by0.9mm (large) or0.8mm (compact), not a general escape/strength proof.
ROI analysis identifies4muscles/5regions to remap in compact variant.
Viewer --AC --compact supports detail/full view and accepted motion;180-frame
GIF excludes static muscles. Current candidate remains unfrozen and undriven.

G1 GH physical stage: built shaped spherical bearing with two real neck slots,
two independently screwed cap sectors, retained ball centre and continuous
bone-side connectors. Initial dangling529mm3 cup patch was bridged rather than
discarded. Existing clavicle/AC cap neutral interferences were explicitly relieved;
saved rejected iteration in GH_connection/iterations/initial_neutral_conflict.
Newly created0.00406mm3 loose clearance sliver received recorded local0.721mm
radius contour trim, not hidden component omission. All12 solids connected.
25 moving/42 fixed solids,13 full forward+return prescribed paths checked with
1deg step and0.1deg boundary; flex145.9375,abd123.9375,coupled humerus164.994deg.
First flex/abd restrictions localized into added engineering geometry.36 compiled
zero-pose geometries align within1.97e-5mm; all12 STEP plus compound validate.
26 straight sphere translations blocked by0.9mm; load/arbitrary escape unknown.
All-pair minimum gap0 at same-owner scapula/ACcapright remains, so no manufacturing
clearance acceptance. Cuff ROI displacement8–12mm needs new muscle-side interfaces.
Compared20mm head/retained tuberosity with real external brace: more anatomy kept,
but movement blocked nearzero and humeral STEP Invalid BRep. Rejected comparison
is preserved separately, not passed by mesh/render validity. G1 review and180frame
accepted-prefix animation generated. Elbow physical joint and TPU/steel design
remain incomplete; full objective stays active.

E1 physical elbow: kept fitted centre/axis and added real fork, 6mm shoulder pin,
8mm OD bearing bush, locating sleeves, large washers/M6 nut. Source ulna/radius
are joined by two6mm bridges in one connected forearm; no pronation/supination.
Initial nine-part neutral overlaps fixed via real holes and end faces. First
moving version still could not assemble (144 straight candidate trials failed),
so saved geometry and made local39mm-radius side access. Sleeves extend through
fork holes and can be installed axially after seating forearm/bush. Checked
sequential insertion; final5mm forearm approach0.1mm, others1mm. Six forearm and
two pin translations block within2mm, not arbitrary retention/load proof.
Final nine parts connected and neutral clear at0.001mm3 threshold;31 moving
vs42 fixed solids,13 forward/return paths including elbow0..150 unchanged GH/SC
limits.42 compiled geometries align max2.21e-5mm. All9 STEP parts+compound pass;
assembled selected source layers into42-solid anatomical_arm_review.step, read
back valid (relative volume discrepancy3.93e-13). No duplicated superseded bone.
Produced128frame MuJoCo assembly GIF,180frame motion GIF and review/hash evidence.
New elbow-local ROI impact: brachialis33points/max6.45mm; medial triceps1/max1.60.
Shoulder cuff remapping, zero gap, real mount/frame and muscle drive still open.
Next deliverable is real shaped TPU antagonist modules with channels/ends and
measured required stroke, not more unconstrained anatomy overlay. Goal active.

M1 shaped module stage: read AIS2026 sections2.3/4.2 and PSM mechanism/figures.
Extracted17 actual sectional polygons perFJ1486/1477; concave brachialis sections
do not admit a simple ellipse/visibility kernel throughout. Cubic/periodic
B-spline envelope and surface diamonds; clipped eight-branch grid is not a full
Eq3 volumetric-coordinate reproduction. Active span180.34/205.94mm. Raw loft
outside anatomy fractions~1.55/3.29%, not an IoU or total retained-volume claim.
Generated3D lattice, guide walls/void, TPU pads, rigid reaction/termination caps,
4M3bolts+4nuts+8washers+crimp barrel:20 actual components each. Manufacturing
STL exceeds MuJoCo's STL-face limit; same-source exact OBJ used for viewing.
Initial per-section max-clearance poles introduced excessive wire turns; replaced
with fixed-end shortest-step interior route and tube axes bisecting straight spans.
Max nominal turns7.87/8.91deg, min nominal guide wall bounds1.045/1.031mm.
Filled only logged tiny generated closed pores and trimmed generated micro-slivers;
actual positive material remains single and source anatomy unchanged.
All40 OBJ/STL components pass; nominal0.6mm steel clears all parts, minimum~.1487mm
at pre-crimp barrel. Final-state cap/metal overlaps absent. Three nut insertion
paths failed despite clear final state; saved snapshot and added explicit lateral
access windows retaining3.8mm pad thickness. Both modules now pass11 sampled
cap/nut/bolt/barrel fitting paths. Actual threading/preload/crimp not simulated.
Four rigid cap STEP readbacks pass and cached cap inputs match current geometry.
Three rotated print candidates each now include source/target hashes and STL
roundtrip checks; no slicing, printer or physical deformation acceptance.
Straight anchor diagnostics: flexor51.89mm excursion,15/151collision poses;
extensor22.07mm,133/151collision poses, torque sign wrong from85deg. No fictitious
drum constraint used. Neutral raw modules collide rigidly10/3pairs with skeleton;
outward14mm anterior/6mm posterior first rigid-clear samples are review-only.
Updated source-vs-module MuJoCo figures, reviewed design and bench template with
unmeasured fields null. Physical bone adapters/output tendon guides, frame and
whole-arm muscle actuation remain incomplete; goal remains active.

M2 recovered completed local processes from current artifact evidence, no duplicate
jobs launched.27 new STEP parts valid; neutral100 has no significant hard/context
or TPU intersection, active nominal wires clear.16 fitting paths rigid-clear but
183 TPU contact flags remain.13 motion paths hash-match currentXML/assembly;
actual quaternion/hinge maximum1.000000000006deg,normerror3.33e-16.78 known rigid
motion instances,22 omitted compliant/free distal components. Positive GH82.25deg
first screw/scapula limit; elbow150deg, flex145.9375,abd123.9375,combined164.994.
Generated review_muscle_mounts.py and MOUNT_REVIEW.md/currenthash status.30 exact
translated original module rigid STEP parts added; assembled98 rigid neutral
pieces, checked STEP/OBJ vertices bidirectionally <0.001mm per instance, overall
98 solids valid and relative volume error1.44e-13. Two TPU bodies stay meshes.
ActualMuJoCo100-piece front/back static preview refreshed and visually inspected.
No physical actuation or whole-module movement claim. Source39 reference copies
rechecked SHA against user library. Two attempted multi-file documentation patches
failed atomically on absent context; corrected to matching file context.

Next M3 screen:40 real grooved coaxial rotor/collar geometries all collide with
current humeral elbow. Explicit0.35mm radial/axial pockets evaluated without
moving joint centre or changing fork ears/pin/sleeves. Least humeral material
removal among radius>=10mm eligible candidates is station+12,r10,756.904mm3.
Saved independent humerus_rotor_pocket and fused forearm_with_rotor; each1
material component. All new evidence remains separate from M2. CAD/neutral/elbow
diagnostic launched with known process handles; still requires inspected results.
Distal tendon guide, attachment/wrapping and actual muscle shortening not solved.

M3 experiment follow-through: direct sewing of union mesh failed OCP at0.001mm
and0.01mm simplification; sewing separate MF rotor also failed. Native24-sided
prismatic construction and BRep union with verified source forearm succeeded.
Native/mesh initial surface check found0.058516mm mismatch at groove shoulder;
did not relax threshold. Re-exported valid CAD as canonical manufacturing mesh,
positive material1, relativevolume correction3.20e-7. Rechecked2 exports,100/78
scene registrations(max2.12e-5mm),neutralallpairs and elbow151forward/151return
to150deg, all clear. Mesh hashes now captured BEFORE motion and verified after,
in addition to XML/assembly path hashes. CAD2 valid; sampled512 vertex-to-CAD
surface distances pass.16 muscle mount fitting paths still clear/183softflags.
However bare elbow assembly fails at28mm approach with fused rotor. Tried actual
0.35mm expanded rotor insertion recess; preserved no mutation on rejection:
humerus split into120549.16 and8.1201mm3 positive components. No fragment removed
to fake connectivity. M3 therefore NOT adopted, not frozen, not a drive pass.
M2 remains current98-rigid neutral CAD. M3 evidence/failed gate in independent
ROTOR_EXPERIMENT_REVIEW.md,candidate_status.json. Next compare mechanically
detachable torque-coupled rotor vs independently verified alternate installation
directions, before passive tendon guide/termination and muscle stroke validation.
Rebuild script records actual sequence; full shoulder13 suite NOT transferred to
M3. Current full objective remains active, not completed on these local checks.

M4 continuation: previous goal turn classified progress(CAD,geometry,failed
assembly evidence changed design direction). Revalidated current files/processes;
no stale running jobs assumed.600 M3 insertion directions/skews/flexions all
failed unchanged bones; documented finite search, not impossibility proof.
Implemented detachable side lever with source-ulna-section spots at42/58mm and
actual M3 hardware. Missing repairedFJ3391 path corrected to inventory source;
original raw OBJ stays untouched. Initial clearance cuts severed straight brace
from ears; archived and added8mm structural detour, then actual contour seats
to resolve86/133mm3 forearm overlap. All6 variants single-component parts and
0..150deg rigid elbow paths clear with all neighboring parts checked.
Added actual dual eyelets and passive crimp-stop envelopes; phase1deg sweep under
explicit unvalidated axial cap translation found monotonic stroke candidates.
Chosen review reference radius15,negative-axis56,phase62/102deg,32.345mm flex
shorten/42.229mm ext stretch. Not force/strain validated, no nonexistent guide
inserted into MuJoCo.12 STEP parts pass;111 neutral meshes/89 prescribed meshes
register max2.04e-5mm.89 includes2 uncrimped stopper packing envelopes, not87+2
grounded pieces.8 local install paths clear, no TPU flags for these specific
paths; M2's whole-module183softflags remain unresolved. Full13 motion paths
recomputed with actual mesh hashes; same SC/GH limits asM2, actual quaternion
step1deg/boundary0.1. Whole109-rigid STEP readback valid, sourceOBJ/STEP vertices
bidirectional<0.001mm, volumeerror1.55e-13. Actual MuJoCo close-up inspected.
Contact diagnosis rejects the assumed direct-tendon family: flex hits forearm
from39deg(112/151),ext hits humerus/spacer/pin/washer atzero(64/151), caps145/151
interfering from6deg. No direct-drive pass or fabricated printable actuation
claim. M4 rigid lever may remain as candidate; physical cap terminals/guides and
revised muscle movement are the next requirement. Added review/rebuild scripts,
SIDE_LEVER_REVIEW and source-hash status. Goal active, whole-arm muscles and
frame/TPU deformation/real actuation/printing still incomplete.

M5 continuation:previous goal turn is progress(M4 CAD/assembly plus failed
tendon diagnostics changed next action), no live process assumed from a file.
Built real split spherical terminals with source-cap M3 attachments retained:
11 rigid pieces each, ballR3.5/socketR3.7/outerR5.2, neckR1.2 and passive crimp
blankR1.6/boreR0.45. Real physical attachment replaces unsupported outlet bends.
M5A neutral131 has no hard/soft contact;24 local angular paths to25deg,28 finite
neutral-oriented pullout rays and20 combined fitting paths pass. Extensor rear
nut initially blocked; lateral nut-before-bolt order resolves it without fake
nut sockets.22 individual STEP plus129-solid neutral assembly pass. Dynamic
local GIF is prescribed mocap only; exact12 pose registrations checked.

Neutral-output projection showed M5A direct to old lever assists active pull
(both151/151), so generated independent M5B distal-load orientation. Both axes
now oppose neutral active pull(+1), first30mm cable corridor clear without
invented support. Lifted bridge root0.7mm to remove0.584mm3 TPU overlap, retained
failed candidates. Rebuilt stepped M2 bolts as a continuous mesh profile after
Boolean sewing failure. Added real0.3mm annular ear stops to carry nominal clamp
closure, because an open shell split alone could let tightening pinch the ball.
MF stop unions retained internal seam vertices and failed BRep. Diagnosed via
solid-angle winding:0.3mm-deep vertices are INSIDE valid CAD, material symmetric
difference~0.0044mm3, reverse exterior deviation~5e-6mm. Used native Boolean
shoulders on verified caps, re-exported canonical mesh, validated no exterior
bound relaxation, and reran actual geometry/assembly/retention/export checks.
M5B22 parts,131-neutral/129rigidCAD verified; all20 fitting paths pass with0TPU
flags, all25deg local sweeps clear, tested pullout obstructions0.2..0.3mm.
Printed stops have real lid contact and block an extra0.01mm closing probe.
129-solid STEP valid, sourceOBJ/STEP vertex registration<0.001mm. No materials,
preload/crimp/strength or full installed-arm motion claims. Four rigid housing
print placements are unsliced candidates. Bench template values remain null.

M5B direct-to-old-lever counterexample requires~141/108deg terminal direction
change and penetrates hard parts; this is rejected, not a reason to force25deg
ports toward invalid targets. Output first branch must connect via actual
turning/receiving geometry. Free cap/body deformation cannot be replaced by an
axial slider. Re-read paper pages7/8/10:empirical printed-sample mapping is
essential; reference compression is not this module's usable stroke. Shoulder/
scapular layout and frame remain open. Added shared current-candidate pointer,
CADs, actual local visualizations and rebuild/review scripts. Goal stays active.
An intermediate multi-file documentation patch failed atomically on stale
context; reapplied only verified contexts. Rendering API enum corrected using
installed MuJoCo enum inspection; shadows disabled visually to remove acne,
and camera direction checked from actual scene coordinates. Hidden anatomy in
detail renders is purely visual; full131 geometry is separately checked.

M6 continuation: previous turn's support-source check yielded39 hash matches,
but no new design state; resumed actual first-turn hardware rather than repeat
status. Recovered completed7-part cassette build fromactualdesign.json. Built
exact tangent/arc cable tube,checked361 wheel angles and3 side/axial local
installation paths. Seven faceted STEP exports valid/readbackpassed.
Initial source JSON read usedWindowsGBK default; fixed explicitUTF8. Raw FJ3178
support is notclosed; switched to existing documented repaired_supports copies,
without modifying originals or inventing repair duringchecking.131 neutral
manufactured parts plus31 additional bone references included. Initial first12
neutral-clear flexor positions only reached31.375deg; expanded to ALL55/68
neutral-clear placements before choosing.16flexor/68extensor guide positions
clear11determinate elbow solids through150deg at1deg; collisionboundary<=.1deg.
Shortest selected distances50/30mm,terminaltilt10.787/15.476deg; mutual two-guide
checks clear. Two finite outgoing40mm corridors remainunsupported. Rendered147
actual MuJoCo meshes; visually inspected elbow_detail.png. No fixed guide frame,
downstream transmission,free muscle motion or actual drive claim. RootM5B
candidate unchanged; new study hasseparateFIRST_GUIDE_REVIEW.md. Goal remains
active including whole shoulder/scapular muscles,frame,material/physical tests.

M7 in progress:made integral extensions ofthe two existing humeral adapters to
positive seated keyed fork receivers,plus8 real nominal mount-fastener solids.
First candidate10 CAD readbacks valid,neutral clear,elbow150 butGH axial severely
limited byframe/scapula/rib.16 changed assembly paths found2front-bolt approaches
blocked byforearm. Tested guide-first order:bareforearm radial entry thenhits
extensor guide/frame; rejected ratherthan omitthis finalconnection. A reused
sleeve metadata direction describedvelocity,notinitialoffset; correctedsign and
reran alternative-order check; sleevesclear, actualforearm obstruction remains.
Low-root secondcandidatefollows humeral axis,rear-entering18mm shafts (20mm shaft
tiphitforearm6.018mm3),frontloosenutsfirst.16changed fittingpaths nowclearwith0TPU
flags; extensor nut needsopposite lateral entry,2otherdirectionsretained asfails.
10CADparts valid. GHpositive limitrestored82.25,butelbow newly hitslever~18.8deg;
notadopted. Thirdseparate swept-root study nowgates actualnewsolids against
elbow0..150andGH+82.25/-90 BEFORE selection. Priorcandidates/sourcefiles preserved.

M7 selectedposterior branch:36front-rail candidate attempts didnotclear elbow;
isolatedreceiver/mount hardwareallpass,obstructionwasrail/lever interaction.
Archivedfrontsearch; posteriorflexcandidate1andextcandidate4pass elbow150 and
GHaxial+82.25/-90 before committing.16fittingpaths nowpass0TPU flags: flexbare
frame needs40mm sideapproach at5mmnormaloffset then5mmnormalclosure; frontnuts
needminusY flex/oppositeX extaccess.10CADSTLroundtrips pass,24changedCAD/source
vertexchecks pass,127unchangedM5sourceproofs reusedwithhashes,151-solid combined
STEP and24-solid interfaceSTEP valid afterreadback.109compiled mesh registration
max2.038289e-5mm.155meshesMuJoCo staticrender visuallyinspected; twofreeballs
tangent-aligned10.787/15.476deg versusnominalCAD pose explicitlydocumented and
recheckedwith2cablebranches against153parts(all4clear).Full13rigid path suite
completedwithnoregressionvsM4. ReviewinitialJSONwriter hitnumpy.bool_ conversion;
magnitude nowPythonfloat,rerunpassed andCURRENTindex updatedtoM7. A designbrief
patch failedatomically onextraneouscontext; reapplied verifiedsinglecontext.
Goal remainsactive:downstream receivingroute,muscledeformation/actuation,whole
shoulder/scapularlayout,motorreaction andphysical/material tests stillopen.

M8 active:13-piece yaw/pitch receiving terminal per muscle, actual shoulder
sleeves and blind-bore crimp blanks. Initial arms penetrated hub/sleeve and
washers; added realcentralfork clearance andtrimmed outercheekfaces. Local360deg
yaw at1deg and-80..80pitch at5 yawheadings pass,3bench installationpaths/module
pass;26localCADparts validreadback. Raisedreceiverworkingplanes6/11mm toclear
leverunderpitch,nominalM3lengths20/25mm. Rootoriginal phase probe showedwrong
moment directions;new203.5099/323.5340deg lever anchor candidates preserve bone
mounts and make provisionalcentre-directed signs correctthrough150deg. Missing
global bone-bolt holes initially filled bynewbraceunion; globally recut real
holes andhead recesses. Newlever/gimbal151pose probe clear; nowadaptivegrid
checks ALL receiver hinges <=1deg andcaptures sourcehashes toavoid stale data.
Current M7maincandidate remainsselected,no completed passive route yet.

User repliedprinter/nozzle/TPU allundecided;askedlocalKB references. Re-read
2026section4.4(page15),refs60/61 andearlyPSMTableI. Recorded TPU95A asstarting
hardnesscandidate,FDM/PrusaXL asliteratureprocessreference,PLA bones/supports.
No nozzle/layerheight claimed:localKB lacksTableS2. Addeddesign/TPU_PRINTING_
REFERENCE.md andmachine-readablecandidate;brand,printer,calibration null/false.
Paper650/480%failure elongations andmotor300N nottransferred tocurrent muscles.

M8 checkpoint delivered:27-piece receiving interfaceSTEP valid,sourcevertices
max9.97e-11mm;27STLroundtrips/27compiledMuJoCo max7.8794e-6mm.6benchand13mounted
paths clear(0TPUflags),includingextrear nutaxialentryafter5lateral optionsfail.
273conditionalposes(allhingesteps<=1deg) clearnewlever/receivers,extra31bone
refs andallactualcoupledyaw/pitch internalpairs. NativeunitCAD26current6/11mm
heights checked. New179meshMuJoCo detail visuallyinspected. Directactualroot
departureprobehitswheel151/151 both;extalsohitsfork73. Rejectedroute,noactuation.
M8review script passedscopedproofs whileM7current indexunchanged. Added
RECEIVER_REVIEW.md andhollowroll-swivel alternative blueprint. Also compared
directTPUdistalend connectionkinematics toavoidlocking goalto externalguides:
padchordflexshorten37.82%,extgrow23.66%,straightforcefirstsignreversal139/109deg.
No actualdirectadapter,softdeformation orrange acceptanceclaimed. Nextpriority
compare directbone-muscle/joint-covering design withrealturningguide option.

Errorscorrected:initialwrite importnotexportedby screenmodule,movedto actual
geometrywriter; forkarm intrusionremovedbyrealcentralclearance/outerface trim.
Newbraceunionrefilledboneholes,globallycutrealboltboresandheadrecesses. Adaptive
report variableangles accidentallyshadowedby module-rowdicts,renamedandreran;
sourcehash guards nowprotectstudies fromconcurrent/stalegeometry changes.
Onefindingspatch failedatomically onwrongcontext;appliedverifiedtailseparately.
Printingreference addedactualM1spacing15x20x10 vs paper15x10x30,axis-mapping caveat.
No full13M8armmotion/TPUphysics/physicaldrive claims;wholearm goalremainsactive.

## M9 direct TPU end attachment checkpoint — 2026-10-31

Created an independent direct-muscle candidate at
`exports/anatomical_concept/direct_muscle_interfaces/mechanical_yoke_shifted`.
The two distal caps are integrated into one physical forearm yoke; active wire
remains inside the shaped TPU. Original candidates and the M7 current index
are preserved. The flexor distal profile shifts by [13.177762,10.605833,-1.364878]
mm (16.97056 mm norm), blended through the distal profile; its actual connected
lattice, skin, channels and end pads were regenerated and checked. This is not
a render-only muscle. Local ulna seating removes 18.814045 mm3 within the two
existing mounting regions; joint centres are unchanged.

The 103 attached rigid parts retain the previous 13-path sampled motion ranges,
including elbow 0–150 degrees. Two TPU bodies and two uncrimped wire stops are
outside that rigid-motion acceptance scope. All 107 neutral manufactured
instances have no positive rigid interference or TPU intersection flags in the
neutral interface check; internal neutral threading also checked. The separate
whole-arm STEP has 105 valid rigid solids, volume roundtrip relative error
4.775e-14; it includes nominal uncrimped stops and excludes TPU. Faceted BRep
does not provide a native feature history or assembly mates.

Assembly-path report recovered: 17 paths, 15 rigid-clear; flexor proximal-left
and distal-right insertion paths have rigid obstructions. There are 109 soft
intersection observations, requiring deformation/physical handling assessment.
Complete assembly is explicitly NOT verified. Resolve insertion sequence or
geometry before selecting this candidate; do not present rigid ROM as muscle
drive, deformation or full assembly validation. No final freeze or main-candidate
replacement has occurred. Printer/nozzle remain unknown; the local-paper TPU95A
reference is already recorded separately in design/TPU_PRINTING_REFERENCE.md.


## M9 access relief / elbow-only scope — 2026-10-31

Previous goal turn yielded evidence and changed the authoritative progress record.
This turn: revised connected yoke with R3.5mm bolt-head entry relief, removing
328.777082mm3 while preserving the cap plate below localz6.15mm. Independent
yoke_access_relief_study assembly107; original shifted candidate preserved.
Changed order: flexor humeral front bolts/washers installed after cap bolts.
19/19 rigid fitting paths clear;109 soft observations remain, no full assembly
claim. Simple lateral staging and preloaded bolt alternatives failed and remain
recorded; no failing geometry overwritten.
Full13 forward/return rigid audit103 parts rerun, elbow150 and prior ranges
retained. STEP105 valid, source vertex maximum1.13835e-10mm, relative volume
4.46334e-13. Six generated STL pairs checked; MuJoCo neutral107/rigid103 source
registration passed; 2TPU+2stops versus32reference bones neutralclear. Revised
flexor standalone11 end fitting paths pass. Two bench-print sets and actual
yoke print candidate prepared, no slicing/printing qualification. Review script
passed. ELBOW_DESIGN.md and elbow_interfaces.json include real reaction shoulder,
steel path, crimp blanks and blank physical test record; mesh is not force proof.
Current candidate index now M9, with M7 historical index preserved.
User says finish elbow then stop. Do not start shoulder/whole-arm muscles.
Elbow drive remains incomplete: extensor contact/wrapping/force direction and
TPU deformation/stroke/friction require evidence. Do not pause merely because
rigid/candidate documents are finished; conditional stop is after elbow design,
not a false powered-elbow completion claim. Goal remains active within elbow.


## Extensor attachment alternative — 2026-10-31

Previous goal turn classified progress: M9 access relief,19 fitting paths andCAD
completed; elbow-only boundary recorded, compliant drive still unverified.
New finite search keeps bones/axis/near end fixed. Initial153 endpoint trials
(rotation phase andaxis shift) found no real-hardware150deg solution;108 met
2mmnegative-chord-moment screen butcollided. Expandedfamily also varies radial
position. Selected phase-44deg, axial0,radialscale1.2,delta[1.09951789179,
-13.1173089648,33.8943682421]mm,norm36.360705887mm. Actual distal hardware
clear151angles includingextra referencebones/crossflexcap; negative chord
moment magnitude>=2.113964948mm. Chord sign is onlya screen,not compliant
force proof. Original yoke cannot be reused after moving endpoint. Real
neutral extensor profile/lattice regeneration launched; main M9unchanged.
No shoulder/scapular muscle expansion authorized.


Real extensor regeneration completed:single connectedTPU volume29528.6742mm3,
21bridges,20manufacturedpart roundtrips andneutral wire clearance pass. New
profile span172.180469mm; proximal unchanged,distalshift36.36mm. Independent
mechanical_yoke_reanchored_pair created,selectedbrace reach10,dodge+18,extra0;
retainedactual flexor boltaccess void. All107neutral hard/softgeometryclear;
bothbone seats andinternal steelpaths pass.19fittingpaths clear109soft flags,
newextensor11standaloneendpaths pass. SixCADvalidandwhole105CADreadback done.
Endpointdemanddiagnostics showflexchord71.309mm/35.90%shorten;ext95.804mm/50.37%
grow. THESE ARE NOT TPU STRAIN OR ACTUAL WIRE STROKE. Straightthinextchord
firsthits actualhumerus132deg; flexchord nobonehits through150. Actual curved
TPU/wire/contact geometry still required,so do NOT promote candidate as full
150degdrive or mark elbowdesigncomplete. MainM9accessrelief stayscurrent.
New neutral107actualmeshes visuallyinspected;full103rigid13suite processlive
session84821 (lastobservedGHaxialpositive); newmeshregistration/printing checks
sessions95966/31781 must be polled to terminal. WholeCADsession5662 terminal0.
Goalactive; onlyelbowcontact/deformation remainsauthorized,not shoulder muscles.


Extensor-reanchor checkpoint:full103-rigid13path suite session84821 terminal0;
elbow150,other continuousprefixes unchanged. Whole105CAD terminal0,valid,
relativevolumeerror5.88176e-13. Meshregistration/STL session95966 terminal0;
printplacements31781 terminal0;19fittingpaths94624 terminal0. Allnewcandidate
outputs stable; no process remainsrunning fromthischeckpoint. MainM9access
candidate unchanged. NewREANCHOR_REVIEW.md documentsactual meshesandlimits.
NEXT:realbone-surface/compliantpathanalysis forstretchingextensor; thinchord
firstbonecontact132deg despitecorrectchordsign,requiredendpointgrowth50.37%.
No achievedmusclerange/TPUstrain/materialcapacity claimed. Elbowdesignnotyet
complete,keepgoalactive untilthatdesignresolved; stopafterelbowperuser.


## Bone corridor study — 2026-10-31

Previous turn was progress: real reanchored extensor and yoke generated, CAD105
and rigid103/full13,neutral107,19fitting andnew11bench paths all checked.
Maincandidate unchanged; corrected chord sign did not solve bone contact.
Implemented numerical free-node corridor study using actual cap-pad endpoints,
real humerus and coupledforearm solids; no guide mount added to physicalmodel.
Nearest-triangle distances with ray parity onlyaid optimization; final radius
2.1mmcapsule-solid intersection checked independently. Initial9angle sparse
optimizer found residual actualbone intersections135/145/150deg despite point
constraints. KeptfailedquickJSON,never classified muscleorwire aspassed.
Dense refinement launched session77028 with7samplesperlink,k128neartriangles,
0.3mmoptimizer margin. No actual lattice deformation, fullhardwareclearance,
friction, drive, or equilibrium proof follows from a clearbonecorridor.


## Real-port corridor and local swivel checkpoint — 2026-10-31

Densebone corridor9angles all exactcapsuleclear, butrealport-boundary steel
hardwarecheck failsat90+ (yoke/cap;nominalstop131+). Maxportturn127.58deg150.
Do notfreeze routingfrombone-only corridor. Allcorridorprocessesterminal0.
Created actualextensorcap splittrunnionrotor,localfork,2M2nominalshoulderpins,
2washers:6parts. Conditionalorientation151angles clearrotor+fork vsremaining
rigidarm;oldintegratedyokeabsent,notcompleteassemblyproof. Needrealbone-rooted
replacementbridge,no shoulderwork. Relativecapangle-.047..-142.555deg,axis
alignmentresidualmax8.308deg. Blindtrunnionboreavoidscentralwire. R1.8stub,
R1.1bore,ear15.2..18.2 vsstub10..15,nominalgap.2. Full6CADroundtripsvalid.
Failedfirstforkroundedwebbulgedaxiallyandhitrotor/pins/washers;real3mmaxial
slabweb fixesit. Firstunioncylinderpinmesh validbutCADinvalid;singlestepped
profilefix,all6CADrechecked. No looseningtolerances. Allprocessesterminal.
SWIVEL_REVIEW.md summarizes scopesandnextaction. ExistingM9mainunchanged.
NEXT:actualboneconnectionforfork; evaluateBOTHdistalcaprotationneeds,then
actualTPUdeformation/internalguidefriction. Passiveangle prescribedonlyfor
conditionalcollisionstudy,notforceequilibrium. Elbowstillnotcomplete;goal
active,stopafterelbowperuser. No printedperformance oractuatedrange claims.


## Paired passive caps grounded — 2026-10-31

Previous turn progress:actual6partswivel andhardbone/portdiagnostics; no blocker.
Flexor ordinaryforkfailsbone141 (front20),backforkhitsforearm atzero. Offset
familyselectsfront16+distal16mmbridge,conditional151clear. No anchors/axis
changes. Newpaired_swivel_support replaces oldintegratedfixedcaps/yoke with
actualconnectedframeon SAMEtwo ulnaseats,2rotors,4pins,4washers.117nominal
instances,113attachedrigid+2TPU+2uncrimpedstops. Commonframechosenreach10,
dodge+18,extra0 after3failedanteriorneutralframecases. Allchangedrigidpairs
againstparentarm/extrabones andtwo cap packssimultaneously clear266adaptive
conditionalstates,all3hingesteps<=1deg. Independentreturnrecompute passed.
All117nominalneutralpairshard0/TPUflags0. Actual117meshMuJoCo neutraldetail
visuallyinspected.11newCADparts valid.115wholeCADexport live session13432;
mesh/compiledregistrationprocesslaunched separately,check handle.
The passiveangles areprescribed,NOT equilibrium. NoTPU deformationor actual
steelguidepath validated; full assembly, bearingretention, physicaldrive still
pending. NominalCADcapsangles0; conditionalzero hascomputed smallcapangles.
Maincurrentindex staysoldM9pendingnewreview. UserstopafterELBOWonly remains.


## M10 grounded paired swivel checkpoint complete as geometry — 2026-10-31

CAD115 session13432 terminal0,valid,maxsourcevertexerror1.13835e-10mm,
relativevolumeerror2.22192e-13. New11partCADandSTL roundtrips pass; actual
117static/113rigid x4poses compiledMuJoCoregistration pass. Return266states
recomputed withallchanged rigidpairs/extrabones,neutral117hard0/TPUflags0.
review_paired_swivel_support.py passed sourcehashes/completecad/geometrygates.
3physicalhinges defined(elbow+2passivecaps),viewerstaticandprescribedrigid
90frameGIFgenerated. Neutral117 andrigid150 visualsinspected; noTPUanimation.
Allprocessesterminal,no runningjob. Maincandidate nowM10geometry;M9history
index preserved. DESIGN_BRIEF andCURRENT_CANDIDATE entry updated. Full
shoulderactuationnotstarted,norisitinremainingauthorizedscope.
NextELBOWonly:realpinretention andassemblypaths,internalsteelpathwithswivels,
actualguide/TPUdeformation/contact andusabletravel. Anglesprescribednotpassive
equilibrium; noforces,strength,friction,physicaldrive orvalidatedmuscleROM.
The userconditionalstopisafterelbowdesign;goalactivebecauseelbowunfinished.


## Captive pivot design — 2026-10-31

Previous goal turn progress:M10 paired groundframe,CAD115,113rigid,266paths
andindependentreturn complete. No activeprocess atstart,goalnotblocked.
Actualcustomhelicalcollar/caps trial capturedheadbutcollided withulna atzero
(frame7.48mm3,flexcover72.59mm3),so didnotselectordiscardbones. Preserved
captive_covers failure report. Alternative Ukeeper percap catchesbothgrooved
pinheads; throughM2clamp andnuts anchorittoforkbridge,notunmodeledblindthread.
Metal sheet candidate0.4mm,headgroove19.05..19.65R1.6,keeper19.15..19.55
innerR1.65. Four capture witnesses:rest/+0.09mmclear,+0.11mmblocked~.025mm3,
inward-.01mm washerblocked~.087mm3,minimumremainingengagement3.2mm.
Newcandidate127instances,15new/revisedparts;selectedbothclampoffsets+8mm.
Stale extensor normal initially shearedflexclampfastener;fixedbyusingmoduleF.
Local keeper relief needed,allchangesretainoneconnectedboneframe. Full266
changedrigidconditionalposespassedwithM2x10clampsforfullnutengagement.
Assembly:4pinbenchentry paths pass;nut beforekeeper fromnormal-plus passes,
trunnion-axisapproachesblocked. Frontclamps andwholeunittoBAREforearm n60
plus4boneboltpaths pass. Bareforearmsubassemblynotcompletehumerus/TPUassembly.
Keeper entrylast2.5/.4mm blocked onframe; retainedasfailure untilactualswept
reliefcorrected. TrueMinkowskisweep5mm with.01clearance launchedgen65884;
clampplanarfacesrestored,boresrecut,positiveframe/removedvolumegates remain.
Rerunlocalassemblycheckaftergenterminal;existingpathsreportwillbestaleuntilthen.
MainM10unchanged,no shoulderwork, no manufacturing/drive/TPUphysics claim.


## M10-R1 keeper checkpoint — 2026-10-31

TwoUkeepers+4groovedpinheads+throughclamps provide actualcapturegeometry.
All4rest/+0.09 clear,+0.11 blocked,reverse-.01 blockedbywasher;remaining
engagement3.2mm. Localplate sweptclearance removes29.710042mm3/frame stays
oneconnectedpart,clampplanarfacesrestored. ActualMinkowski5mmsweep+.01 margin
solvedfinal2.5/.4mm keeperentryobstructions. Backnut preinsertfromnormalplus
BEFOREkeepers;twotrunnion-axisentries percapfailretained.15selectedrigid
paths pass outof19trials;4rejected. IncludesunittoBAREforearm and4boneboltpaths,
notcompletearm/TPUinstallation. M2x10 clampscoverfullnutthickness.
CADfirstFJ1477keeperfailed:2nearzeroareatriangles(~8e-16mm2). ActualMFsource
simplify.001 removesdegenerates(292to268verts),volumechange~.007mm3;reran
266paths/assembly/CAD,notlooserCADtolerance. All15STEPvalid. Whole125CADvalid
relativevolumeerror1.60538e-12,maxvertexerror1.13835e-10mm.15STLroundtrips,
MuJoCo127static/123rigid x4poses pass. Return266 independentlyrechecked.
Sourceguardedreview_captive_keeper.py passed.127meshneutralvisuallyinspected,
rigidGIF90framesgenerated. WholeCAD70509/export47126 andallotherjobs terminal0.
Mainentry nowM10-R1 withM10historypreserved;DESIGN_BRIEF/KEEPER_REVIEW updated.
Thin0.4mmMETALkeeperforming/bendradius,actualfastenerdimensions,preload/loosening
andphysicalpulloutstrength remainunverified. Printed95ATPUcandidateunchanged.
NEXT strictlyelbow:internalwire/actualguidegeometryunder BOTH swivelcaps,
TPUlattice deformation/contact andusabletravel/passiveangles. Avoidfurther
rigidredesignunlessgeometry evidence requiresit. No actualdriverange claimed;
goalactive becauseelbowdesignstillunfinished. Stopafterelbowperuser.


## Real wire/guide study — 2026-10-31

Previous goalturn classifiedprogress:actualUkeeper geometry,capture,15selected
subassemblypaths andCAD125; alljobs terminal atstart. Currentstudyusesactual
M10-R1swivelcap portaxes fromcap pose matrices,notinherited onepiece-yoke
boundarydescription. Extracted19actual short-guide materialpieces perbody
fromfull~154kvertex/310ktriangleTPU meshes; originalthreadingsteelzerooverlap.
Free-centrelinequick16cases:extall8 clearbone/metal/guide;flex132/140 hitsactual
pairedframe andkeeper andsomeguidepieces.150clear doesNOT provecontinuous
range. Bone+frame/keeperconstraints variantstilloptimizerfails132,otherdata
retained. No fabricatedfixedspatialguides oractualmodel guideconstraints added.
Fullmanufacturingmesh diagnostic(noelasticity) underoriginalprescribedskinning
witharea compensation showslocalmappinginversions insomeposesincludingzero;
otherposes contactmetal. Thus thismapping isnot usable asTPUdeformationproof,
negative resultdoesNOT proveprintedmuscle cannotdeform. Do notanimate asreal.
Rest-shape-regularized/normal-bow initialization criticalflexor study launched
session79417. Actualsourcesunchanged;mainM10-R1notfrozen. Newregularizer is
numerical only,notmeasured stiffness/externalguideanchors. Baseline41030,
metalaware44112,fullmesh53161 terminal0. Scope ELBOWonly,stopafterdesign.


## Paper-method correction — 2026-10-31

Re-read actual2026section3.1p7–8:continuous wire-channel walls intentionally
stifferthan sparse lattice; parameters selectedEMPIRICALLY,n=3,not theoretical
prediction. Currentdiscreteguide model isengineering simplification, notthat
continuouswall. Do notspendfurtherfullscan asif itvalidatedthatstructure.
Own75798full532shortguidehypothesisanalysis intentionallyinterruptedviaits
actualtoolsession(ctrlC,terminalexit1) formodelreview,notobservationtimeout.
Partialcheckpointsavedstopped_for_continuous_channel_review_partial.json.
CIMprocesslookupdeniedclientaccess; notneeded,usedsameownsessioncontrol.
Regularizedfine16casesclearsteel/hard/guidebores,butindependentverifierfound
flexguide-self1pairat140/3pairsat150. No continuousphysicalrangeclaimed.
Addedconservativeguidebounding-sphere noninterpenetration onlytodiagnostic
fullscan,notactualTPUmodel. Fullmeshgraphprobe36159 terminal0:mathpatcherrors
~3.7e-12mm pass,but80iterationsnotconverged,materialvolume-18.69%,8metal
contactreviews+2.986mm3steelpenetration. Rejectasvalidatedmaterialresponse.
DoNOTconvertthis toprintedmusclefailureeither. Prior simpleskinningfailed
localorientation in5/6cases(includingzero),so no suchbodyanimation delivered.
Next alignswithpaper:actualcontinuousTPUwall1.6outer/.85inner R(.75nominal
wall) candidate,baseline shortguidesretainedasreinforcement; originalbodies
remainunchanged. Generator32463launched,needsverifiedterminal. Prepareactual
comparative lattice/guide/channelcoupons andexperimentalrecords,printerstill
unknown95A literaturecandidate. No newrigidhardware unlessevidencerequiresit.


## Elbow-only design handoff — 2026-10-31
Continuous-channel job completed: both bodies one material component, nominal wire intersection0; additions653.882/602.128mm3, neutral changed soft-pair flags0, rigid files unchanged. Prepared four actual connected empirical coupons, STL roundtrip pass, and12 blank specimen records; no physical measurements. Current127geom neutral MuJoCo compiled and rendered, image inspected. Indexed ELBOW_HANDOFF.md, actual interface coordinates and hashes; DESIGN_BRIEF updated. Prior hypotheses remain diagnostics, not validated TPU response. Actual driven range unknown, not frozen. Stop at elbow design-stage handoff per user; no shoulder continuation, no whole-arm/physical-performance completion claim.
