# Independent anatomical concept

## Current elbow-only handoff (2026-10-31)

User scope: finish elbow design, then stop; no shoulder-muscle continuation.
Current entry: `D:/XJTLU/PSP/mujoco_arm/exports/anatomical_concept/ELBOW_HANDOFF.md`.
Elbow candidate geometry, continuous channels, rigid CAD, neutral rendering,
empirical coupon meshes and blank measurements are delivered. This does not
complete physical actuation, nonlinear deformation, full assembly, or whole-arm work.
Earlier unchecked whole-arm tasks below are historical, not authorized next work.
Actual TPU-driven range remains unknown. No manufacturing freeze or physical pass.

User correction 2026-10-30: new concept, no required compatibility with V2/A1.
The whole musculoskeletal arm remains the objective, not a single coupon.

- [x] Identify real bone and muscle inputs, metadata, units and shared frame.
- [x] Preserve sources and prepare explicit closed copies of muscle volumes.
- [x] Render full anatomy and document available/missing labeled muscles.
- [x] Map all27 labeled muscles to sourced anatomical structures and extract review-only proximity regions.
- [ ] Label attachment regions/fiber paths; define new feasible joint constraints.
- [x] Build physical SC/AC/GH/elbow candidate connectors; check the elbow's nominal assembly sequence.
- [ ] Resolve remaining gaps/attachments and freeze manufacturable skeleton interfaces.
- [ ] Generate shaped TPU lattices with wire channels and real end interfaces.
- [ ] Assemble coordinated whole-arm geometry and prescribed motion paths.
- [ ] Export/reload CAD and report achievable geometry/drive limits separately.

Previous SC2/AC1 or K1-K7 is optional historical information, not a hard gate.
No actuator, joint architecture or anchor is frozen from old experiments.
Do not infer fiber directions from volume PCA alone. Do not infer anatomical
attachments solely from nearest bone distance. Extra model inputs can supplement
or replace this data without independently recentering their coordinates.

Checkpoint: both clipped-BCC and BCC+surface-cage screens cover all14 muscles.
All14 cage OBJ/STL roundtrips pass topology and relative-volume checks. Cage
component counts match source counts, retaining original closed micro-islands.
This does NOT complete shaped-muscle design: fiber direction, wire channel,
end transitions and compliant actuation remain open. Mesh-derived triangular
cage is exploratory and too dense in some thin regions; replace/refine with
controlled surface parameter patterns rather than freezing it on connectivity.

Current reference inventory now27 muscles+5 bones. The14-mesh lattice screen is
a limited earlier screen, not a27-muscle drive acceptance. Next: choose actual
engineering attachment regions, resolve trunk/neck support interfaces and
define a new feasible joint/kinematic chain; do not inherit old SC/AC limits.

S1 checkpoint: independently fitted SC/AC/GH/elbow centres and built a10DOF
mathematical candidate. Full torso-aware audit rejected original neutral
clavicle/manubrium overlap319.292mm3. Added a real SC local ball/socket/split-cap
candidate with12 solids including fasteners; zero pair overlaps and STEP/OBJ/STL
roundtrip checks pass. Remaining joints are still ideal constraints. Corrected
GH probe axes to be perpendicular to humerus; verified quaternion sample step.
Next: shoulder-girdle coupled paths and physical AC/GH/elbow interfaces, then
engineering muscle channels/endpoints and drive checks. No whole-arm completion.

AC checkpoint: two complete local ball/socket candidates now exist without
moving AC/GH centres. Large10mm ball introduced a socket/humerus restriction
at70.9375deg abduction; localized by exact intersection into added material.
Compact6mm ball/9mm outer socket restores119.3125deg and coordinated164.049deg,
while isolated AC upward probe reduces47.125->35.5deg. Prefer compact as next
work candidate, not frozen.12 connected solids each, OBJ/STL/STEP checks pass;
15 moving solids versus42 fixed references/accessories, no moving pair masked.
Four muscles/five source ROI sets need remapping after AC revisions.
Next: actual GH and elbow retention/connectors, feasible coordinated girdle
paths, then fiber-informed muscle geometry, steel channels/endpoints and drives.

G1 checkpoint: GH now has a real slotted spherical bearing candidate, headR16,
innerR16.4/outerR19.5, neckR4.5, two cap+M3x12 sets. All12 solids connected.
Resolved initial conflicts against existing clavicle/AC caps by explicit local
reliefs, saved initial rejected geometry.26 straight pullout directions block
within0.9mm; not general retention proof.36 scene registration checks and all12
individual+compound STEP roundtrips pass. Full13 paths: nominal flex145.9375,
abd123.9375, coupled humerus164.994deg; elbow150 still ideal, axial±90 clear.
First flex/abd limitations are in ADDED GH engineering material, not retained
original scapula. Significant shoulder-head change affects cuff insertion ROIs
by8–12mm: new bone-side muscle interfaces are REQUIRED. Zero minimum gap still
appears at scapula/AC cap right; no positive manufacturing clearance claim.
Larger anatomical-head comparison improves some ROIs but is blocked almost at
zero and fails humerus STEP. It is rejected, not substituted for verified G1.
Next: elbow physical hinge/retention and GH-side muscle interfaces; retain
whole-arm TPU/steel drive as the goal, not completion on a joint-only candidate.

E1 checkpoint: real clevis, shoulder pin, bearing bush, locating sleeves, washers
and nut now complement G1. Nine connected parts and actual ulna/radius bridges;
radius is included inside the merged forearm, not deleted to dodge collision.
Initial fit and then assembly failure were recorded in two iterations. Added
local side opening and through-sleeve bores: forearm+bush insertion, sequential
sleeve insertion and pin insertion now pass sampled geometry. Last5mm forearm
installation refined0.1mm. Six forearm translations/two pin directions blocked.
31 moving/42 fixed solids;13 forward+return paths pass expected limits, elbow150.
42 compiled geometry checks pass; nine STEP parts+compound and selected whole
42-solid skeleton STEP all reload valid. No native feature tree/mates, no frame.
E1 added brachialis origin changes33samples/max6.45mm and medial triceps1/max1.60.
Source shoulder zero gap and significant cuff remapping still not fixed.
Next priority: FJ1486/FJ1477-informed elbow antagonist muscle modules: actual
steel channel, hard end/interface, proximal reaction and tendon path; compute
required excursion and moment arm before claiming TPU contraction capability.
Then extend27-muscle layout and real frame/motor support. Do not stop at bone CAD.

M1 checkpoint: FJ1486/FJ1477 anatomy-derived variable-section muscle bodies now
exist, using17 extracted profiles, periodic/cubic B-spline surfaces, actual eight-
branch internal cells, controlled surface diamonds, short tubular wire guides,
TPU transition pads and hard end caps. Each20-component module has one connected
TPU material body. Four M3x12 bolts/four nuts/eight washers and real crimp-stop
envelope included; crimp/preload/threads still unvalidated. Source shell scraps
and generated closed pores are distinguished from positive material components.
Fixed-end minimum-step interior routing removes artificial pole jumps; nominal
guide wall lower bounds>1.03mm, wire0.6mm geometry minimum clearance~0.1487mm.
Nut insertion initially failed3paths total; explicit lateral access windows now
let both modules pass11 ideal fitting paths.40 exported parts roundtrip,4 cap
STEPs pass;3 print placements each are hash/roundtrip checked. Not yet printed.
Raw source straight-anchor probes fail: flexor excursion51.89mm/15colliding
poses; extensor excursion22.07mm/133colliding poses and wrong torque sign at85deg.
Ideal drums are only dimension comparisons; actual output guides remain to build.
Neutral module placement has10/3hard interference pairs; sampled outward14/6mm
clears rigid parts, NOT adopted or validated over motion. Next build real bone-
side mounting/reaction adapters and tendon/rotor interfaces, then integrate motion
with explicit TPU deformation assumptions/tests. Whole shoulder/scapular muscles,
frame and whole-arm active drive remain part of the objective, not completed.

M2 checkpoint: humeral mounting seats and two proximal grounding adapters exist.
Selected placements: flexor16mm anterior/extensor14mm posterior, centres fixed.
27 new CAD/STL roundtrips pass;100 neutral/78 known rigid meshes register within
2.22e-5mm.13 paths hash-match current XML/assembly; actual step<=1deg,boundary0.1.
PositiveGH axial82.25deg limited by extensor screw/scapula; full round flange
had limited63.19deg.16 rigid fitting paths clear;183 soft-contact sample flags
remain, full assembly unverified.22 TPU/free distal parts excluded from motion.
Next real distal transmission, grounded conduit/motor reaction, printed tests;
shoulder/scapula muscle layout, frame and whole-arm drive remain incomplete.

M2 integration:98 rigid neutral CAD pieces combined and read back, individual
sourceOBJ/STEP vertices matched bidirectionally<0.001mm, overall volume error
1.44e-13. Two flexible bodies remain separate mesh manufacturing candidates.
M3 actual rotor/pocket experiment passes neutral/elbow150/CAD2 geometry but FAILS
bare elbow insertion at28mm approach. Expanding that entrance disconnects8.12mm3
bone material, rejected without deleting fragment. Not adopted or frozen. CAD
microfacet issue fixed via native prismatic Boolean then canonical re-export;
all geometry checks rerun after0.0585mm correction. Current M3 is a documented
failed assembly experiment, not the selected arm. Next mechanical output route
needs a removable torque-coupled rotor or proven alternative assembly path.

M4 checkpoint:600 unchanged-geometry M3 straight-path candidates all obstructed.
Built six actual detachable forearm levers(radius10/15/20,axis+/-56), two ulna
spotfaces at42/58mm, two M3 bolt/nut/washer sets and two uncrimped cable stops.
Initial brace disconnected by bolt-head recess, then real detour restored one
piece. Local actual brace-seat cut removed overlap without breaking forearm.
Selected15mm negative-side reference:phase62/102deg holes.12 positive-connected
new solids,8 rigid install paths pass;109-neutral-rigid compoundSTEP/111 meshes
and all13 prescribed bone paths verified.89 prescribed meshes includes2 loose
stop positions, not fixed-drive proof.22 TPU/free distal omitted. Whole arm goal
remains open. Axial-cap diagnostic21..62mm excursions compared, chosen32.345mm
flex shortening/42.229mm ext stretching only hypothetical. Direct path fails:
flex wire112/151 samples,ext64/151; cap positions145/151. NOT a selected tendon
route or muscle strain target. Next physical distal terminal/guide and revised
cap/body motion, then measure printed contraction. Shoulder/scapular whole-arm
muscle layout and frame remain incomplete. Do not stop at a rigid lever/CAD gate.

M5B checkpoint:physical spherical output terminals(11 pieces per muscle) now
replace both distal caps. M5A direct-to-lever reference is preserved and rejected
for work-direction/routing reasons. M5B outlet axis opposes the neutral active
wire pull(projection+1); the first30mm free output corridor is clear, not a guide.
Added printed annular0.3mm clamp stops so tightening does not rely on an open
shell gap.24 radial0..25deg forward/return paths,28 neutral-oriented straight
pullout probes and10 local fitting paths per endpoint pass. These are geometric,
not strength/crimp/thread/TPU or full installed-arm proof. Native Boolean stop
CAD plus canonical mesh removed invalid internal MF seams without exterior
tolerance relaxation; all22 parts/129-rigid CAD/131 neutral meshes rechecked.
Four rigid-shell print placements and a null bench-data template prepared.
M5B is the current geometry candidate, not a frozen or actuated arm; unified
entry is exports/anatomical_concept/CURRENT_CANDIDATE.md. Next establish actual
turning/receiving structures or revise muscle/end layout, check swept bones and
free cap movement without invented rails, then compute/meaure usable stroke.
Full shoulder/scapular muscle layout, frame/motor reaction and real arm drive
remain in scope. Do not substitute local endpoint verification for this goal.

M6 checkpoint:actual seven-piece first-sheave cassette generated;7 valid STEP
readbacks,361 local wheel angles,one finite tangent/arc cable example and3
local fitting paths pass.80 positions per muscle checked against131 neutral
parts plus31 additional repaired bone references.55/68 neutral-clear positions
then ALL swept against11 determinate bolted forearm solids:16/68 clear to150deg.
Selected flexorL50,beta20,az135 andextensorL30,beta0,az0; paired guides/corridors
mutually clear.147 actual meshes rendered; review inexports/anatomical_concept/
elbow_guides/FIRST_GUIDE_REVIEW.md. Guides NOT grounded,40mm outgoing ends FREE,
TPU static. No drive,cap-motion,stroke or installed whole-arm motion pass.
Next build actual bone-ground mounting frames and receiving termination/routing,
then evaluate free cap/TPU deformation and empirically usable stroke. Need not
repeat completed M5 checks unless geometry changes. Source originals unchanged.

M7 checkpoint COMPLETE as rigid interface candidate only:chosen posterior swept
frames inelbow_guides/ground_frame_swept,notfront candidateorlow-root failedcase.
2integral adapter/frame revisions,8key fasteners,14existing M6guidepieces=24
interfaces;153neutral instances,109determinate rigid instances,44free/TPU omitted
withoutinvented motion constraints.16 changed fittingpaths andpositivekey seats
pass,includingreal side-then-normal adapter entry andloose nuts before rear bolts.
New solids elbow150 andGH+82.25/-90 gated before selection; full13rigid paths
rerun after selection,noregressionvsM4.151-solid neutralCAD readback valid,24new
freshvertexchecks,127unchanged verifiedsourceproofs reused.109 MuJoCo registered
max2.0383e-5mm. Fourneutral posed ball/first-branch checks against153pieces pass.
CURRENT_CANDIDATE nowM7; reviewGROUNDED_GUIDE_REVIEW.md includesactualranges and
failures. CAD nominalfreeballpose differsfromtangent-aligned render bydocumented
10.787/15.476deg articulation,not geometry. No frozen/drive/fullTPUassembly claim.
Next:actual forearm receiving termination anddownstream cable support,derive
usable compliant muscle motion without fictitiousrails; whole shoulder/scapular
muscles,motorreaction andempirical material/printing tests stillrequired. Asked
printer/nozzle/TPUhardness asynchronously; no assumption or dependency blockyet.

M8 checkpoint:two13-piece physical yaw/pitchreceiving units +rephasedlever=27
newrigid interfaces. Local360yaw/±80pitch sampled,6benchpaths and13changedmount
paths pass(0TPUflags);27CADsourcevertex checks and27STL/OBJ/compiledMuJoCo checks
pass(maxcompiled7.8794e-6mm). Rephased203.5099/323.5340deg candidates clear273
conditionalposes,allhingesteps<=1deg,including31additionalbone refs andactual
coupledyaw/pitchinternalpairs. Old62/102phasegeometry/moments rejected. Units
have6/11mmstand heights,M3x20/25 nominal. NONE ofthisproves fulltransmission:
actualfixedgroove-to-receiver departurehitswheel151/151 eachroute(flex4.33..12.06
deg outofplane,ext15.08..29.58). MainM7retained, M8study_status marksdrivefalse.
SeeRECEIVER_REVIEW.md andNEXT_ROLL_SWIVEL_GUIDE.md.

NEXT PRIORITY:compare direct bone-to-TPU-end connection andrealjoint-covering
muscle geometry beforeaddingmore externalcablehardware. Userdoesnotrequire
compatibilitywithourlocalprototypes. Endpointchordstudy derivedactualpadmotion:
flexchordmaxshortening37.82%,extmaxgrowth23.66%; straightchordforcesignfirstwrong
at139/109deg. NOT actualmusclestrain,wirestroke,clearance orachievedrange. Need
actualbone-sideadapters andcompliantpath/contact design; externalroll-swiveling
lastsheave isalternativeifneeded,notautomaticcurrentgoal. Preservewholearm
shoulder/scapulamusclelayout,motorreaction andTPUactive drive requirements.

Userprintingreply:allundecided,useKB. RecordedTPU95A ashardnesscandidate from
2026§4.4(page15) andPSMTableI; PrusaXL/FDM andPLA asliteratureprocessreferences.
TableS2notlocal:printer/nozzle/temperature/layerheight remainnull. ActualM1cells
15x20x10 differfrompaper15x10x30,keepperformanceunverified. See design/TPU_PRINTING_
REFERENCE.md and JSON. Do notaskuseragainforalready-statedunknown equipment.


## User scope change — 2026-10-31

Finish the elbow design and then stop, as explicitly requested by the user.
Do not start shoulder/scapular or whole-arm muscle design after this handoff.
Package the two actual lattice muscles, internal steel paths, real bone/cap
connections, editable CAD, assembly evidence, print candidates and test protocol.
Report TPU deformation, usable drive range, crimp grip and physical performance
as unverified unless evidence actually establishes them. This request does not
turn rigid elbow150deg clearance into an actuated-muscle pass.

Remaining: finish M9 access-relief CAD/motion checks, prepare shifted flexor
bench-print placements, provide consistent assembly/test instructions and review.
Stop after the elbow handoff; no completion claim for the original whole-arm goal.

Two attempted plan patches failed on a findings.md context; no file changed.
Recovered by appending this verified scope record directly.


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
