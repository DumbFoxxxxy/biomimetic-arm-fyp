# Findings

Source directory already contains 14 muscle OBJ meshes. Their headers identify
right-side BodyParts3D compatibility4.0, matching five labeled right-side bones.
All share raw mm coordinates; no independent alignment transform was needed for
the anatomical overlay. Old assembled STL offsets must not be mixed into this.

Closed edge check initially suggested35 triangular holes across11 muscles.
OBJ roundtrip disproved that interpretation: all35 boundaries are detached
single-face islands. Capping would create zero-thickness two-sided surfaces.
Cleanup instead removes these zero-thickness triangle islands, largest edge
0.356651mm, leaving the actual closed muscle volume unchanged. Original vertex
positions and source files retained. Validate cleanup through OBJ re-import,
not only an in-memory Manifold status. This does not infer fibers/attachments.

Available: subclavius; biceps long/short; triceps lateral/long/medial; brachialis;
coracobrachialis; infraspinatus; subscapularis; supraspinatus; teres major/minor;
levator scapulae. Deltoid, pectoralis major, trapezius, serratus anterior and
rhomboids are absent from these filenames/labels, not proved absent from all
possible user data. Need complement or explicitly design substitutes later.

Prior M2 ferrule STL and STEP export failed topology checks. Export process53773
was re-polled and confirmed terminal(exit1), not restarted after user correction.
That historical result is not a valid deliverable for this independent concept.

All14 real muscle volumes now have clipped BCC and BCC+surface-cage screens.
Uniform10x10x15mm cells / 2mm struts alone cause boundary fragmentation.
Adding a mesh-derived triangular surface cage restores each source volume's
component count. Three source inputs themselves contain tiny closed islands:
FJ1460 ~0.000206mm3; FJ1506 ~0.296928mm3; FJ1507 ~0.000366mm3. Retained explicitly.
PCA axes are grid coordinates, not muscle fibers. No channels or ends invented.

Legacy STL reader removes duplicate vertex-set triangles regardless of winding;
after float32 rounding this can falsely create boundary holes. Roundtrip tests
must also preserve all facets and their orientation, then validate closed solid
and relative volume. Do not treat legacy-reader NotManifold alone as proof of
invalid exported lattice, nor topology success as manufacturing approval.

M2 proximal-grounding adapters force a new GH-positive axial limit82.25deg,
not the former skeleton90deg. Remaining hard limit is extensor screw/scapula;
full-round flange had limited63.19deg.16 rigid insertion paths still produce183
TPU contact sample flags, not a fully proved assembly. Neutral100 physical pieces
vs78 known rigid moving and98 rigid neutral CAD are distinct scopes.
Coaxial8..14mm-radius rotors at +/-12..28mm axial stations all intersect original
M2 humerus. A real local rotor pocket can preserve one connected humerus and the
fork ears, pin and sleeves; candidate r10,+12 removes756.904mm3. This is geometry
only and does not supply a passive tendon route or muscle output constraint.

M4:detachable side lever avoids the integral rotor's insertion trap.600 tested
M3 straight approaches failed; this is not global path impossibility.8 finite
M4 install paths pass,12 standalone STEP and109-rigid combined CAD valid. Actual
new bone shape remains connected after42/58mm ulna spotfaces and brace seating.
Monotonic assumed axial muscle-end stroke can be achieved by changing real
eyelet phase, but this does not prove a route: selected15mm negative-side holes
62/102deg require32.345/42.229mm axial motion in that family; steel penetrates
hard parts112/64 samples, caps145/151 samples. Outlet-reference bends lack real
termination/guiding; never treat their polyline as an imposed physical guide.
Two passive stop barrels remain uncrimped and unretained:89 prescribed rigid
geometries means87 nominal bone/bolted pieces plus2 packing envelopes, not89
confirmed attachments.87/109/111/22 scopes must stay distinguished.

M5B:real ball joint cap terminals are necessary to avoid unsupported free cable
kinks. Nominally retaining a split lid0.3mm away is insufficient under tightening:
positive annular stops now set this gap; contact and further-closing interference
explicitly tested, while compression/preload remain unvalidated. Native Boolean
CAD removes MF internal zero-volume seam artifacts; winding numbers classified
deep rejected vertices as inside valid material, not exterior shape loss.
Source-neutral active terminal pull is[-.1726,.01418,-.9849] reversed for output;
direct old-lever tension assists it. M5B first output points along the opposite
direction and has a30mm neutral corridor, but this endpoint is NOT a guide.
Any actual fixed guide must be checked against moving forearm swept solids;
neutral corridor clearance does not prove installation or entire elbow motion.
M5A/M5B comparison, local25deg,20 fitting paths,131-neutral and129rigidCAD are
different scopes from full-arm drive. M5B whole installed motion/TPU/strength/
crimp still false. Root CURRENT_CANDIDATE index points to the current prototype.
Future source provenance must hash actual meshes: relativeXML strings can be
identical between sibling variants with different geometry. Animation inputs
are now captured before compilation and checked unchanged through rendering.

M6 exact guide placement:neutral clearance is weak evidence ofmotion clearance.
Of55 neutral-clear flexor first-sheave poses only16 clear the determinate150deg
forearm sweep;68extensor poses clearboth. Selected actualradius9.5mm flanged
wheel cases atsource distances50/30mm require10.787/15.476deg M5Bball articulation.
The first geometric route canthus be turned using real tangent/wheel contact,
but its keyed fork stillneeds an actual humeral groundframe. The finiteoutgoing
cable branch hasno endpoint support; neither branch is a valid imposed muscle
motion constraint. Local sleeve+wheelmust enterfork TOGETHER fromside: sleeve
R2.5 cannot pass theR1.7 cheek bore. Three actual installation sequenceschecked.
Axial clamp loadcan bypasswheel via4.8mm sleeve between cheeks while4.5mm wheel
hasnominal.15mm per-side clearance; thisis geometry only,notpreload/wear proof.

M7: groundedguide adds finite-volume support,so M6clear wheelposition doesnot
imply supports fit. Front verticalrails cross theside-lever sweptannulus even
whenkey receiverandfasteners themselves pass150deg. Tested36 front-axis-following
paths; noneclearfull elbow. Posterior route solves this: selectedflexXoffset20/
Yoffset+25,extXoffset20/Yoffset+35 relativehumeralaxis,lowerbrace roots,3mmrails.
Bothclear elbow150,GH+82.25/-90; full13suite verifiesno M4 range regression.
Frame fittingalsoneedsrealpath: flexnormal30mm approach fails,but40mm sideentry
at5mm normalstaging then5mm normalclosure passes. Original20mmfront mountshaft
couldnotinsertpastforearm. Rear18mm shaft andfrontloosenutfirst solveaccess;
20mmreverse tipintersectedforearm6.0177mm3. Nut access differsbetweenmodules.
151rigid CAD includesfreeend blanks,not151 load-connected pieces.109knownrigid
motion excludes44free/TPU instances. CAD keepsneutral freeballaxes while static
routing render articulatesballs10.787/15.476deg; bothareexplicitposesofsameparts.

M8 local receiver topology:real yaw shoulder sleeve plusorthogonal pitch sleeve
canaccommodatefull yaw andsampled±80deg pitch,notjusta loose steelbend through
old axialeyelet.13pieces each;20/25mm M3shaftsand6/11mmworkingstandoffclearlever.
Actual273combined conditionalposes clear,alljointsteps<=1deg,31extra bone refs
included; local3benchpaths each and13mountedpaths pass.27CADsourcevertexchecks
max~1e-10mm;compiled27MuJoCovertexerror7.8794e-6mm. Thisis NOT routedarmdrive.
Old62/102deg phasesgiveopposite moment signs.203.5099/323.5340 candidates keep
conditionalcentre-directed moment signs butareprovisional. Realrootdeparture
spanshitflangedwheel151/151 each(flexoutplane4.33..12.06deg;ext15.08..29.58deg).
Thusreceivinggimbals alonecannotrepairfixed-firstwheel plane incompatibility.
Prioritize directTPU-to-boneend interfaces andrealjoint-covering muscle path
comparisonbeforemoreexternalhardware. Existingpadchordtrajectory gives37.82%
flexshortening and23.66%extgrowth,withstraightchord moment signflip139/109deg.
Theseare endpoint geometry diagnostics,NOT achievedrange,strainorwirestroke.

PrintingKB:2026§4.4p15 TPU95A,FDMPrusaXL,PLAbones/supports;TableS2notinlocalKB.
Do notinvent nozzle/temp/layerheight. PaperYoungmoduli9.8±.7/7.4±.6MPa andbreak
elong650/480%notcurrent latticecalibration. EarlyTableI300Nbelongs tomotor-wire
unit,notcurrentmuscle output. M1cells15x20x10 differfrompaper15x10x30;nominal2mm
struts alone donotmakeperformance transferable. Userprinter/materialactualunset.


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
