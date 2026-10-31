"""A2: shoulder-clavicle ball/socket attached to actual clavicle/scapula.
Preserves SC candidate; dynamically owned AC caps are not treated as world-fixed.
"""
from build_SC_connection import unit,rod,write
from repair_inputs import *
from extract_attachment_regions import projector
import os
sys.path.insert(0,str(ROOT/'mujoco_arm/blade_motion_design'))
from geometry_distance import segment_mesh_distance
BASE=OUT/'joint_model/SC_connection';DEST=BASE/'AC_connection'
COMPACT='--compact' in sys.argv
if COMPACT:DEST=BASE/'AC_compact'

def load(path,raw=False):
    v,f=read_obj(str(path))
    if raw:v,f=weld_stl(v[f]);v=v-SHIFT
    s=mf.Manifold(mf.Mesh64(v.astype(np.float64),f.astype(np.uint64)))
    if s.status()!=mf.Error.NoError:raise RuntimeError('Invalid source '+str(path))
    return s,v,f

def run():
    DEST.mkdir(exist_ok=True);definition=json.loads((OUT/'joint_model/definition.json').read_text(encoding='utf8'));inv=json.loads((OUT/'inventory.json').read_text(encoding='utf8'));byid={r['id']:r for r in inv['models']};C=np.array(definition['centres_source_mm']['AC'])-SHIFT;S=np.array(definition['centres_source_mm']['SC'])-SHIFT;H=np.array(definition['centres_source_mm']['GH'])-SHIFT
    oldclav,cv,cf=load(BASE/'clavicle_SC_ball.obj');oldscap,sv,sf=load(ROOT/byid['FJ3384']['source'],True);hum,hv,hf=load(ROOT/byid['FJ3368']['source'],True)
    hint=unit(S-C);possible=cv[np.linalg.norm(cv-C,axis=1)>23];root=possible[np.argmin(np.linalg.norm(possible-(C+hint*25),axis=1))];z=unit(root-C);x=unit(np.cross(z,H-C));y=np.cross(z,x);F=np.stack([x,y,z],1)
    def local(s):return s.transform(np.c_[F,C])
    def box(lo,hi):return local(mf.Manifold.cube(np.array(hi)-lo).translate(lo))
    ball_radius=6. if COMPACT else 10.;inner=ball_radius+.35;outer=9. if COMPACT else 13.;neck_radius=3. if COMPACT else 4.;capdepth=3. if COMPACT else 4.;relief=16. if COMPACT else 22.
    clav=(oldclav-mf.Manifold.sphere(relief,32).translate(C))+mf.Manifold.sphere(ball_radius,32).translate(C)+rod(C,root,neck_radius)
    shell=mf.Manifold.sphere(outer,32).translate(C)-mf.Manifold.sphere(inner,32).translate(C);rear=shell^box([-30,-30,-30],[30,30,-.1]);caps=shell^box([-30,-30,.1],[30,30,capdepth]);base=oldscap-mf.Manifold.sphere(12 if COMPACT else 18,32).translate(C)
    candidates=sv[(np.linalg.norm(sv-C,axis=1)>19)&(np.linalg.norm(sv-C,axis=1)<45)];_,hd,_=projector(hv,hf)(candidates);roots=[]
    for sign in [-1,1]:
        tip=C-z*(outer-1.5)+x*sign*2;target=C-z*25+x*sign*4;score=np.linalg.norm(candidates-target,axis=1)+100*np.maximum(0,6-hd)**2;chosen=None
        for k in np.argsort(score)[:80]:
            if segment_mesh_distance(tip,candidates[k],hv[hf])>=3.3:chosen=candidates[k];break
        if chosen is None:raise RuntimeError('No GH-clear AC support root found')
        roots.append(chosen.tolist());rear=rear+rod(tip,chosen,2.8)
    base=base+rear;hardware={}
    for side,sign in [('left',-1),('right',1)]:
        c=C+x*sign*(outer+2);caps=caps+rod(c+z*.1,c+z*capdepth,3.7);base=base+rod(c-z*4,c-z*.1,3.7)+rod(C+x*sign*(outer-2)-z*2,c-z*2,3);bore=rod(c-z*14,c+z*12,1.7);caps=caps-bore;base=base-bore
        head=capdepth+.75;bolt=rod(c+z*(head-(12 if COMPACT else 16)),c+z*head,1.5)+rod(c+z*head,c+z*(head+3),2.75);nut=mf.Manifold.cylinder(2.4,5.5/np.sqrt(3),5.5/np.sqrt(3),6).transform(np.c_[F,c-z*7.15])-rod(c-z*9,c-z*3,1.6);top=rod(c+z*(capdepth+.15),c+z*(capdepth+.65),3.5)-rod(c+z*(capdepth-.5),c+z*(capdepth+1.5),1.7);bottom=rod(c-z*4.65,c-z*4.15,3.5)-rod(c-z*5,c-z*3,1.7)
        hardware.update({'AC_bolt_'+side:bolt,'AC_nut_'+side:nut,'AC_top_washer_'+side:top,'AC_bottom_washer_'+side:bottom});base=base-(rod(c-z*7.4,c-z*4.5,3.45)+rod(c-z*4.9,c-z*3.95,3.75))
    clearance=mf.Manifold.batch_boolean([caps]+[caps.translate(d*s*.15) for d in [x,y,z] for s in [-1,1]],mf.OpType.Add);base=base-clearance;moving_clearance=mf.Manifold.batch_boolean([clav]+[clav.translate(d*s*.15) for d in [x,y,z] for s in [-1,1]],mf.OpType.Add);base=base-moving_clearance
    left=caps^box([-30,-30,-30],[-.1,30,30]);right=caps^box([.1,-30,-30],[30,30,30]);parts={};cleanup={}
    for name,s in {'clavicle_SC_AC_balls':clav,'scapula_AC_socket':base,'AC_cap_left':left,'AC_cap_right':right,**hardware}.items():
        s=s.as_original().simplify(.01);shells=s.decompose();meaningful=[a for a in shells if abs(a.volume())>1e-10]
        if len(meaningful)==1 and len(shells)>1:cleanup[name]=len(shells)-1;s=meaningful[0]
        parts[name]=s
    rows=[write(DEST/(n+'.stl'),s) for n,s in parts.items()];overlaps=[]
    for (an,aa),(bn,bb) in itertools.combinations(parts.items(),2):
        vol=max(0,float((aa^bb).volume()))
        if vol>.001:overlaps.append({'parts':[an,bn],'volume_mm3':vol})
    hum_hits=[{'part':n,'volume_mm3':float((s^hum).volume())} for n,s in parts.items() if (s^hum).volume()>.001];owners={n:'SC' if n=='clavicle_SC_AC_balls' else 'AC' for n in parts}
    report={'scope':'AC physical geometry candidate on SC assembly, not drive/strength acceptance','centre_source_mm':(C+SHIFT).tolist(),'ball_radius_mm':ball_radius,'socket_inner_radius_mm':inner,'socket_outer_radius_mm':outer,'neck_radius_mm':neck_radius,'capdepth_mm':capdepth,'bone_relief_radius_mm':relief,'socket_axis_toward_neck':z.tolist(),'support_roots_display_mm':roots,'components':{r['name']:r['components'] for r in rows},'part_owners':owners,'neutral_overlaps':overlaps,'neutral_humerus_overlaps':hum_hits,'zero_volume_shell_cleanup':cleanup,'surface_simplification_tolerance_mm':.01,'geometry_delta':{'clavicle_removed_mm3':float((oldclav-parts['clavicle_SC_AC_balls']).volume()),'clavicle_added_mm3':float((parts['clavicle_SC_AC_balls']-oldclav).volume()),'scapula_removed_mm3':float((oldscap-parts['scapula_AC_socket']).volume()),'scapula_added_mm3':float((parts['scapula_AC_socket']-oldscap).volume())},'frozen':False,'pending':['full-body SC/AC sweeps','assembly access and printed tolerances','muscle attachment remapping','loads and strength']};(DEST/'assembly.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf8');(DEST/'cad_input.json').write_text(json.dumps({'parts':rows,'units':'mm'}),encoding='utf8')
    tree=ET.parse(BASE/'arm.xml');r=tree.getroot();r.find('compiler').set('meshdir','.')
    for mesh in r.findall('./asset/mesh'):mesh.set('file',os.path.relpath((BASE/mesh.get('file')).resolve(),DEST).replace('\\','/'))
    for fid,name in [('FJ3362','clavicle_SC_AC_balls'),('FJ3384','scapula_AC_socket')]:r.find(f'./asset/mesh[@name="{fid}"]').set('file',name+'.stl')
    A=np.array(definition['centres_source_mm']['AC'])-SHIFT
    for name in parts:
        if name in ['clavicle_SC_AC_balls','scapula_AC_socket']:continue
        ET.SubElement(r.find('asset'),'mesh',name=name,file=name+'.stl',scale='.001 .001 .001');ET.SubElement(r.find('.//body[@name="AC"]'),'geom',name=name,type='mesh',mesh=name,pos=' '.join(str(v/1000) for v in -A),rgba='.3 .6 .75 1',group='3',contype='0',conaffinity='0',mass='.001')
    ET.indent(r);tree.write(DEST/'arm.xml',encoding='utf-8',xml_declaration=True);print('AC components',report['components'],'overlaps',overlaps,'humerus',hum_hits,flush=True)

if __name__=='__main__':run()
