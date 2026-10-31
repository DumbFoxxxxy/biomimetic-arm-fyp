"""Check actual quaternion sampling and compiled zero-pose mesh placement."""
from audit_joint_model import *
from scipy.spatial import cKDTree

def run():
    dest=OUT/'joint_model/SC_connection'
    if WITH_AC:dest=dest/('AC_compact' if '--compact' in sys.argv else 'AC_connection')
    if WITH_GH:dest=DEST
    report=json.loads((dest/'bone_motion_audit.json').read_text(encoding='utf8'));paths=json.loads((dest/'motion_paths.json').read_text(encoding='utf8'));maximum=0.;norm_error=0.
    for path in paths.values():
        q=np.array([r['qpos'] for r in path])
        for a in [0,4,8]:
            norm_error=max(norm_error,float(np.max(abs(np.linalg.norm(q[:,a:a+4],axis=1)-1))));dot=np.sum(q[:-1,a:a+4]*q[1:,a:a+4],axis=1);maximum=max(maximum,float(np.max(np.degrees(2*np.arccos(np.clip(abs(dot),0,1))))))
        maximum=max(maximum,float(np.max(abs(np.degrees(np.diff(q[:,-1]))))))
    assert maximum<=1+1e-6 and norm_error<1e-8
    m=mujoco.MjModel.from_xml_path(str(dest/'arm.xml'));d=mujoco.MjData(m);mujoco.mj_forward(m,d);inv=json.loads((OUT/'inventory.json').read_text(encoding='utf8'));byid={r['id']:r for r in inv['models']};assembly=json.loads((dest/'assembly.json').read_text(encoding='utf8'));rows=[]
    part_paths={n:dest/(n+'.obj') for n in assembly['components']}
    if WITH_AC:
        sc=json.loads((SC_DEST/'assembly.json').read_text(encoding='utf8'))
        part_paths.update({n:SC_DEST/(n+'.obj') for n in sc['components'] if n!='clavicle_SC_ball'})
    if WITH_GH:
        ac=json.loads((AC_DEST/'assembly.json').read_text(encoding='utf8'))
        part_paths.update({n:AC_DEST/(n+'.obj') for n in ac['components'] if n!='scapula_AC_socket'})
    if WITH_ELBOW:
        gh=json.loads((GH_DEST/'assembly.json').read_text(encoding='utf8'))
        part_paths.update({n:GH_DEST/(n+'.obj') for n in gh['components'] if n!='humerus_GH_ball'})
    aliases={'FJ3362':'clavicle_SC_AC_balls' if WITH_AC else 'clavicle_SC_ball','FJ3290':'manubrium_socket_base'}
    if WITH_AC:aliases['FJ3384']='scapula_AC_socket'
    if WITH_GH:aliases.update(FJ3362='clavicle_SC_AC_balls',FJ3384='scapula_AC_GH_sockets',FJ3368='humerus_GH_ball')
    if WITH_ELBOW:aliases.update(FJ3368='humerus_GH_elbow_fork',FJ3391='forearm_elbow_hub')
    bones=['FJ3362','FJ3384','FJ3368','FJ3391','FJ3290'] if WITH_ELBOW else ['FJ3362','FJ3384','FJ3368','FJ3391','FJ3349','FJ3290']
    for geom in [*bones,*[n for n in part_paths if n not in aliases.values()]]:
        gid=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_GEOM,geom);mid=int(m.geom_dataid[gid]);start=int(m.mesh_vertadr[mid]);n=int(m.mesh_vertnum[mid]);v=m.mesh_vert[start:start+n];compiled=v@d.geom_xmat[gid].reshape(3,3).T+d.geom_xpos[gid]
        name=aliases.get(geom,geom)
        if name in part_paths:expected,_=read_obj(str(part_paths[name]));expected=expected/1000
        else:
            expected,f=read_obj(str(ROOT/byid[geom]['source']));expected=(expected[np.unique(f)]-SHIFT)/1000
        err=max(float(cKDTree(expected).query(compiled)[0].max()),float(cKDTree(compiled).query(expected)[0].max()));rows.append({'geometry':geom,'compiled_mesh_vs_source_max_deviation_mm':err*1000,'pass':err<1e-6})
    assert all(r['pass'] for r in rows);assert not report['neutral']['hits']
    result={'maximum_sampled_joint_increment_deg':maximum,'maximum_quaternion_norm_error':norm_error,'zero_pose_geometry':rows,'neutral_bone_and_hardware_overlap_free':True,'scope':'sampling and scene registration; not mechanics or physical-drive validation'};(dest/'audit_checks.json').write_text(json.dumps(result,indent=2),encoding='utf8');print('actual joint step',maximum,'zero pose meshes pass',len(rows),'max mm',max(r['compiled_mesh_vs_source_max_deviation_mm'] for r in rows))

if __name__=='__main__':run()
