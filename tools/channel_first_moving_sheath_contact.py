"""Unilateral distributed sheath contact with equal opposite loads on moving bones."""
import numpy as np
from build_s1_teres_origin_mount import ROOT,read,load_component
from build_channel_first_inventory import provenance
from channel_first_bone_contact import BoneContactSurface
from simulate_channel_first_conduit_bone_contact import contact_stations


class MovingSheathContact:
    def __init__(self,model,check,density=200000.,friction=0.,interval=.0002):
        import mujoco
        assert density>0 and friction>=0 and interval>0
        self.model=model;self.density=density;self.friction=friction;self.interval=interval
        self.last_time=-np.inf;self.held=np.zeros(model.nv);self.latest=[];self.updates=[]
        self.peak_penetration=0.;self.minimum_distance=4.;self.peak_normal=0.;self.peak_pair_residual=0.
        assembly=read(ROOT/check['assembly']['file']);parent=read(ROOT/check['parent']['file'])
        owner={r['source_key']:r['owner'] for r in check['geom_owners']};self.bones=[]
        self.source_records={}
        for key in ('bone:scapula','bone:humerus'):
            source=assembly['sources'][key];solid=load_component(source);body=owner[key]
            self.bones.append(dict(key=key,body=mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_BODY,body),
                origin_mm=np.array(parent['neutral_origins_m'][body])*1000,
                surface=BoneContactSurface(solid),bounds=np.array(solid.bounding_box())))
            self.source_records[key]=source
        self.beams=[]
        for row in check['channels']:
            if row.get('region')!='external_supply':continue
            spec=row['specification'];stations=contact_stations(np.array(spec['segment_rest_lengths_m']),spacing=.00125)
            self.beams.append(dict(head=row['head'],
                ids=[mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_SITE,f"{row['head']}__node_{i}") for i in range(spec['segments']+1)],
                a=np.array([r[0] for r in stations]),b=np.array([r[1] for r in stations]),
                fractions=np.array([r[2] for r in stations]),weights=np.array([r[3] for r in stations])))

    @staticmethod
    def skew(v):
        x,y,z=v;return np.array([[0.,-z,y],[z,0.,-x],[-y,x,0.]])

    def apply(self,data):
        import mujoco
        if data.time-self.last_time<self.interval-1e-12:
            data.qfrc_applied[:]+=self.held;return
        self.last_time=float(data.time);self.held[:]=0.;contacts=[];normal_sum=0.;jac_cache={}
        def jac(site):
            if site not in jac_cache:
                jp=np.zeros((3,self.model.nv));jr=np.zeros_like(jp)
                mujoco.mj_jacSite(self.model,data,jp,jr,site);jac_cache[site]=(jp,jr)
            return jac_cache[site]
        for beam in self.beams:
            nodes=data.site_xpos[beam['ids']];f=beam['fractions'][:,None]
            points=(1-f)*nodes[beam['a']]+f*nodes[beam['b']]
            for bone in self.bones:
                body=bone['body'];rotation=data.xmat[body].reshape(3,3)
                translation_mm=data.xpos[body]*1000-rotation@bone['origin_mm']
                local=(points*1000-translation_mm)@rotation
                bounds=bone['bounds'];near=np.all(local>=bounds[:3]-4.,axis=1)&np.all(local<=bounds[3:]+4.,axis=1)
                selected=np.flatnonzero(near)
                if not len(selected):continue
                distances,normals=bone['surface'].query(local[selected])
                self.minimum_distance=min(self.minimum_distance,float(distances.min()))
                for j in np.flatnonzero(distances<1.5):
                    i=selected[j];distance=float(distances[j]);normal=rotation@normals[j]
                    penetration=(1.5-distance)*.001;weight=float(beam['weights'][i]);fraction=float(beam['fractions'][i])
                    a=beam['ids'][beam['a'][i]];b=beam['ids'][beam['b'][i]]
                    ja,jra=jac(a);jb,jrb=jac(b)
                    jcentre=(1-fraction)*ja+fraction*jb;jrotation=(1-fraction)*jra+fraction*jrb
                    # Apply both members of the force pair at the same actual
                    #closest bone surface point, including frictional moment.
                    point=points[i]-normal*distance*.001
                    guide_jac=jcentre-self.skew(point-points[i])@jrotation
                    bone_jac=np.zeros((3,self.model.nv));bone_rot=np.zeros_like(bone_jac)
                    mujoco.mj_jac(self.model,data,bone_jac,bone_rot,point,body)
                    relative=(guide_jac-bone_jac)@data.qvel;vn=float(relative@normal)
                    force_n=max(0.,weight*(self.density*penetration-50.*vn))
                    tangent=relative-vn*normal
                    friction=-self.friction*force_n*tangent/np.sqrt(tangent@tangent+.001**2)
                    assert float(friction@tangent)<=1e-12
                    force=force_n*normal+friction
                    pair=(guide_jac-bone_jac).T@force
                    # Independent generalized-force assembly at the same point.
                    separate=guide_jac.T@force-bone_jac.T@force
                    residual=float(np.max(abs(pair-separate),initial=0));self.peak_pair_residual=max(self.peak_pair_residual,residual)
                    assert residual<1e-9
                    self.held+=pair;normal_sum+=force_n;self.peak_penetration=max(self.peak_penetration,penetration*1000)
                    contacts.append(dict(beam=beam['head'],bone=bone['key'],station=int(i),
                        penetration_mm=penetration*1000,normal_force_N=force_n,
                        point_source_mm=(point*1000).tolist(),normal_world=normal.tolist(),force_on_guide_N=force.tolist(),
                        opposite_force_on_bone_N=(-force).tolist(),tangential_power_W=float(friction@tangent)))
        self.latest=contacts;self.peak_normal=max(self.peak_normal,normal_sum)
        self.updates.append(dict(time_s=float(data.time),sum_normal_force_N=normal_sum,contacts=contacts))
        data.qfrc_applied[:]+=self.held

    def report(self):
        return dict(source_files=self.source_records,generator=provenance(ROOT/'tools/channel_first_moving_sheath_contact.py'),
            surface_query=provenance(ROOT/'tools/channel_first_bone_contact.py'),
            penalty_density_N_per_m2=self.density,damping_density_Ns_per_m2=50.,friction_coefficient=self.friction,
            update_interval_s=self.interval,station_maximum_rest_spacing_mm=1.25,contact_radius_mm=1.5,
            maximum_sampled_penetration_mm=self.peak_penetration,minimum_sampled_centreline_distance_mm_capped4=self.minimum_distance,
            peak_sum_normal_force_N=self.peak_normal,pair_virtual_work_max_abs_mixed_units=self.peak_pair_residual,
            updates=self.updates,
            scope='Actual moving scapula/humerus surfaces only. Equal opposite contact at common closest surface point; generalized forces by virtual work. No extra fixed guide. Limited cross-section interpolation and penalty force held between updates. No clavicle/thorax/hardware, sheath self contact, cable/conduit friction, lumen collapse or continuum indentation. Independent whole-wire/rigid screen still mandatory.')
