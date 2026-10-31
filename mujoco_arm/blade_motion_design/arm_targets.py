"""Full target scope retained after the user's GH concentric-mechanism choice."""
import json
import numpy as np
from shoulder_model import inputs,EXPORT

def target_actions():
    cfg,fit,anchors=inputs();zero=np.array(fit['poses'][0]['q_rad']);actions={}
    for name in ['UP_ROT','DOWN_ROT','PROTRACT','RETRACT','ELEVATE','DEPRESS']:
        seq=[zero]+[np.array(p['q_rad']) for p in fit['poses'] if p.get('action')==name]
        actions['SCAP_'+name]=np.array([np.r_[q,np.zeros(4)] for q in seq])
    for name,k,end in [('GH_FLEX',0,120.),('GH_EXT',0,-40.),('GH_ABD',1,120.),('GH_INT',2,70.),('GH_EXT_ROT',2,-90.),('ELBOW',3,145.)]:
        seq=[]
        for f in np.linspace(0.,1.,21):
            values=np.zeros(4);values[k]=np.radians(end*f);seq.append(np.r_[zero,values])
        actions[name]=np.array(seq)
    protract=actions['SCAP_PROTRACT'];cross=[]
    for i,f in enumerate(np.linspace(0,1,21)):
        cross.append(np.r_[protract[i,:5],np.radians([40*f,-20*f,0.,90*f])])
    actions['CROSS_BODY']=np.array(cross)
    endpoint=json.loads((EXPORT/'combined_raise_endpoint_search.json').read_text())
    goal=np.radians(endpoint['q_deg']);combined=[]
    for f in np.linspace(0.,1.,21):combined.append(np.r_[zero*(1-f)+goal[:5]*f,goal[5:7]*f,0.,0.])
    actions['RAISE_170']=np.array(combined)
    refined=EXPORT/'combined_raise_path_candidate.json'
    if refined.exists():
        candidate=json.loads(refined.read_text())
        if candidate['pass']:actions['RAISE_170']=np.array(candidate['points_rad'])
    return actions
