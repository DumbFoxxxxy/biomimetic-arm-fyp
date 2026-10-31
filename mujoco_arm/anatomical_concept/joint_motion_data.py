"""Accepted prescribed rigid paths only; no actuator or muscle simulation."""
import json
import numpy as np

def controls_q(c):
    c=np.array(c);q=np.zeros(13)
    for i,a in [(0,0),(3,4),(6,8)]:
        v=np.radians(c[i:i+3]);angle=np.linalg.norm(v);q[a]=np.cos(angle/2)
        if angle>1e-15:q[a+1:a+4]=np.sin(angle/2)*v/angle
    q[-1]=np.radians(c[9]);return q

def accepted_paths(folder):
    report=json.loads((folder/'bone_motion_audit.json').read_text(encoding='utf8'));paths=json.loads((folder/'motion_paths.json').read_text(encoding='utf8'));out={}
    for row in report['actions']:
        seq=[]
        for p in paths[row['probe']]:
            if p['hits']:break
            seq.append(np.array(p['qpos']))
        if seq:
            last=controls_q(row['last_clear_controls_deg'])
            if np.max(abs(last-seq[-1]))>1e-9:seq.append(last)
            out[row['probe']]=seq+seq[-2::-1]
    return out
