"""Translate historical record locations to real project directories.

Historical evidence is byte-preserved. Newly generated paths are canonical.
The only filesystem operations here are normal pathlib operations: no junctions,
symlinks, environment rewrites, or monkey-patching of the standard library.
"""
import os
import sys
from pathlib import Path as NativePath

PROJECT_ROOT = NativePath(__file__).resolve().parents[1]
OLD = '01_工程骨架方案'
NEW = '02_人体骨骼晶格肌肉方案'
ROUTES = {
    'final_model': OLD+'/SolidWorks', 'model': OLD+'/早期模型',
    'archive': OLD+'/历史归档', '骨骼': NEW+'/输入/骨骼',
    'muscle': NEW+'/输入/muscle',
}
for name in ['anatomical_concept','anatomical_arm']:
    ROUTES['mujoco_arm/exports/'+name] = NEW+'/输出/'+name
for name in ['blade_motion_design','cable_design_v2','cable_design_v3',
             'cable_design_v4','cable_design_v5','cable_design_v6',
             'cable_design_v7','cable_design_v8','k6_baseline',
             'skeleton_from_v2','tendon_demo','tpu_muscle_design','wire_psm_from_v2s1']:
    ROUTES['mujoco_arm/exports/'+name] = OLD+'/输出/'+name


def _location(parts):
    value = os.fspath(NativePath(*parts)).replace('\\','/')
    workspace = PROJECT_ROOT.as_posix()
    absolute = False
    # Historical artifacts from this workstation are rebased in any clone.
    for prefix in [workspace, 'D:/XJTLU/PSP']:
        if value.casefold() == prefix.casefold(): return str(PROJECT_ROOT)
        if value.casefold().startswith(prefix.casefold()+'/'):
            value=value[len(prefix)+1:]; absolute=True; break
    for source,target in sorted(ROUTES.items(),key=lambda row:-len(row[0])):
        if value==source or value.startswith(source+'/'):
            value=target+value[len(source):];break
    return str(PROJECT_ROOT/NativePath(value)) if absolute else value


class WorkspacePath(type(NativePath())):
    """Concrete pathlib path that resolves archived project addresses on join.

    Usable under CAD Python3.10 and simulation Python3.14. Plain labels and
    unrelated external absolute paths are not changed.
    """
    def __new__(cls,*parts,**kwargs):
        if sys.version_info < (3,12):
            return super().__new__(cls,_location(parts),**kwargs)
        return super().__new__(cls,*parts,**kwargs)

    def __init__(self,*parts,**kwargs):
        if sys.version_info >= (3,12): super().__init__(_location(parts),**kwargs)

    def with_segments(self,*parts): return type(self)(*parts)

    def __truediv__(self,other):
        # Python3.10's internal parsed-parts constructors bypass __new__.
        return type(self)(str(super().__truediv__(other)))

    def __rtruediv__(self,other): return type(self)(other,str(self))

    def joinpath(self,*parts): return type(self)(str(super().joinpath(*parts)))


def project_path(value):
    p=WorkspacePath(value)
    return p if p.is_absolute() else WorkspacePath(PROJECT_ROOT)/p
