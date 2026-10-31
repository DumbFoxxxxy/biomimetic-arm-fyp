"""Read actual GIF metadata; no claim of force-driven simulation."""

# Locate the project package from this script, independent of working directory.
import sys as _project_sys
from pathlib import Path as _NativeProjectPath
_project_root = next(p for p in _NativeProjectPath(__file__).resolve().parents if (p / 'project_runtime').is_dir())
if str(_project_root) not in _project_sys.path:
    _project_sys.path.insert(0, str(_project_root))
from project_runtime import load_mujoco_model as _load_project_model, load_mujoco_xml as _load_project_xml
from project_runtime import WorkspacePath as Path
from PIL import Image
import json,hashlib,sys
ROOT=Path(__file__).resolve().parents[2]
DEST=ROOT/'mujoco_arm/exports/anatomical_concept/elbow_muscles/mounted_elbow_candidate/side_lever_candidate/axis_negative_r15'/('swivel_terminal_distal_load' if '--distal-load' in sys.argv else 'swivel_terminal_candidate')
def run():
    path=DEST/'terminal_articulation.gif';im=Image.open(path);durations=[]
    for i in range(im.n_frames):im.seek(i);durations.append(im.info.get('duration',0))
    report={'scope':'prescribed local ball-terminal rotation via mocap; TPU held neutral, anatomy hidden in close-up, no forces/drive/deformation simulated','file':str(path.relative_to(ROOT)),'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'actual_frame_count':im.n_frames,'size_pixels':list(im.size),'total_duration_ms':sum(durations),'nominal_articulation_angles_deg':[-25,25],'sampled_axes':['selected_retainer_v','cross(neutral_terminal_axis,v)'],'source_cutaway_XML_sha256':hashlib.sha256((DEST/'neutral_cutaway_view.xml').read_bytes()).hexdigest()}
    if (DEST/'terminal_animation_inputs.json').exists():report['animation_inputs_sha256']=hashlib.sha256((DEST/'terminal_animation_inputs.json').read_bytes()).hexdigest()
    (DEST/'terminal_articulation.json').write_text(json.dumps(report,indent=2),encoding='utf8');print(report)
if __name__=='__main__':run()
