"""MuJoCo kinematic display of the partial scapula installation step.

One mocap bone, no articulated arm dynamics. Three liners are preplaced;
the GH/humerus/clavicle bone and muscles are explicitly not installed yet.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[2]
DEST = ROOT/'02_人体骨骼晶格肌肉方案/输出/shoulder_design/s1_fixture_assembly_review'
sys.path.insert(0, str(ROOT))
from project_runtime import load_mujoco_model


def hashed(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(render=False, check=False):
    import mujoco
    import numpy as np

    dense_file = DEST/'selected_path_dense_check.json'
    dense = json.loads(dense_file.read_text(encoding='utf8'))
    audit_file = DEST/'scapula_insertion_audit.json'
    if not dense['selected_partial_path_samples_clear']:
        raise ValueError('Dense path failed; do not display it as a passed installation step')
    if dense['audit_sha256'] != hashed(audit_file):
        raise ValueError('Changed insertion audit')
    audit = json.loads(audit_file.read_text(encoding='utf8'))
    if dense['assembly_sha256'] != hashed(ROOT/dense['source_assembly']):
        raise ValueError('Changed source assembly')
    scene_file = DEST/'selected_stage_scene.json'
    scene = json.loads(scene_file.read_text(encoding='utf8'))
    expected = dense['moving_keys']+dense['fixed_keys']
    if set(scene['sources']) != set(expected):
        raise ValueError('Partial scene membership differs from verified stage')
    for key, row in scene['sources'].items():
        if row != audit['source_manifest'][key] or hashed(ROOT/row['file']) != row['sha256']:
            raise ValueError('Changed exact scene mesh or placement: '+key)
    # Reuse the exact-triangle loader/large-STL conversion, not simplified
    # display meshes. It writes only this partial stage's view assets.
    subprocess.run([sys.executable, str(ROOT/'mujoco_arm/tools/view_shoulder_manufactured.py'),
                    '--report', str(scene_file), '--check'], check=True)
    xml = ET.parse(DEST/'available_parts.xml')
    root = xml.getroot()
    root.set('model', 'scapula_partial_assembly_mocap_NOT_actuation')
    ET.SubElement(root, 'option', gravity='0 0 0')
    world = root.find('worldbody')
    keys = list(scene['sources'])
    bone_index = keys.index('bone:scapula')
    bone_name = f'part_{bone_index:03d}'
    bone_geom = world.find(f"geom[@name='{bone_name}']")
    if bone_geom is None:
        raise ValueError('Exact scapula geometry missing')
    world.remove(bone_geom)
    body = ET.SubElement(world, 'body', name='assembly_scapula', mocap='true', pos='0 0 0')
    body.append(bone_geom)
    for i, key in enumerate(keys):
        if key.endswith('/front_liner'):
            world.find(f"geom[@name='part_{i:03d}']").set('rgba', '.12 .55 .85 1')
    visual = root.find('visual/global')
    visual.set('offwidth', '1200')
    visual.set('offheight', '900')
    model_file = DEST/'scapula_insertion_demo.xml'
    xml.write(model_file, encoding='utf-8', xml_declaration=True)
    m = load_mujoco_model(model_file)
    d = mujoco.MjData(m)
    direction = np.array(dense['selected_direction_source'])
    maximum = dense['maximum_translation_mm']
    gid = mujoco.mj_name2id(m, mujoco.mjtObj.mjOBJ_GEOM, bone_name)

    def pose(distance):
        d.mocap_pos[0] = direction*distance*.001
        mujoco.mj_forward(m, d)

    pose(0.)
    origin = d.geom_xpos[gid].copy()
    pose(maximum)
    observed = d.geom_xpos[gid]-origin
    if not np.allclose(observed, direction*maximum*.001, atol=1e-12):
        raise ValueError('MuJoCo mocap displacement differs from verified source path')

    def camera(cam):
        cam.lookat[:] = [-.12, -.01, 1.29]
        cam.distance = .52
        cam.azimuth = 130
        cam.elevation = 8

    metadata = dict(scope=__doc__, partial_scene_instances=m.ngeom, source_assembly_instances=audit['source_count'],
                    mujoco_version=mujoco.__version__, nq=m.nq, nmocap=m.nmocap,
                    exact_scene_sources_hash_verified=True, selected_dense_path_sha256=hashed(dense_file),
                    script_sha256=hashed(Path(__file__).resolve()), model_sha256=hashed(model_file),
                    observed_scapula_displacement_m=observed.tolist(),
                    displacement_matches_verified_path=True, display_collision_disabled=True,
                    used_as_collision_or_dynamic_proof=False, stage_pass=False)
    if render:
        from PIL import Image, ImageDraw
        width, height = 640, 650
        renderer = mujoco.Renderer(m, height=height, width=width)
        cam = mujoco.MjvCamera()
        camera(cam)
        canvas = Image.new('RGB', (width*3, height+85), 'white')
        draw = ImageDraw.Draw(canvas)
        draw.text((12, 10), f'SCAPULA INSERTION STEP | {m.ngeom} stage objects / {audit["source_count"]} source objects | MOCAP, NOT ACTUATION', fill='black')
        draw.text((12, 32), 'Blue: 3 preplaced loose liners (holding unverified). GH, humerus, clavicle bone and muscles not installed.', fill='black')
        for i, distance in enumerate([maximum, maximum/2, 0.]):
            pose(distance)
            renderer.update_scene(d, camera=cam)
            renderer.scene.flags[mujoco.mjtRndFlag.mjRND_SHADOW] = False
            canvas.paste(Image.fromarray(renderer.render()), (i*width, 85))
            draw.text((i*width+12, 61), f'Bone offset {distance:.1f} mm | assembly translation only', fill='black')
        renderer.close()
        output = DEST/'scapula_insertion_demo.png'
        canvas.save(output)
        metadata['render_sha256'] = hashed(output)
    (DEST/'view_check.json').write_text(json.dumps(metadata, indent=2)+'\n', encoding='utf8', newline='\n')
    print('Partial scene', m.ngeom, 'geoms; nq', m.nq, '; mocap', m.nmocap, flush=True)
    if render or check:
        return
    from mujoco import viewer
    state = dict(playing=False, distance=maximum, start=time.monotonic())

    def on_key(key):
        if key == 32:
            state['playing'] = not state['playing']
            state['start'] = time.monotonic()
        elif key in [48, 49]:
            state['playing'] = False
            state['distance'] = 0. if key == 48 else maximum

    print('SPACE: play/pause; 0: nominal installed pose; 1: withdrawn pose. No dynamics stepping.', flush=True)
    with viewer.launch_passive(m, d, key_callback=on_key) as handle:
        camera(handle.cam)
        while handle.is_running():
            if state['playing']:
                t = (time.monotonic()-state['start']) % 8.
                state['distance'] = maximum*(1-t/4) if t <= 4 else maximum*(t-4)/4
            pose(state['distance'])
            handle.sync()
            time.sleep(.02)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--render', action='store_true')
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    run(args.render, args.check)
