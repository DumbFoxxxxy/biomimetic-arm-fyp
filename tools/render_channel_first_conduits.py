"""Render actual static guide meshes without optional image libraries."""
import struct
import sys
import zlib

import numpy as np

from build_s1_teres_origin_mount import ROOT, read, save
from build_channel_first_inventory import provenance


def png(path, pixels):
    height, width, _ = pixels.shape
    def chunk(kind, data):
        return struct.pack('>I', len(data))+kind+data+struct.pack('>I', zlib.crc32(kind+data)&0xffffffff)
    data = b''.join(b'\0'+row.tobytes() for row in pixels)
    path.write_bytes(b'\x89PNG\r\n\x1a\n'+chunk(b'IHDR', struct.pack('>IIBBBBB', width, height, 8, 2, 0, 0, 0))+
                     chunk(b'IDAT', zlib.compress(data))+chunk(b'IEND', b''))


def main():
    import mujoco
    sys.path.insert(0, str(ROOT))
    from project_runtime import load_mujoco_model
    folder = ROOT/'02_人体骨骼晶格肌肉方案/输出/channel_first'
    report_file = folder/'minimal_shoulder_conduits_scene.json'
    report = read(report_file)
    model_file = folder/'available_parts.xml'
    model = load_mujoco_model(model_file)
    assert model.ngeom == len(report['sources'])
    for i, (key, source) in enumerate(report['sources'].items()):
        assert provenance(ROOT/source['file'])['sha256'] == source['sha256']
        model.geom_rgba[i] = [1., .09, .04, 1.] if key.startswith('channel:') else [.82, .83, .78, .35]
    data = mujoco.MjData(model)
    mujoco.mj_forward(model, data)
    renderer = mujoco.Renderer(model, height=900, width=1200)
    outputs = []
    for azimuth in (55., 145.):
        cam = mujoco.MjvCamera()
        cam.lookat[:] = [-.14, -.045, 1.29]
        cam.distance, cam.azimuth, cam.elevation = .43, azimuth, 8.
        renderer.update_scene(data, camera=cam)
        renderer.scene.flags[mujoco.mjtRndFlag.mjRND_SHADOW] = False
        destination = folder/f'minimal_conduits_static_{int(azimuth)}.png'
        png(destination, renderer.render())
        outputs.append(dict(image=provenance(destination), azimuth_deg=azimuth,
                            distance_m=cam.distance, lookat_m=cam.lookat.tolist()))
    renderer.close()
    save(folder/'minimal_conduits_render.json', dict(scene=provenance(report_file),
        model=provenance(model_file), images=outputs,
        scope='Actual four hollow meshes in static neutral thorax/arm reference. Anatomy translucent, channels red. No mounted fittings, deformation or dynamic approval.'))
    print('Rendered two static views; anatomy translucent, actual guides red')


if __name__ == '__main__':
    main()
