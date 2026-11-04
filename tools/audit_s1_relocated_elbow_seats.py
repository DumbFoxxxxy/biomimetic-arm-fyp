"""Measure bone stock along the relocated elbow adapter screw axes."""
from pathlib import Path

import numpy as np

from build_s1_teres_origin_mount import A, ROOT, load_component, read, save, sha


def centre(solid):
    return np.array(solid.bounding_box()).reshape(2, 3).mean(axis=0)


def main():
    folder = ROOT / '02_人体骨骼晶格肌肉方案/输出/integration_review/s1_elbow_shared_humerus'
    manifest = folder / 'assembly_tpu_readback_candidate.json'
    channels = folder / 'elbow_relocated_channel_cores/channel_core_check.json'
    frame_file = A.DESIGN / 's1_retained_joint_comparison/shaft_collar/station_60/interface_check.json'
    sources = read(manifest)['sources']
    selected = read(channels)
    if selected['assembly_sha256'] != sha(manifest):
        raise ValueError('Channel route was screened against a different assembly')
    choice = selected['choice']['pair']
    frame = read(frame_file)['frame_source']
    axial = np.array(frame['basis_columns'])[:, 2]
    bone = load_component(sources['bone:humerus'])
    rows = []
    for muscle, head in [('flexor', 'FJ1486'), ('extensor', 'FJ1477')]:
        transform = A.I.rotation_about(frame['origin_mm'], axial, choice[f'{muscle}_angle_deg'])
        transform[:, 3] += axial * choice[f'{muscle}_shift_mm']
        for index in range(2):
            prefix = f'elbow:{head}_bone_{index}'
            front, back = [
                centre(load_component(sources[f'{prefix}_{side}_washer']).transform(transform))
                for side in ('front', 'back')
            ]
            normal = A.D.unit(front - back)
            tangent = A.D.unit(np.cross(normal, axial))
            bitangent = np.cross(normal, tangent)
            rays = []
            for radius in (0., 1.7, 4.5):
                for angle in np.linspace(0., 2*np.pi, 16 if radius else 1, endpoint=False):
                    offset = radius * (np.cos(angle)*tangent + np.sin(angle)*bitangent)
                    origin = front + offset
                    hits = bone.ray_cast(origin + normal*100, origin - normal*100)
                    crossings = [float((np.array(hit.position) - origin) @ normal) for hit in hits]
                    rays.append(dict(radius_mm=radius, angle_deg=float(np.degrees(angle)),
                                     crossings_relative_to_front_washer_mm=crossings))
            # Preserve screw direction. Allow the proximal extensor seat to
            # move axially away from the shoulder backing plate.
            candidates = []
            axial_steps = range(-12, 13) if muscle == 'extensor' and index == 0 else [0]
            offsets = sorted(((x, y) for x in range(-50, 51) for y in axial_steps),
                             key=lambda pair: pair[0]**2 + pair[1]**2)
            along = bitangent if bitangent @ axial > 0 else -bitangent
            for step, axial_step in offsets:
                delta = .5 * (step * tangent + axial_step * along)
                bands = []
                for ray in rays:
                    angle = np.radians(ray['angle_deg'])
                    offset = ray['radius_mm'] * (np.cos(angle)*tangent + np.sin(angle)*bitangent)
                    point = front + delta + offset
                    hits = bone.ray_cast(point + normal*100, point - normal*100)
                    if len(hits) != 2:
                        break
                    bands.append(sorted(float((np.array(hit.position)-point) @ normal) for hit in hits))
                if len(bands) == len(rays):
                    back_faces, front_faces = np.array(bands).T
                    candidates.append(dict(
                        transverse_shift_mm=step*.5, axial_shift_mm=axial_step*.5,
                        axis_point_source_mm=(front+delta).tolist(),
                        minimum_sampled_bone_thickness_mm=float(min(front_faces-back_faces)),
                        front_boss_face_relative_mm=float(max(front_faces)+.5),
                        back_boss_face_relative_mm=float(min(back_faces)-.5)))
            candidates = [item for item in candidates if item['minimum_sampled_bone_thickness_mm'] >= 8.]
            rows.append(dict(muscle=muscle, anchor=index, front_washer_source_mm=front.tolist(),
                             back_washer_source_mm=back.tolist(), axis=normal.tolist(),
                             washer_spacing_mm=float((front-back) @ normal), samples=rays,
                             nearest_footprint_candidate=candidates[0] if candidates else None,
                             footprint_candidates=candidates))
    save(folder / 'elbow_relocated_seat_audit.json', dict(
        assembly_sha256=sha(manifest), channel_report_sha256=sha(channels),
        frame_sha256=sha(frame_file), script_sha256=sha(Path(__file__)), anchors=rows,
        scope='M3 bore and boss ray samples. Candidate requires 8 mm sampled bone thickness as a geometry screen only. Full contact, drilling and assembly unverified.'))
    for row in rows:
        print(row['muscle'], row['anchor'], round(row['washer_spacing_mm'], 2),
              [round(value, 2) for value in row['samples'][0]['crossings_relative_to_front_washer_mm']],
              sum(bool(sample['crossings_relative_to_front_washer_mm']) for sample in row['samples']),
              row['nearest_footprint_candidate'])


if __name__ == '__main__':
    main()
