"""Smooth centreline and indexed hollow sweep for channel prototypes (mm)."""
import numpy as np
from scipy.interpolate import CubicSpline

from build_s1_teres_origin_mount import A


def smooth_route(points, iterations, step=.5):
    arc = np.r_[0., np.cumsum(np.linalg.norm(np.diff(points, axis=0), axis=1))]
    parameter = np.linspace(0., arc[-1], int(np.ceil(arc[-1]/step))+1)
    original = np.column_stack([np.interp(parameter, arc, points[:, j]) for j in range(3)])
    smooth = original.copy()
    for _ in range(iterations):
        smooth[1:-1] = .25*smooth[:-2]+.5*smooth[1:-1]+.25*smooth[2:]
    curve = CubicSpline(parameter, smooth)
    fine = np.linspace(0., arc[-1], int(np.ceil(arc[-1]/.1))+1)
    velocity, acceleration = curve(fine, 1), curve(fine, 2)
    radii = np.linalg.norm(velocity, axis=1)**3/np.maximum(
        np.linalg.norm(np.cross(velocity, acceleration), axis=1), 1e-14)
    tangents = curve(parameter, 1)
    tangents /= np.linalg.norm(tangents, axis=1)[:, None]
    endpoint_directions = np.array([points[1]-points[0], points[-1]-points[-2]])
    endpoint_directions /= np.linalg.norm(endpoint_directions, axis=1)[:, None]
    return smooth, tangents, dict(iterations=iterations, parameter_step_mm=step,
        curvature_sample_parameter_step_mm=float(fine[1]-fine[0]),
        minimum_sampled_bend_radius_mm=float(radii.min()),
        maximum_same_parameter_displacement_mm=float(np.linalg.norm(smooth-original, axis=1).max()),
        endpoint_tangent_change_deg=np.rad2deg(np.arccos(np.clip(
            np.sum(endpoint_directions*tangents[[0, -1]], axis=1), -1., 1.))).tolist())


def hollow_sweep(points, tangents, outer_radius=1.6, inner_radius=.8, sides=48):
    normals = []
    reference = np.eye(3)[np.argmin(abs(tangents[0]))]
    for tangent in tangents:
        normal = reference-(reference@tangent)*tangent
        normal /= np.linalg.norm(normal)
        normals.append(normal)
        reference = normal
    normals = np.array(normals)
    binormals = np.cross(tangents, normals)
    theta = np.arange(sides)*2*np.pi/sides
    offsets = np.cos(theta)[None, :, None]*normals[:, None]+np.sin(theta)[None, :, None]*binormals[:, None]
    origin = points.mean(axis=0)
    local = points-origin
    vertices = np.vstack([(local[:, None]+outer_radius*offsets).reshape(-1, 3),
                          (local[:, None]+inner_radius*offsets).reshape(-1, 3)])
    count = len(points)*sides
    faces = []
    for i in range(len(points)-1):
        for j in range(sides):
            a, b = i*sides+j, i*sides+(j+1)%sides
            c, d = b+sides, a+sides
            faces.extend([(a, b, c), (a, c, d),
                          (a+count, c+count, b+count), (a+count, d+count, c+count)])
    for j in range(sides):
        a, b = j, (j+1)%sides
        faces.extend([(a, a+count, b+count), (a, b+count, b)])
        a, b = (len(points)-1)*sides+j, (len(points)-1)*sides+(j+1)%sides
        faces.extend([(a, b+count, a+count), (a, b, b+count)])
    faces = np.asarray(faces, dtype=np.uint32)
    solid = A.I._manifold().Manifold(A.I._manifold().Mesh64(vertices, faces))
    assert str(solid.status()).endswith('NoError') and solid.volume() > 0
    return solid, origin
