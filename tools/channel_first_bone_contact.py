"""Narrow-band nearest-triangle contact against checked, closed bone meshes.

The voxel cache stores candidate triangle indices, never an interpolated SDF.
For each radius bucket, centroid search includes band + cell circumradius +
maximum triangle circumradius. This contains every triangle within the band
of any query point in the cell. Near-surface distances are evaluated afresh.
"""
import sys

import numpy as np
from scipy.spatial import cKDTree

from build_s1_teres_origin_mount import ROOT, mesh_geom

sys.path.insert(0, str(ROOT/'mujoco_arm/anatomical_concept'))
from extract_attachment_regions import closest


class BoneContactSurface:
    def __init__(self, solid, band_mm=4., cell_mm=2.):
        self.solid, self.band, self.cell = solid, band_mm, cell_mm
        vertices, faces = mesh_geom.manifold_mesh_to_arrays(solid.to_mesh64())
        self.vertices, self.faces = vertices, faces
        self.triangles = vertices[faces]
        self.centres = self.triangles.mean(axis=1)
        radii = np.linalg.norm(self.triangles-self.centres[:, None], axis=2).max(axis=1)
        cross = np.cross(self.triangles[:, 1]-self.triangles[:, 0],
                         self.triangles[:, 2]-self.triangles[:, 0])
        norm = np.linalg.norm(cross, axis=1)
        assert norm.min() > 1e-16
        self.face_normals = cross/norm[:, None]
        self.vertex_normals = np.zeros_like(vertices)
        self.edge_normals = {}
        for corner in range(3):
            a = self.triangles[:, (corner+1)%3]-self.triangles[:, corner]
            b = self.triangles[:, (corner+2)%3]-self.triangles[:, corner]
            angle = np.arctan2(np.linalg.norm(np.cross(a, b), axis=1), np.sum(a*b, axis=1))
            np.add.at(self.vertex_normals, faces[:, corner], angle[:, None]*self.face_normals)
            for triangle, normal in zip(faces, self.face_normals):
                key = tuple(sorted((int(triangle[corner]), int(triangle[(corner+1)%3]))))
                self.edge_normals[key] = self.edge_normals.get(key, np.zeros(3))+normal
        self.vertex_normals /= np.maximum(np.linalg.norm(self.vertex_normals, axis=1)[:, None], 1e-30)
        self.edge_normals = {key: normal/max(np.linalg.norm(normal), 1e-30)
                             for key, normal in self.edge_normals.items()}
        self.groups = []
        bins = np.digitize(radii, [.5, 1., 2., 4., 8., 16., 32., 64.])
        for bucket in np.unique(bins):
            indices = np.where(bins == bucket)[0]
            self.groups.append((indices, cKDTree(self.centres[indices]), float(radii[indices].max())))
        self.cache = {}

    def _cell(self, point):
        key = tuple(np.floor(point/self.cell).astype(int))
        if key not in self.cache:
            centre = (np.array(key)+.5)*self.cell
            candidates = []
            for indices, tree, radius in self.groups:
                local = tree.query_ball_point(centre, self.band+np.sqrt(3)*self.cell/2.+radius+1e-7)
                if local:
                    candidates.append(indices[local])
            self.cache[key] = dict(centre=centre,
                indices=np.concatenate(candidates) if candidates else np.array([], dtype=int))
        return self.cache[key]

    def ray_inside(self, point):
        diagonal = np.linalg.norm(self.vertices.max(axis=0)-self.vertices.min(axis=0))
        direction = np.array([1., .371, .193])
        direction /= np.linalg.norm(direction)
        hits = self.solid.ray_cast(point, point+direction*(diagonal*4.+1000.))
        parameters = np.unique(np.round([hit.distance for hit in hits], 10))
        return len(parameters)%2 == 1

    def nearest(self, point):
        cell = self._cell(point)
        candidates = cell['indices']
        if len(candidates):
            q, distance, local_index = closest(point, self.triangles[candidates])
        else:
            distance = np.inf
        if distance > self.band:
            if 'inside' not in cell:
                cell['inside'] = self.ray_inside(cell['centre'])
            if cell['inside']:
                raise ValueError('Bone penetration deeper than the contact query band; reject this run')
            return self.band, np.zeros(3)
        index = int(candidates[local_index])
        triangle, face = self.triangles[index], self.faces[index]
        pseudo = self.face_normals[index]
        vertex_distance = np.linalg.norm(triangle-q, axis=1)
        if vertex_distance.min() < 1e-7:
            pseudo = self.vertex_normals[face[vertex_distance.argmin()]]
        else:
            for a, b in ((0, 1), (1, 2), (2, 0)):
                delta = triangle[b]-triangle[a]
                t = np.clip((q-triangle[a])@delta/(delta@delta), 0., 1.)
                if np.linalg.norm(q-(triangle[a]+t*delta)) < 1e-7:
                    pseudo = self.edge_normals[tuple(sorted((int(face[a]), int(face[b]))))]
                    break
        inside = (point-q)@pseudo < 0.
        sign = -1. if inside else 1.
        normal = sign*(point-q)/distance if distance > 1e-10 else pseudo
        return sign*distance, normal

    def query(self, points):
        result = [self.nearest(point) for point in points]
        return np.array([item[0] for item in result]), np.array([item[1] for item in result])
