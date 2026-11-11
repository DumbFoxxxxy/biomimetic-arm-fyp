"""Search sampled SC/AC coordination using the current rigid assembly."""

from itertools import combinations
from pathlib import Path

import numpy as np
from scipy.spatial.transform import Rotation

from audit_s1_series_elbow_integration import MOBILE, compose, find_hits
from build_channel_first_inventory import DEST, ROOT, provenance, read
from build_s1_teres_origin_mount import A, export, load_component, save, shoulder_pads


def affine_at(matrix, pivot):
    return np.c_[matrix, pivot - matrix @ pivot]


class ShoulderAssembly:
    def __init__(self, sources, owners, origins):
        self.keys = list(sources)
        self.solids = [load_component(sources[key]) for key in self.keys]
        self.groups = [2 if owners[key] in MOBILE | {"scapula"} else int(owners[key] == "clavicle")
                       for key in self.keys]
        self.pairs = np.asarray([(i, j) for i, j in combinations(range(len(self.keys)), 2)
                                 if self.groups[i] != self.groups[j]])
        self.sc, self.ac, self.gh = [np.asarray(origins[key]) * 1000.
                                     for key in ("clavicle", "scapula", "humerus")]

    def kinematics(self, target, sc):
        clavicle = Rotation.from_rotvec(np.deg2rad(sc))
        target_rotation = Rotation.from_rotvec(np.deg2rad([0., target, 0.]))
        ac = clavicle.inv() * target_rotation
        clavicle_world = affine_at(clavicle.as_matrix(), self.sc)
        scapula_world = compose(clavicle_world, affine_at(ac.as_matrix(), self.ac))
        return clavicle_world, scapula_world, ac

    def evaluate(self, angle, sc):
        clavicle, scapula, ac = self.kinematics(angle, sc)
        if np.linalg.norm(sc) > 30. or np.rad2deg(ac.magnitude()) > 40.:
            return [], "search_bound"
        frames = (None, clavicle, scapula)
        posed = [solid if group == 0 else solid.transform(frames[group])
                 for solid, group in zip(self.solids, self.groups, strict=True)]
        return find_hits(posed, self.keys, self.pairs), None

    def transition(self, start, end, previous, current):
        sc1, sc2 = (Rotation.from_rotvec(np.deg2rad(vector)) for vector in (previous, current))
        ac1 = self.kinematics(start, previous)[2]
        ac2 = self.kinematics(end, current)[2]
        step_bound = max(np.rad2deg((sc1.inv() * sc2).magnitude()), np.rad2deg((ac1.inv() * ac2).magnitude()))
        subdivisions = max(2, int(np.ceil(step_bound / .5)))
        points = []
        for f in np.linspace(0., 1., subdivisions + 1)[1:]:
            angle = start + f * (end - start)
            sc = previous + f * (current - previous)
            hits, limit = self.evaluate(angle, sc)
            if hits or limit:
                return points, dict(target_deg=float(angle), SC_rotvec_deg=sc.tolist(), hits=hits, limit=limit)
            points.append(self.record(angle, sc))
        return points, None

    def record(self, angle, sc):
        clavicle, scapula, ac = self.kinematics(angle, sc)
        return dict(scapular_coronal_deg=float(angle), SC_rotvec_deg=sc.tolist(),
                    AC_relative_rotvec_deg=np.rad2deg(ac.as_rotvec()).tolist(),
                    AC_joint_displacement_mm=(clavicle[:, :3] @ self.ac + clavicle[:, 3] - self.ac).tolist(),
                    GH_displacement_mm=(scapula[:, :3] @ self.gh + scapula[:, 3] - self.gh).tolist())


def follow(model, limit, direction, coupled):
    angle = 0.
    sc = np.zeros(3)
    path = [model.record(angle, sc)]
    candidates = [np.zeros(3)]
    if coupled:
        values = [np.array([i, j, k], dtype=float) for i in (-1, 0, 1)
                  for j in (-1, 0, 1) for k in (-1, 0, 1) if i or j or k]
        candidates += sorted(values, key=np.linalg.norm)
    trial_count = 0
    for step in range(1, limit + 1):
        target = direction * step
        failures = []
        for increment in candidates:
            trial_count += 1
            next_sc = sc + increment
            samples, failure = model.transition(angle, target, sc, next_sc)
            if failure is None:
                path.extend(samples)
                angle, sc = target, next_sc
                break
            failures.append(failure)
        else:
            return dict(path=path, attempted_transitions=trial_count, rejected_next_step=failures)
    return dict(path=path, attempted_transitions=trial_count, rejected_next_step=[])


def main():
    folder = DEST.parent / "integration_review/s1_series_elbow"
    source_file = folder / "assembly_FJ1504_local_interface_bidirectional.json"
    source = read(source_file)
    origin_file = ROOT / "design/S1_CURRENT_CANDIDATE.json"
    _ = read(origin_file)
    origins = read(DEST / "elbow_passive_model/model_check.json")["neutral_origins_m"]
    key = "origin:FJ1504/bolt_0"
    routes = read(A.DESIGN / "s1_port_path_revision/routes_compact_distal_candidate.json")["routes"]
    port = next(row["ports"]["origin"] for row in routes if row["head"] == "FJ1504")
    outer = np.asarray(port["outer_face_source_mm"])
    normal = np.asarray(port["normal_source"])
    head_base = outer + (-9.5 + .75) * normal
    block = A.I._manifold().Manifold.cube([12., 12., 6.]).translate(
        [-6., -6., 1.1]).transform(np.c_[shoulder_pads.pad_basis(normal), head_base])
    existing = load_component(source["sources"][key])
    shortened, mesh = export((existing - block).as_original().simplify(1e-6),
                             folder / "FJ1504_bolt0_low_head_motion_candidate.obj", head_base)
    models = {
        "existing": ShoulderAssembly(source["sources"], source["owners"], origins),
        "low_head": ShoulderAssembly({**source["sources"], key: mesh}, source["owners"], origins),
    }
    cases = {}
    for name, model in models.items():
        if model.evaluate(0., np.zeros(3)) != ([], None):
            raise RuntimeError(f"Neutral rigid cross-body interference: {name}")
        cases[name] = {}
        for direction in (-1, 1):
            direct = follow(model, 36, direction, False)
            coordinated = follow(model, 36, direction, True)
            cases[name][str(direction)] = dict(AC_only=direct, SC_AC=coordinated)
            print(name, direction, "AC only", direct["path"][-1]["scapular_coronal_deg"],
                  "coupled", coordinated["path"][-1]["scapular_coronal_deg"], flush=True)
    save(folder / "coupled_scapula_screen.json", dict(
        generator=provenance(Path(__file__).resolve()), assembly=provenance(source_file),
        low_head=dict(source=mesh, height_mm=1.1, removed_mm3=float((existing - shortened).volume())),
        cases=cases, coronal_axis_source=[0., 1., 0.], search_sc_rotvec_limit_deg=30.,
        search_ac_rotation_limit_deg=40., step_target_deg=1.,
        scope="Greedy rigid kinematic path search. SC and AC share a prescribed net coronal rotation "
        "around source Y. The AC centre follows the actual rotated clavicle; no free translations. "
        "Uses ideal spherical joints, holds GH and elbow neutral, excludes flexible collision and "
        "unmodelled pads. Positive Y raises the lateral GH landmark at neutral if AC rotates alone. "
        "This path is not a K6 or drive acceptance. Low-head CAD is a geometric envelope only."))


if __name__ == "__main__":
    main()
