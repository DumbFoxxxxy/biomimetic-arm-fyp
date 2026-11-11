"""Screen a passive shortening-to-sideways TPU restoring law for FJ1500."""

import numpy as np
from scipy.spatial.transform import Rotation

from build_channel_first_inventory import DEST, ROOT, provenance, read
from build_s1_teres_origin_mount import load_component, save
from check_shoulder_front_liner_keepers import is_compliant
from study_s1_infraspinatus_cap_relief import HEAD
from study_s1_infraspinatus_tpu_bias import collision_state, rest_geometry
from study_s1_local_force_geometry import load_reference
from study_s1_y_terminal_equilibrium import solve_branch


def main():
    folder = DEST.parent / "s1_local_muscle_load/FJ1500_cap_relief"
    path_file = folder / "FJ1500_flexion_path.json"
    reference = read(path_file)
    rest = load_reference(reference["rest_geometry"])
    screen = load_reference(rest["bias_screen"])
    previous = load_reference(screen["reference_equilibrium"])
    geometry = load_reference(previous["local_geometry"])
    assembly = load_reference(previous["shoulder_assembly"])
    model = next(row for row in load_reference(previous["kinematics"])["models"] if row["head"] == HEAD)
    ports = next(row["ports"] for row in load_reference(geometry["port_geometry"])["routes"] if row["head"] == HEAD)
    sleeve = load_reference(reference["sleeve"])
    sources = dict(assembly["sources"])
    key = f"link:{HEAD}/muscle_uncrimped_sleeve"
    if sources[key] != sleeve["original"]:
        raise ValueError("Compact sleeve and S1 assembly do not match")
    sources[key] = sleeve["selected"]["source"]
    solids = {key: load_component(source) for key, source in sources.items()
              if not is_compliant(key) and not key.startswith(("bench:", "fixture:", "clamp:", "clavicle:"))
              and not key.endswith("/tendon_wire")}
    centre = np.array(geometry["gh_centre_source_mm"])
    flexion = np.array(geometry["axes_source"][geometry["axes"].index("flexion")])
    pin_a, pin_b = (np.array(model[f"neutral_{end}_pin_source_mm"]) for end in ("cap", "bone"))
    axis_a, axis_b = (np.array(model[f"{end}_pin_axis_source"]) for end in ("cap", "bone"))
    offset = np.array(rest["rest_offset_source_mm"])
    direction = offset / np.linalg.norm(offset)
    adjusted, adjusted_ports = rest_geometry(model, ports, offset)
    probes = (0., 30., 60., 64., 68., 75., 85., 96.)
    loads = [0., 5., 10.]
    rows = []

    def evaluate(angle, gain, input_loads):
        gh = Rotation.from_rotvec(flexion * np.deg2rad(angle))
        equilibria = solve_branch(adjusted, adjusted_ports, centre, gh, screen["stiffness"], input_loads,
                                  retain_leads=True, lateral_coupling=(direction, gain))
        result = []
        for state in equilibria:
            pin, hits = collision_state(solids, assembly["intended_body_owners"], state, offset, gh,
                                        centre, pin_a, pin_b, axis_a, axis_b)
            result.append(dict(coupling_gain=gain, GH_deg=angle, cap_pin_source_mm=pin.tolist(), hits=hits, **state))
        return result

    statistics = []
    for gain in (0., -.25, .25, .5, .75, 1.):
        case = [row for angle in probes for row in evaluate(angle, gain, loads)]
        if gain == 0.:
            baseline = {(r["angle_deg"], r["input_N"]): r for r in reference["outward"]}
            for row in case:
                prior = baseline[row["GH_deg"], row["commanded_muscle_tension_N"]]
                if (not np.allclose(row["cap_pin_source_mm"], prior["cap_pin_source_mm"], atol=1e-7, rtol=0.)
                        or bool(row["hits"]) != bool(prior["hits"])):
                    raise ValueError("Uncoupled equilibrium did not reproduce previous baseline")
        rows.extend(case)
        failure = next((row["GH_deg"] for row in case if row["hits"]), None)
        metrics = dict(gain=gain, failing_count=sum(bool(row["hits"]) for row in case), first_failing_probe_deg=failure)
        statistics.append(metrics)
        print(metrics, flush=True)

    candidate = max(statistics, key=lambda row: (row["first_failing_probe_deg"] or 180.,
                                                 -row["failing_count"], -abs(row["gain"])))
    gain = candidate["gain"]
    # The parameter is shared across poses. Check an independent load cycle at
    # the representative angles near and beyond the former 10 N first hit.
    cycles = [row for angle in (60., 64., 68., 75., 85., 96.)
              for row in evaluate(angle, gain, [0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 8., 6., 4., 2., 0.])]
    save(folder / "FJ1500_tpu_coupling_screen.json", dict(
        generator=provenance(ROOT / "tools/study_s1_infraspinatus_tpu_coupling.py"),
        solver=provenance(ROOT / "tools/study_s1_y_terminal_equilibrium.py"),
        prior_path=provenance(path_file), coupling_direction_source=direction.tolist(),
        candidate=candidate, trials=statistics, rows=rows, cycles=cycles,
        scope="Passive coupled transverse strain is lateral displacement minus gain times axial shortening."
              " A single law applies to all poses. Stiffness is provisional, and the law is not tied to"
              " a manufactured lattice. The earlier lead length lower bound is still used."
              " Prescribed GH and rigid-only checks do not measure free-drive range."))
    print("Selected comparison gain", gain, "cyclic failures", sum(bool(row["hits"]) for row in cycles), flush=True)


if __name__ == "__main__":
    main()
