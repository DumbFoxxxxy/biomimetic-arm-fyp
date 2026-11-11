"""Generate and screen a diagonal TPU cell for the FJ1500 stiffness target."""

import numpy as np
from scipy.optimize import minimize_scalar

from build_channel_first_inventory import DEST, ROOT, provenance
from build_s1_teres_origin_mount import A, export, rod, save
from study_s1_infraspinatus_cap_relief import union
from study_s1_local_force_geometry import load_reference


def leaf(start, end, width, thickness):
    axis = end - start
    length = np.linalg.norm(axis)
    axis /= length
    out = np.array([0., 0., 1.])
    transverse = np.cross(out, axis)
    frame = np.c_[axis, transverse, out, start]
    block = A.I._manifold().Manifold.cube([length + 2., width, thickness])
    return block.translate([-1., -width / 2., -thickness / 2.]).transform(frame)


def section_clearance(cell, x, gap):
    bore = rod(np.array([x - .1, 0., 0.]), np.array([x + .1, 0., 0.]), gap)
    return float(A.I.intersection_volume(cell, bore))


def main():
    folder = DEST.parent / "s1_local_muscle_load/FJ1500_cap_relief"
    model_file = folder / "FJ1500_tpu_coupling_screen.json"
    screen = load_reference(provenance(model_file))
    target_gain = screen["candidate"]["gain"]
    if target_gain <= 0.:
        raise ValueError("The diagonal cell requires a positive selected coupling")

    # Coordinates are axial shortening, desired lateral displacement, and
    # the remaining transverse direction. Mirrored pairs cancel z moments.
    axial, transverse = .2, .4
    target = np.array([[axial + transverse * target_gain**2, -transverse * target_gain],
                       [-transverse * target_gain, transverse]])
    slope = np.sqrt(target[1, 1] / target[0, 0])
    gaps = np.array([10., -10. * slope])
    weak_width = .8
    thickness = .8
    base = target[0, 0] * slope
    width_ratio = (base - target[0, 1]) / (base + target[0, 1])
    families = [(gaps.copy(), weak_width * width_ratio, 1.8),
                (gaps * np.array([1., -1.]), weak_width, 3.4)]
    terms = []
    solids = []
    manifold = A.I._manifold()
    for vector, width, elevation in families:
        for sign in (-1., 1.):
            start = np.r_[0., -.5 * vector[1], sign * elevation]
            end = np.r_[vector[0], .5 * vector[1], sign * elevation]
            solids.append(leaf(start, end, width, thickness))
            terms.append((vector, width * thickness, width, thickness))
    rails = [manifold.Manifold.cube([2.4, 17., 8.6]).translate([x - 1.2, -8.5, -4.3])
             for x in (0., gaps[0])]
    bore = rod(np.array([-2., 0., 0.]), np.array([gaps[0] + 2., 0., 0.]), 1.1)
    solid = union(solids + rails) - bore
    actual, meta = export(solid, folder / "FJ1500_diagonal_TPU_cell.obj", np.zeros(3))
    wire_overlap = float(A.I.intersection_volume(actual, bore))
    if wire_overlap > .001 or any(section_clearance(actual, x, 1.09) > .001 for x in (0., gaps[0])):
        raise RuntimeError("The internal straight channel is obstructed")

    matrix_per_modulus = np.zeros((2, 2))
    for vector, area, _, _ in terms:
        length = np.linalg.norm(vector)
        matrix_per_modulus += (area / length**3) * np.outer(vector, vector)
    effective_modulus = float(np.sum(matrix_per_modulus * target) / np.sum(matrix_per_modulus**2))
    fitted = effective_modulus * matrix_per_modulus
    if not np.allclose(fitted, target, atol=1e-10):
        raise RuntimeError("Pin-jointed cell does not match its tangent design")
    cycle = []
    for shortening in (0., .25, .5, 1., 2., 3., 4., 2., 1., 0.):
        def response(sideways):
            energy = 0.
            force = np.zeros(2)
            data = []
            for vector, area, width, depth in terms:
                deformed = vector - np.array([shortening, sideways])
                length = np.linalg.norm(vector)
                actual_length = np.linalg.norm(deformed)
                strain = (actual_length - length) / length
                link_force = effective_modulus * area * strain
                energy += .5 * effective_modulus * area * length * strain**2
                force -= link_force * deformed / actual_length
                weak_I = min(width * depth**3, depth * width**3) / 12.
                euler = np.pi**2 * effective_modulus * weak_I / length**2
                data.append(dict(axial_strain=float(strain), axial_link_force_N=float(link_force),
                                 pin_ended_buckling_reference_N=float(euler),
                                 compression_to_reference_ratio=float(max(0., -link_force / euler))))
            return energy, force, data

        optimum = minimize_scalar(lambda y: response(y)[0], bounds=(-6., 6.), method="bounded",
                                  options={"xatol": 1e-10})
        if not optimum.success:
            raise RuntimeError("Failed to solve the free lateral response")
        _, force, data = response(optimum.x)
        if abs(force[1]) > 1e-6:
            raise RuntimeError("Lateral force residual exceeds tolerance")
        cycle.append(dict(shortening_mm=shortening, free_sideways_mm=float(optimum.x),
                          axial_force_N=float(force[0]), lateral_force_residual_N=float(force[1]), links=data))
    needs = []
    for angle in sorted({row["GH_deg"] for row in screen["cycles"]}):
        rows = [row for row in screen["cycles"] if row["GH_deg"] == angle]
        first = np.array(rows[0]["cap_pin_source_mm"])
        excursion = max(np.linalg.norm(np.array(row["cap_pin_source_mm"]) - first) for row in rows)
        needs.append(dict(GH_deg=angle, max_from_unloaded_cap_travel_mm=excursion))
    save(folder / "FJ1500_diagonal_TPU_cell.json", dict(
        generator=provenance(ROOT / "tools/build_s1_infraspinatus_coupled_tpu_cell.py"),
        coupling_model=provenance(model_file), mesh=meta, mesh_bounding_box_mm=list(actual.bounding_box()),
        axial_rail_pitch_mm=gaps[0], nominal_thickness_mm=thickness, nominal_weaker_width_mm=weak_width,
        nominal_stronger_width_mm=weak_width * width_ratio, bore_radius_mm=1.1,
        channel_overlap_mm3=wire_overlap, target_stiffness_N_per_mm=target.tolist(),
        ideal_stiffness_per_modulus_mm=matrix_per_modulus.tolist(),
        fitted_ideal_modulus_MPa=effective_modulus, modulus_is_fit_parameter=True,
        pin_jointed_cycle=cycle, previously_estimated_cap_excursion=needs,
        scope="Standalone geometry coupon with two mirrored diagonal ligament families, not an installed part."
              " Truss model ignores fixed root bending, large-strain TPU, self-contact and sleeve loads."
              " Euler compression values are screening references, not tested capacity."
              " The steel wire is represented by a clear nominal neutral bore. Travel is not validated."))
    print("Mesh", meta["file"], "fitted ideal modulus", effective_modulus,
          "cycle", [(row["shortening_mm"], round(row["free_sideways_mm"], 3),
                     round(max(link["compression_to_reference_ratio"] for link in row["links"]), 2))
                    for row in cycle], flush=True)


if __name__ == "__main__":
    main()
