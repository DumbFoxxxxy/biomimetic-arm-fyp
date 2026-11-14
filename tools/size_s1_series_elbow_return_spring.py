"""Size a helical extension spring for the checked FJ1477 mounting lane."""

import math

import numpy as np

from audit_s1_series_elbow_timed_muscle_geometry import checked_trace
from build_channel_first_inventory import DEST, ROOT, provenance, read
from build_s1_teres_origin_mount import save
from simulate_elbow_series_slider import SLIDE_DAMPING, SLIDE_STIFFNESS
from study_s1_local_force_geometry import load_reference


def main():
    folder = DEST.parent / "integration_review/s1_series_elbow"
    audit_path = folder / "interface_relief_elbow_muscle_geometry.json"
    audit = read(audit_path)
    extended = load_reference(audit["extended_drive"])
    budget = load_reference(extended["budget"])
    assembly = load_reference(budget["assembly"])
    mount = load_reference(assembly["mounts"])["spring_mount"]
    trace = checked_trace(audit["trace"])
    eye_spacing = float(np.linalg.norm(np.array(mount["moving_eye_neutral_source_mm"])
                                       - mount["fixed_eye_source_mm"]))
    travel = extended["rail_stroke_mm"]
    slide = np.array([float(row["slide_mm"]) for row in trace])
    time = np.array([float(row["time_s"]) for row in trace])
    # Lengths match the checked capsule envelope used by spring_envelope.
    # Reserve 0.25 mm radial clearance inside the nominal 6.5 mm lane.
    coil_length = eye_spacing - 8.
    max_od = 6.0
    preload = .25
    nominal_rate = SLIDE_STIFFNESS / 1000.
    force_end = preload + nominal_rate * travel
    fitted_velocity = np.gradient(slide / 1000., time)
    apparent_peak = np.maximum(0., preload + nominal_rate * slide + SLIDE_DAMPING * fitted_velocity).max()
    moduli = [72000., 79000., 82000.]
    rows = []
    for wire in (.40, .45, .50, .55, .60):
        for outer in (4.5, 5., 5.5, 6.):
            mean = outer - wire
            index = mean / wire
            correction = (4. * index - 1.) / (4. * index - 4.) + .615 / index
            for turns in range(6, math.floor(coil_length / wire) - 1):
                length = (turns + 2) * wire
                rates = np.array(moduli) * wire ** 4 / (8. * mean ** 3 * turns)
                if abs(rates[1] / nominal_rate - 1.) > .05:
                    continue
                peak_stress = correction * 8. * (preload + rates[-1] * travel) * mean / (np.pi * wire ** 3)
                rows.append(dict(wire_diameter_mm=wire, coil_OD_mm=outer, active_turns=turns,
                                 closed_body_mm=length, neutral_body_room_mm=coil_length - length,
                                 spring_index=index, rate_N_per_mm=rates.tolist(),
                                 full_stroke_force_N=(preload + rates * travel).tolist(),
                                 estimated_body_shear_at_full_stroke_MPa=float(peak_stress)))
    rows.sort(key=lambda row: row["estimated_body_shear_at_full_stroke_MPa"])
    result = dict(generator=provenance(ROOT / "tools/size_s1_series_elbow_return_spring.py"),
                  geometry=provenance(audit_path), spring_mount=assembly["mounts"],
                  eye_spacing_rest_mm=eye_spacing, eye_spacing_at_measured_trace_peak_mm=eye_spacing + float(slide.max()),
                  eye_spacing_at_full_rail_stroke_mm=eye_spacing + travel,
                  full_rail_stroke_mm=travel, spring_body_length_available_at_rest_mm=coil_length,
                  nominal_lane_OD_mm=6.5, screened_body_OD_limit_mm=max_od,
                  pin_diameter_mm=1.8, nominal_hook_axial_allowance_mm_each=4.,
                  installed_preload_N=preload, target_rate_N_per_mm=nominal_rate,
                  target_static_force_full_rail_stroke_N=force_end,
                  effective_model_force_peak_N_with_5ms_difference_velocity=float(apparent_peak),
                  assumed_shear_modulus_MPa=moduli, rates_within_5_percent_at_central_modulus=rows,
                  scope="First-pass closed-coil rate and Wahl-corrected body shear estimates at full rail travel."
                        " Two inactive turns assumed. Hook shape and stress, material grade, installed preload,"
                        " fatigue, damping source and actual fit on the pins require separate verification.")
    save(folder / "return_spring_size_screen.json", result)
    print("Geometric rate matches:", len(rows), "lowest body-stress case:", rows[0] if rows else None, flush=True)


if __name__ == "__main__":
    main()
