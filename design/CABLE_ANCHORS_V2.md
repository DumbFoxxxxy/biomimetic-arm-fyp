# Cable anchors V2

Generated from the force-aware V2 simulation. This file is provisional until every staged action passes.

| Action | Selected cables | 80% | Minimum motion | Tension solution |
| --- | --- | --- | --- | --- |
| K1 | E01, E02 | PASS | PASS | PASS |
| K2 | G02, G05 | PASS | PASS | PASS |
| K3 | G01, G02, G05, G06 | PASS | PASS | PASS |
| K4 | G01, G02, G03, G05 | PASS | PASS | PASS |
| K5 | G01, G02, G03, G05 | PASS | PASS | PASS |
| K6 | S01, S02, S03, S04, S06 | PASS | PASS | PASS |
| K7 | S01, S02, S03, S04, G01, G03, S06 | PASS | PASS | PASS |
| SLIDE_FWD | S03, S04, S01 | PASS | PASS | PASS |
| SLIDE_BACK | S03, S04 | PASS | PASS | PASS |
| SCAP_EXT | S04, S05 | PASS | PASS | PASS |
| SCAP_INT | S04, S05, S06 | PASS | PASS | PASS |
| PECTORAL | P01, P02, G05 | PASS | PASS | PASS |

## Route coordinates

| Cable | Role | Element | Part | World (x,y,z) mm | Radius mm | Wrap deg |
| --- | --- | --- | --- | --- | ---: | ---: |
| E01 | elbow_flexion | 0:anchor | humerus | (162.55, 147.69, -7.70) |  |  |
| E01 | elbow_flexion | 1:pass_through | humerus | (183.03, -86.41, -12.70) |  |  |
| E01 | elbow_flexion | 2:anchor | ulna | (194.80, -220.89, -8.18) |  |  |
| E02 | elbow_extension | 0:anchor | humerus | (162.55, 147.69, 22.30) |  |  |
| E02 | elbow_extension | 1:pass_through | humerus | (183.03, -86.41, 27.30) |  |  |
| E02 | elbow_extension | 2:anchor | ulna | (194.80, -220.89, 21.82) |  |  |
| G01 | gh_flexion | 0:anchor | scapula | (160.11, 231.37, -20.09) |  |  |
| G01 | gh_flexion | 1:guide | scapula | (162.36, 186.69, -16.02) | 8.0 | 15.540207191066381 |
| G01 | gh_flexion | 2:anchor | humerus | (173.00, 28.15, -7.70) |  |  |
| G02 | gh_extension | 0:anchor | scapula | (121.37, 250.84, 42.67) |  |  |
| G02 | gh_extension | 1:guide | scapula | (152.39, 184.84, 27.36) | 8.0 | 29.248467790733837 |
| G02 | gh_extension | 2:anchor | humerus | (173.00, 28.15, 22.30) |  |  |
| G03 | gh_abduction | 0:anchor | scapula | (164.15, 239.15, -18.29) |  |  |
| G03 | gh_abduction | 1:guide | scapula | (173.12, 186.82, -0.44) | 8.0 | 27.46156144987087 |
| G03 | gh_abduction | 2:anchor | humerus | (187.95, 29.46, 7.30) |  |  |
| G04 | gh_adduction | 0:anchor | scapula | (120.43, 152.53, 54.07) |  |  |
| G04 | gh_adduction | 1:guide | scapula | (140.59, 171.87, 14.31) | 8.0 | 120.42480208088934 |
| G04 | gh_adduction | 2:anchor | humerus | (158.06, 26.84, 7.30) |  |  |
| G05 | gh_internal_rotation | 0:anchor | scapula | (108.36, 150.21, 64.55) |  |  |
| G05 | gh_internal_rotation | 1:guide | scapula | (156.03, 132.89, 4.49) | 8.0 | 136.72516047292876 |
| G05 | gh_internal_rotation | 2:anchor | humerus | (153.29, 131.82, -3.31) |  |  |
| G06 | gh_external_rotation | 0:anchor | scapula | (51.10, 112.86, 114.24) |  |  |
| G06 | gh_external_rotation | 1:guide | scapula | (154.10, 129.19, 25.55) | 8.0 | 138.44634370069102 |
| G06 | gh_external_rotation | 2:anchor | humerus | (153.29, 131.82, 17.91) |  |  |
| S01 | upper_scapula_pull_down_medial | 0:anchor | thorax | (60.21, 117.80, 103.55) |  |  |
| S01 | upper_scapula_pull_down_medial | 1:anchor | scapula | (70.98, 273.53, 83.75) |  |  |
| S02 | upper_scapula_pull_down_lateral | 0:anchor | thorax | (85.53, 173.59, 80.74) |  |  |
| S02 | upper_scapula_pull_down_lateral | 1:anchor | scapula | (120.03, 231.75, 41.18) |  |  |
| S03 | lower_scapula_pull_up_medial | 0:anchor | thorax | (74.23, 260.55, 91.39) |  |  |
| S03 | lower_scapula_pull_up_medial | 1:anchor | scapula | (62.14, 126.65, 91.42) |  |  |
| S04 | lower_scapula_pull_up_lateral | 0:anchor | thorax | (77.12, 235.77, 88.74) |  |  |
| S04 | lower_scapula_pull_up_lateral | 1:anchor | scapula | (136.25, 80.14, 27.09) |  |  |
| S05 | auxiliary_scapula_rotation_positive | 0:anchor | world | (160.00, 220.00, -30.00) |  |  |
| S05 | auxiliary_scapula_rotation_positive | 1:anchor | scapula | (114.97, 199.52, 58.81) |  |  |
| S06 | antagonistic_scapula_rotation | 0:anchor | world | (180.00, 220.00, 70.00) |  |  |
| S06 | antagonistic_scapula_rotation | 1:anchor | scapula | (110.62, 192.67, 62.59) |  |  |
| P01 | pectoral_clavicular_head | 0:anchor | clavicle_body | (90.00, 234.00, -90.70) |  |  |
| P01 | pectoral_clavicular_head | 1:guide | thorax | (116.41, 213.13, -63.07) | 8.0 | 32.58062866468013 |
| P01 | pectoral_clavicular_head | 2:anchor | humerus | (167.59, 133.08, -6.71) |  |  |
| P02 | pectoral_sternocostal_head | 0:anchor | thorax | (25.00, 200.00, -107.00) |  |  |
| P02 | pectoral_sternocostal_head | 1:guide | thorax | (114.40, 193.22, -77.17) | 8.0 | 59.95500362577589 |
| P02 | pectoral_sternocostal_head | 2:anchor | humerus | (164.77, 107.73, -7.14) |  |  |
