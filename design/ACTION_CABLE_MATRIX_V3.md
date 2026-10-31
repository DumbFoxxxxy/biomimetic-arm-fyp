# Action-to-cable responsibility matrix V2

| Action | Allowed | Selected | 80% | Minimum motion | Tension |
| --- | --- | --- | --- | --- | --- |
| K1 | E01, E02 | E01, E02 | PASS | PASS | PASS |
| K2 | G01, G02, G03, G04, G05, G06, P01, P02 | G02, G03 | PASS | PASS | FAIL |
| K3 | G01, G02, G03, G04, G05, G06 | G04, G05 | PASS | PASS | FAIL |
| K4 | G01, G02, G03, G04, G05, G06, P01, P02 | G01, G02 | PASS | PASS | PASS |
| K5 | G01, G02, G03, G04, G05, G06, P01, P02 | G01, G03, G04, G05 | PASS | PASS | PASS |
| K6 | S01, S02, S03, S04, S05, S06 | S01, S02, S03, S04, S05, S06 | PASS | PASS | PASS |
| K7 | G01, G02, G03, S01, S02, S03, S04, S05, S06, P01, P02 | S01, S02, S03, S04, G02, S05 | PASS | PASS | FAIL |
| SLIDE_FWD | S03, S04, S01, S02 | S03, S02 | PASS | PASS | FAIL |
| SLIDE_BACK | S03, S04, S01, S02 | S03, S04 | PASS | PASS | PASS |
| SCAP_EXT | S01, S02, S03, S04, S05, S06 | S01, S03, S04, S05 | PASS | PASS | PASS |
| SCAP_INT | S01, S02, S03, S04, S05, S06 | S04, S05, S06 | FAIL | PASS | PASS |
| PECTORAL | P01, P02, G05, G06 | P01, P02, G05 | PASS | PASS | PASS |
