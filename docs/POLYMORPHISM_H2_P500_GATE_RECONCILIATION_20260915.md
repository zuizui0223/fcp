# P500 H2 prospective gate reconciliation

Date: 2026-09-15 JST
Status: pre-outcome, no P500 image/colour/H2 opening authorized by this document

## Finding

Two already-frozen prospective documents encode different H2 vector-support floors.

1. `docs/POLYMORPHISM_H2_PROSPECTIVE_U100_SELECTION_PROTOCOL_20260913.md` froze the prospective primary decision as requiring **at least 100 P500 species** to survive H1 plus the primary continuous-mode gate. It explicitly states that the N>=100 floor is an identifiability gate and that it must not be lowered post hoc.
2. `docs/POLYMORPHISM_H2_P500_PROSPECTIVE_MEASUREMENT_PROTOCOL_20260915.md` later encoded **20 primary H2 vectors** as the minimum for evaluability, and the corresponding authorization JSON also records 20.

The later document does not explicitly retire or justify relaxation of the earlier N>=100 confirmatory gate. Treating N>=20 as sufficient for the final prospective confirmation would therefore create an avoidable ambiguity about whether a preregistered decision floor was weakened.

## Conservative resolution before outcome opening

The earlier, stricter preregistered gate controls the confirmatory decision:

- **N < 20 H2 vectors:** computationally/structurally not evaluable under the later implementation floor.
- **20 <= N < 100:** H2 may be computed for audit/descriptive purposes, but the prospective confirmatory verdict is **UNDERIDENTIFIED**. A p-value in this band cannot produce `H2_PROSPECTIVE_WHITE_AXIS_CONFIRMED` or `...NOT_CONFIRMED`.
- **N >= 100:** the primary 10% fixed-q_white structured-null test is confirmatorily evaluable; p<0.05 supports the white-axis target and p>=0.05 does not support it.

The 20% threshold remains a sensitivity analysis and cannot rescue the primary result.

This reconciliation changes no species, image, palette, q_white loading, W statistic, structured null, random seed, acquisition rule, measurement-support gate, or outcome. It only resolves a conflict between two pre-outcome support-floor statements by retaining the stricter earlier commitment.

## Chronology boundary

This document makes no retrospective claim that historical non-access has been independently proven. The separate execution-epoch receipt establishes a new prospective chronology from its durable Git freeze forward. No P500 image request, decoded image, flower-colour outcome, D value, H2 vector, W statistic, or H2 p-value is opened by this reconciliation.

## Decision labels for the new execution epoch

The execution layer must preserve these states:

- `H2_PROSPECTIVE_NOT_EVALUABLE_MEASUREMENT_SUPPORT` when fewer than 250 species have >=40 classifiable photos;
- `H2_PROSPECTIVE_NOT_EVALUABLE_VECTOR_SUPPORT` when fewer than 20 primary H2 vectors exist;
- `H2_PROSPECTIVE_WHITE_AXIS_UNDERIDENTIFIED` when 20--99 primary H2 vectors exist;
- `H2_PROSPECTIVE_WHITE_AXIS_CONFIRMED` when N>=100 and primary p<0.05;
- `H2_PROSPECTIVE_WHITE_AXIS_NOT_CONFIRMED` when N>=100 and primary p>=0.05.

No threshold change or replacement is permitted after the execution epoch advances beyond `PREOPEN_FROZEN`.
