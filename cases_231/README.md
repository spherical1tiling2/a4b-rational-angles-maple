# The 257 three-vertex combinations and one two-vertex case

The manifest contains **230 primary records + 27 OR branches = 257 distinct three-vertex combinations**. The separate `unique_2vertex_case.txt` makes **258 cases in total**. Names containing `231` are retained for compatibility. See [CASE_COUNT_CORRECTION.md](CASE_COUNT_CORRECTION.md).

## Programs

| File | Purpose |
|---|---|
| `a4b_rational_core.mpl` | Three-vertex linear models and affine elimination |
| `canonical_phase_engine.mpl` | Integral-exponent phase coordinates |
| `recursive_torsion_engine.mpl` | Factorization, 7/15 sign and square transformations, recursive resultants, cyclotomic diagnostics, and coefficient norms |
| `root_backsubstitution.mpl` | Recovery of eliminated roots and checks in the original factors |
| `phase_backsubstitution.mpl` | Recovery of six normalized variables and two integer phase lifts |
| `run_recursive_231_solution_audit.mpl` | Solution search and comparison with the appendix |
| `run_recursive_231.mpl` | Recursive trees and phase maps |
| `run_unique_2vertex_case.mpl` | The separate two-vertex case |
| `validate_case_manifest.mpl` | Exact count, identity, affine-model, and added-branch checks |

The three-variable demonstration is in `../three_variable_demo/`.

## Computation order

For each factor, compute the seven two-dimensional resultant branches and obtain univariate root-of-unity candidates. Substitute each retained root into the original factor and its transformation, compute their gcd, and recover the second root. Check both roots in the original factor to reject extraneous candidates introduced by coefficient norms. Affine back-substitution then recovers the independent angle variables, five angles, and $q=1/f$.

Check angle ranges, nondegeneracy, even integral $f$, and equation (2.6), then compare exact case records and unique solutions with the appendix. Fixing $f$ fixes only $q$; it does not identify the two remaining phase variables.

## Outputs

The solution-audit entry point writes `recursive_231_solution_audit.csv`, `recursive_231_solution_tree.log`, and `recursive_231_solution_comparison.md`. The CSV contains candidates satisfying equation (2.6); the log contains case, factor, and branch statistics; the comparison report records missing and additional candidates.

By default, the program covers all 230 primary cases and skips the labelled OR records. An optional `solution_audit_range.mpl` can override `first_case` and `last_case`; these are positions in the expanded manifest, not APP identifiers. Use a fresh output directory after updating: raw positions have changed, so do not resume previous checkpoints or generated reports across manifest versions.

The PowerShell extractors accept both `Case` and `Case.` and preserve APP identifiers by matching vertex combinations against `-IdReferencePath` (the checked-in manifest by default). They skip the separately managed two-vertex case and refuse to overwrite output when extraction has missing or duplicate cases. Heading order no longer determines identity. Extracted phase bases and polynomials reflect the supplied manuscript and still require mathematical comparison before replacing verified algebraic data.

## Current results

See [current_formula_solution_audit_v2.md](current_formula_solution_audit_v2.md) for the archived summary. The historical squarefree-comparison table contains 231 records, including the retired duplicate APP-158; the historical solution table contains 91 algebraic candidates. They retain their original IDs. Interpret old APP-158 results as APP-133. These tables have not been relabelled as a fresh computation. See [CASE_COUNT_CORRECTION.md](CASE_COUNT_CORRECTION.md) for this correction and its targeted validation.
