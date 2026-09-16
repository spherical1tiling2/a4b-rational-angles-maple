# The 231 primary cases

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

The three-variable demonstration is in `../three_variable_demo/`.

## Computation order

For each factor, compute the seven two-dimensional resultant branches and obtain univariate root-of-unity candidates. Substitute each retained root into the original factor and its transformation, compute their gcd, and recover the second root. Check both roots in the original factor to reject extraneous candidates introduced by coefficient norms. Affine back-substitution then recovers the independent angle variables, five angles, and $q=1/f$.

Check angle ranges, nondegeneracy, even integral $f$, and equation (2.6), then compare exact case records and unique solutions with the appendix. Fixing $f$ fixes only $q$; it does not identify the two remaining phase variables.

## Outputs

The solution-audit entry point writes `recursive_231_solution_audit.csv`, `recursive_231_solution_tree.log`, and `recursive_231_solution_comparison.md`. The CSV contains candidates satisfying equation (2.6); the log contains case, factor, and branch statistics; the comparison report records missing and additional candidates.

By default, the program covers all 231 primary cases. An optional `solution_audit_range.mpl` can override `first_case` and `last_case` for a smaller run.

## Current results

See [current_formula_solution_audit_v2.md](current_formula_solution_audit_v2.md) for the current summary. The squarefree-comparison table contains 231 passing records; the solution-audit table contains 91 algebraic candidates. The latter includes candidates excluded by subsequent combinatorial conditions, as explained in the summary. Superseded comparison reports are omitted.
