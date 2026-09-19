# Formula and solution audit: version 2

> Historical computation report: this run predates the case-count correction. The active manifest has 230 primary records and 27 OR branches, plus a separate two-vertex case. APP-158 was an exact duplicate of APP-133 and is retired; APP-177-OR1 was added. The original CSV files retain their original IDs and counts. Read archived APP-158 results as APP-133. This is not a fresh full recursive run of the corrected manifest. See [CASE_COUNT_CORRECTION.md](CASE_COUNT_CORRECTION.md).

This report summarizes the current supplied results. It concerns the 231 primary three-variable cases. The four-variable, five-variable, and six-angle experiments are separate from the manuscript classification.

## Formulas

The archived Maple check reports symbolic difference zero between its trigonometric identity and manuscript equation (2.6). All 231 active cases satisfy their three vertex equations and

\[
\alpha+\beta+\gamma+\delta+\varepsilon=3+\frac4f.
\]

After removal of nonzero scalar factors, invertible Laurent monomials, and repeated factors, all 231 appendix polynomials are reported squarefree-equivalent. The final comparison has no common-factor-only or genuinely different cases. APP-048 was resolved by correcting an extractor that missed multiplication of adjacent parentheses. APP-119 uses the Maple-generated $\zeta$ expression.

## Solutions

Parsing active `Case` blocks, excluding struck-out rows, `%%` comments, and later experimental `Subcase` content, gives 77 appendix solution records. The full recursive run processed 231 cases with zero failures and initially produced 90 accepted records.

The earlier coefficient-field candidate tree missed the root $x=-1$ in APP-084. Repeating that case with a coefficient-field norm followed by back-substitution recovered

\[
f=20,\qquad(\alpha,\beta,\gamma,\delta,\varepsilon)
=\left(\frac25,\frac45,\frac25,1,\frac35\right),
\]

with reported residual approximately $3\times10^{-49}$. Combining the full run and this correction gives 91 program records, including all 77 active appendix records. No active appendix record is missing.

## Additional algebraic candidates

The additional 14 records represent seven distinct angle tuples at $f=12$, repeated across cases. They fail equation (2.1): the assumed vertices have degree at least four, and sometimes five, requiring $f\ge14$ or $f\ge16$ respectively.

| $f$ | $(\alpha,\beta,\gamma,\delta,\varepsilon)$, in units of $\pi$ |
|---:|---|
| 12 | $(1/2,1,1/3,1,1/2)$ |
| 12 | $(1,1/2,2/3,2/3,1/2)$ |
| 12 | $(5/6,5/6,1/3,5/6,1/2)$ |
| 12 | $(1/3,2/3,1,1/3,1)$ |
| 12 | $(1/3,2/3,1,1/2,5/6)$ |
| 12 | $(2/3,4/9,8/9,1,1/3)$ |
| 12 | $(1/2,1/3,1,1/2,1)$ |

The former APP-184 row at $f=28$ and the former representative values of the APP-202 infinite family remain struck out in the manuscript and are excluded from the active solution set. APP-231 is marked `No rational solution`; later experimental subcases must not be attributed to it.

Thus the archive reports agreement for all 231 equations and recovery of all 77 active appendix solution records. The seven additional angle tuples are excluded candidates, not additional tilings.
