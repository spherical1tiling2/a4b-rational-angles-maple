# Computation report

Execution date: 20 September 2026 (Asia/Shanghai).

All results in this release were produced from empty working directories using the included programs. The final run applies the Euler vertex-degree bound during candidate acceptance.

## Manuscript and case coverage

Manuscript: `a4br.tex`, SHA-256 `8c173ca840a28db073dd08d224ca4583121371ae19bf82fefbfaf1883fdb11c2`.

The extractor identified 230 primary three-vertex headings and 27 alternative combinations, giving 257 distinct three-vertex records. The separate case `{alpha delta epsilon, beta^2 gamma}` is B2V-1. Thus the total is 258 cases. Input datasets were extracted from this manuscript and compared with the included files before execution.

All 257 records completed, with no reported case failures. Each record occurs exactly once in `results/case_coverage.csv`. All affine vertex equations and angle sums passed exact verification. Extraction tests checked the accepted heading punctuation, case identities, missing and duplicate headings, and preservation of outputs when extraction fails.

## Polynomial comparison

All 230 printed primary equations agree with the recomputed phase polynomials after taking squarefree parts, up to a nonzero Laurent monomial. There were zero disagreements and zero deficient primary phase bases. See `results/appendix_polynomial_squarefree_comparison.csv`.

The 27 alternative combinations are separate inputs to the recursive search. They use canonical free coordinates; a literal comparison with the parent's printed phase coordinates is not imposed on an alternative whose phase coordinates are dependent.

## Finite candidates and the Euler bound

The recursive search retained 85 case-specific finite candidate records and excluded 17 algebraic records using the Euler bound. All 17 exclusions have `f=12` and a required vertex of degree at least four. Euler's identity gives `f>=14` in these cases. The rejection is implemented in `AcceptAffineCandidate`, and the verifier independently recomputes the bound from the vertex vectors.

The APP-202 factor `x+y` gives `delta=1/2+2/f+k`, where `k` is an integer. Solving the vertex equations gives

`(alpha,beta,gamma,delta,epsilon) = (8/f,1-12/f,1-4/f,1/2+2/f+k,1/2+10/f-k)`.

Angle positivity and even `f>=12` allow `k=0` for even `f>=14`, or `k=1` for `f=14,16,18`. The three finite lifts for `k=1` are derived by the program and supply the remaining three printed finite records.

The resulting **88 case-specific finite records match all 88 extracted manuscript records**, with zero missing and zero additional admissible records. These counts include repeated angle tuples attached to different case identifiers; they are not counts of distinct tilings.

All 105 algebraic records (85 retained search records, 17 degree-bound exclusions, and three derived lifts) passed exact vertex, angle-sum, and equation (2.6) checks. For the trigonometric equation, a root-of-unity substitution reduces the numerator modulo its cyclotomic polynomial. Numerical screening in the search uses 50 digits and tolerance `10^(-28)`; the subsequent zero checks are exact.

The exclusion list is in `results/degree_bound_exclusions.csv`, all exact checks are in `results/finite_candidate_checks.csv`, and the manuscript-matching list is in `results/manuscript_finite_candidates.csv`.

## Parameter families

Both printed families satisfy the vertex equations, angle sum, and equation (2.6) identically:

- APP-011 and APP-011-OR1: `(4,4,f-2,f-2,f)/f`, even `f>=12`.
- APP-202: `(16,2f-24,2f-8,f+4,f+20)/(2*f)`, even `f>=14`.

This gives three case-specific family records in `results/parameter_families.csv`. The second APP-202 phase lift is also checked symbolically before its three finite instances are evaluated.

## Two-vertex case

The program recalculated all 15 sign/square-transform resultants of the paper's reduced three-variable polynomial. Its local variables `y,z` correspond to the paper's `u=exp(i*(delta+epsilon))` and `w=exp(i*(delta-epsilon))`; `x=exp(4*i*pi/f)`.

All four listed candidates have `f=20` and pass the exact polynomial, vertex, and equation (2.6) checks. Applying the exclusions stated in the paper rejects two by the angle-order lemma and one by the stated self-intersection case. The retained candidate is `(10,12,6,5,15)/15`. The code records these stated geometric exclusions; it does not provide a new geometric proof of them.

See `results/two_vertex_resultants.mpl` and `results/two_vertex_candidates.csv`.

## Execution and scope

Backend: Maple 2024.2, X86 64 WINDOWS, Build ID 1872373. Four independent search kernels used one Maple CPU thread each. Total driver elapsed time was 190.26 seconds. Peak aggregate working set of the search process trees was 1027461120 bytes. The memory measurement covers the search workers, not the separate verification stages. Resource thresholds, timestamps, batch timings, and input hashes are recorded under `results/`.

Factor classifications in the recursive-search logs:

- `ok`: 255
- `one_variable_factor`: 442
- `rank1_family`: 163

Rank-one and one-variable factors are recorded explicitly. The two printed parameter families are checked by the separate family program. The computations establish the reported identities and candidate comparisons; the unrestricted classification of positive-dimensional components and the geometric and combinatorial tiling arguments remain in the paper. The two-vertex candidate check starts from its four printed candidates and is not an independent exhaustive three-variable torsion search.
