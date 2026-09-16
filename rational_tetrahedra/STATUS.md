# Rational dihedral tetrahedra: prototype status

## Implemented checks

The prototype constructs the six-variable Laurent polynomial from the $4\times4$ Gram determinant of outward face normals. The expansion has 105 monomials; multiplying by 16 gives integral coefficients with constant term $-20$.

Candidates are certified exactly in $\mathbb Z[t]/(\Phi_{2N}(t))$. The code verifies the two one-parameter families of Theorem 1.8 and the 59 isolated tetrahedra of Table 3 in *Space vectors forming rational angles* by Kedlaya, Kolpakov, Poonen, and Rubinstein. It checks four positive principal cofactors and consistent signs in the Gram nullspace to exclude non-tetrahedral solutions. A five-angle search recovers the sixth angle from a quadratic equation and applies exact certification to each candidate.

## Archived benchmarks

| $N$ | Five-angle iterations | Exact Gram solutions | Ordered tetrahedron representatives | Time |
|---:|---:|---:|---:|---:|
| 12 | 25,366 | 692 | 7 | 0.12 s |
| 24 | 900,956 | 5,023 | 19 | 2.12 s |
| 30 | 2,805,895 | 11,509 | 53 | 6.18 s |

These are machine-dependent timings from the original archive. The ordered representatives are not fully quotiented by $S_4$, normal-sign changes, and Regge symmetries. They include specializations of the parameter families and cannot be directly compared with 59 isolated similarity classes.

## Angle convention

The public interface uses interior dihedral angles $\alpha$. The outward-normal Gram matrix has off-diagonal entries $\cos(\pi-\alpha)=-\cos\alpha$. The original 105-term equation is first stated for angles between unoriented lines; omitting this convention change causes the Table 3 geometry checks to fail.

## Remaining work

The prototype checks consistency of known answers and performs finite-denominator searches. A complete classification still requires the mod-2 classification of root-of-unity relations with at most 12 terms (Section 6), substitution into the 105-term equation, torsion closure for the remaining systems with at most three variables (Sections 7 and 9), and consistent deduplication under $W(D_6)$, $S_4$, sign changes, and Regge symmetries.

It also requires a no-missed-candidate error bound for numerical candidate generation analogous to Lemma 8.2, and durable certificates with an independent verifier. Direct Python search at $N=420$ is not practical under the five-loop growth estimate. Parallel candidate generation with separate exact certification is a possible extension.

The classification route reduces 105 terms to mod-2 relations with at most 12 terms, then to families with at most three parameters and torsion closure. Direct recursive resultants on all six variables are not the proposed first step. The Maple lattice and back-substitution routines may serve as later cross-checks.
