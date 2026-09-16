# Rational dihedral tetrahedra - independent audit prototype

This directory contains a dependency-free Python prototype for the problem
studied in *Space vectors forming rational angles* by Kedlaya, Kolpakov,
Poonen, and Rubinstein.

The prototype deliberately does not import or copy the authors' SageMath,
C++, or Magma implementation.  It currently provides four independently
checkable pieces:

1. construct the six-variable Laurent polynomial obtained from the `4 x 4`
   Gram determinant and verify that it has 105 monomials;
2. certify root-of-unity substitutions exactly in
   `Z[t]/(Phi_(2N)(t))`;
3. symbolically certify the two one-parameter tetrahedron families;
4. reproduce the low-order numerical search pattern from Section 8, followed
   by exact cyclotomic certification and the tetrahedron realizability test.

This is an audit/feasibility prototype, not yet an independent proof of the
complete classification.  A complete proof also needs the mod-2 relation
enumeration and torsion-closure stages from Sections 6, 7, and 9 of the paper.

## Run

From this directory's parent:

```powershell
python -m rational_tetrahedra audit
python -m rational_tetrahedra search 12
python -m unittest discover rational_tetrahedra/tests
```

The angle order throughout is the paper's order

```text
(alpha12, alpha34, alpha13, alpha24, alpha14, alpha23).
```

An integer tuple `(m12,m34,m13,m24,m14,m23)` at level `N` represents
angles `(pi/N)` times that tuple.

The public API and CLI use **interior dihedral angles**.  Internally, the
Gram matrix is built from outward face normals, so its off-diagonal entries
are `-cos(alpha_ij)`.  This convention is intentionally explicit because the
paper's polynomial is first stated for angles between unoriented lines.

Use `--include-all-gram-solutions` only for diagnostics.  By default the
search prints counts and the tuples that pass the nondegenerate tetrahedron
test; most exact Gram solutions are coplanar or have the wrong null-vector
orientation and are therefore not tetrahedra.
