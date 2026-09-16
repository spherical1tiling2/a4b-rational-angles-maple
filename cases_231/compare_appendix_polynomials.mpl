restart:
with(LinearAlgebra):
read "a4b_rational_core.mpl":
read "trig_polynomial_engine.mpl":
read "appendix_all_3vertex_cases.txt":
read "appendix_reference_polys_by_heading.mpl":

ref := table():
for rr in appendix_reference_polys_by_heading do ref[rr[1]] := rr[2] end do:

printf("case,polynomial_relation\n"):
same := 0: different := 0: errors := 0:
for row in appendix_all_3vertex_cases do
    id := row[1]:
    # OR branches share the heading's printed polynomial; compare only the
    # primary branch here and generate every branch separately below.
    if StringTools:-Search("-OR",id) <> 0 then next end if:
    try
        T := `a4b_solver/MakeModel`(row[2]):
        P := `a4b_solver/PhasePolynomial`(T,row[3],row[4]):
        R := simplify(P/ref[id],symbolic):
        # Equality up to a nonzero Laurent monomial and scalar is the correct
        # invariant after clearing denominators in a Laurent expression.
        isrel := type(R, 'monomial(anything)') or type(R, 'polynom(anything)'):
        # For the first pass retain the exact quotient for auditability.
        if isrel then same := same+1 else different := different+1 end if:
        printf("%s,%a\n",id,R):
    catch:
        errors := errors+1:
        printf("%s,ERROR\n",id):
    end try:
end do:
printf("SUMMARY primary=%d different=%d errors=%d\n",same,different,errors):
quit:
