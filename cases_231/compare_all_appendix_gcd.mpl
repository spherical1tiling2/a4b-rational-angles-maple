restart:
with(LinearAlgebra):
read "a4b_rational_core.mpl":
read "trig_polynomial_engine.mpl":
read "appendix_all_3vertex_cases.txt":
read "appendix_reference_polys_by_heading.mpl":

ref := table():
for rr in appendix_reference_polys_by_heading do ref[rr[1]] := rr[2] end do:

`a4b_solver/SingleTerm` := proc(q)
    local z;
    z := simplify(expand(q),symbolic):
    if z=0 then return true end if:
    if type(z,`+`) then return false end if:
    return true
end proc:

lf := fopen("appendix_polynomial_gcd_comparison.csv",WRITE):
fprintf(lf,"case,status,gcd_terms,gen_residual_terms,ref_residual_terms,reason\n"):
same := 0: extra := 0: shared := 0: different := 0: deficient := 0:
for row in appendix_all_3vertex_cases do
    id := row[1]:
    if StringTools:-Search("-OR",id) <> 0 then next end if:
    T := `a4b_solver/MakeModel`(row[2]):
    if T["status"] <> "ok" then
        fprintf(lf,"%s,linear_failure,,,,%a\n",id,T["status"]):
        next
    end if:
    PM := Matrix([[coeff(add(row[3][j]*T["affine_expressions"][j],j=1..6),'_u1'),
                   coeff(add(row[4][j]*T["affine_expressions"][j],j=1..6),'_u1')],
                  [coeff(add(row[3][j]*T["affine_expressions"][j],j=1..6),'_u2'),
                   coeff(add(row[4][j]*T["affine_expressions"][j],j=1..6),'_u2')]]):
    if Determinant(PM)=0 then
        deficient := deficient+1:
        fprintf(lf,"%s,phase_rank_deficient,,,,printed_x_y_dependent\n",id):
        next
    end if:
    P := `a4b_solver/PhasePolynomial`(T,row[3],row[4]):
    Gp := expand(numer(P)):
    Rp := expand(ref[id]):
    gg := gcd(Gp,Rp):
    A := normal(Gp/gg): B := normal(Rp/gg):
    at := nops([op(expand(A))]): bt := nops([op(expand(B))]): dt := nops([op(expand(gg))]):
    if `a4b_solver/SingleTerm`(A) and `a4b_solver/SingleTerm`(B) then
        same := same+1:
        fprintf(lf,"%s,same_up_to_nonzero_Laurent_factor,%d,%d,%d,\n",id,dt,at,bt):
    elif dt > 1 then
        shared := shared+1:
        fprintf(lf,"%s,shared_factor_but_different_core,%d,%d,%d,common_factor_not_enough\n",id,dt,at,bt):
    else
        different := different+1:
        fprintf(lf,"%s,different_core,%d,%d,%d,no_nonconstant_common_factor\n",id,dt,at,bt):
    end if:
end do:
fclose(lf):
printf("SUMMARY same=%d shared_factor=%d different=%d phase_deficient=%d\n",same,shared,different,deficient):
quit:
