restart:
with(LinearAlgebra):
read "a4b_rational_core.mpl":
read "trig_polynomial_engine.mpl":
read "appendix_all_3vertex_cases.txt":
read "appendix_reference_polys_by_heading.mpl":

`a4b_solver/IsMonomialExpression` := proc(e)
    local z,i;
    z := factor(e):
    if type(z,`+`) then return false end if:
    if type(z,`*`) then
        for i from 1 to nops(z) do
            if not `a4b_solver/IsMonomialExpression`(op(i,z)) then return false end if
        end do:
        return true
    end if:
    if type(z,`^`) then return `a4b_solver/IsMonomialExpression`(op(1,z)) end if:
    return true
end proc:

`a4b_solver/IsLaurentMonomial` := proc(q)
    local z;
    z := factor(simplify(q,symbolic)):
    if numer(z)=0 then return false end if:
    return `a4b_solver/IsMonomialExpression`(numer(z)) and `a4b_solver/IsMonomialExpression`(denom(z))
end proc:

ref := table():
for rr in appendix_reference_polys_by_heading do ref[rr[1]] := rr[2] end do:

lf := fopen("appendix_polynomial_gcd_comparison_v2.csv",WRITE):
fprintf(lf,"case,status,gcd_terms,ratio,reason\n"):
same := 0: shared := 0: different := 0: deficient := 0:
for row in appendix_all_3vertex_cases do
    id := row[1]:
    if StringTools:-Search("-OR",id) <> 0 then next end if:
    T := `a4b_solver/MakeModel`(row[2]):
    if T["status"] <> "ok" then
        fprintf(lf,"%s,linear_failure,,,""%s""\n",id,T["status"]): next
    end if:
    PM := Matrix([[coeff(add(row[3][j]*T["affine_expressions"][j],j=1..6),'_u1'),
                   coeff(add(row[4][j]*T["affine_expressions"][j],j=1..6),'_u1')],
                  [coeff(add(row[3][j]*T["affine_expressions"][j],j=1..6),'_u2'),
                   coeff(add(row[4][j]*T["affine_expressions"][j],j=1..6),'_u2')]]):
    if Determinant(PM)=0 then
        deficient := deficient+1:
        fprintf(lf,"%s,phase_rank_deficient,,,printed_x_y_dependent\n",id): next
    end if:
    P := `a4b_solver/PhasePolynomial`(T,row[3],row[4]):
    RpRaw := ref[id]:
    ratio := simplify(P/RpRaw,symbolic):
    if `a4b_solver/IsLaurentMonomial`(ratio) then
        same := same+1:
        fprintf(lf,"%s,same_up_to_nonzero_Laurent_factor,,""%a"",\n",id,ratio):
        next
    end if:
    Gp := expand(numer(P)):
    Rp := expand(RpRaw):
    gg := gcd(Gp,Rp):
    gt := nops([op(expand(gg))]):
    if not `a4b_solver/IsLaurentMonomial`(gg) then
        shared := shared+1:
        fprintf(lf,"%s,shared_factor_but_different_core,%d,,common_factor_not_enough\n",id,gt):
    else
        different := different+1:
        fprintf(lf,"%s,different_core,%d,,no_nonconstant_common_factor\n",id,gt):
    end if:
end do:
fclose(lf):
printf("SUMMARY same=%d shared_factor=%d different=%d phase_deficient=%d\n",same,shared,different,deficient):
quit:
