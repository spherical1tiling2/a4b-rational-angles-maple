restart:
with(LinearAlgebra):
with(PolynomialTools):
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

`a4b_solver/SquareFreeProduct` := proc(P)
    local sf,prod,item,f;
    sf := sqrfree(expand(P)):
    prod := 1:
    for item in sf[2] do
        f := item[1]:
        prod := prod*f
    end do:
    return factor(prod)
end proc:

ref := table():
for rr in appendix_reference_polys_by_heading do ref[rr[1]]:=rr[2] end do:


lf := fopen("appendix_polynomial_squarefree_comparison.csv",WRITE):
fprintf(lf,"case,status,ratio,reason\n"):
same:=0: shared:=0: different:=0: deficient:=0:
for row in appendix_all_3vertex_cases do
    id:=row[1]:
    if StringTools:-Search("-OR",id)<>0 then next end if:
    T:=`a4b_solver/MakeModel`(row[2]):
    if T["status"]<>"ok" then fprintf(lf,"%s,linear_failure,,%s\n",id,T["status"]): next end if:
    PM:=Matrix([[coeff(add(row[3][j]*T["affine_expressions"][j],j=1..6),'_u1'),coeff(add(row[4][j]*T["affine_expressions"][j],j=1..6),'_u1')],[coeff(add(row[3][j]*T["affine_expressions"][j],j=1..6),'_u2'),coeff(add(row[4][j]*T["affine_expressions"][j],j=1..6),'_u2')]]):
    if Determinant(PM)=0 then deficient:=deficient+1: fprintf(lf,"%s,phase_rank_deficient,,printed_x_y_dependent\n",id): next end if:
    P:=`a4b_solver/PhasePolynomial`(T,row[3],row[4]):
    Gsf:=`a4b_solver/SquareFreeProduct`(numer(P)):
    Rsf:=`a4b_solver/SquareFreeProduct`(ref[id]):
    ratio:=factor(simplify(Gsf/Rsf,symbolic)):
    if `a4b_solver/IsLaurentMonomial`(ratio) then
        same:=same+1: fprintf(lf,"%s,same_squarefree_up_to_Laurent_factor,\"%a\",\n",id,ratio):
    else
        Gp:=expand(numer(P)): Rp:=expand(ref[id]): gg:=factor(gcd(Gp,Rp)):
        if not `a4b_solver/IsLaurentMonomial`(gg) then shared:=shared+1: fprintf(lf,"%s,shared_squarefree_factor_but_different_core,,common_factor_not_enough\n",id)
        else different:=different+1: fprintf(lf,"%s,different_squarefree_core,,no_nonconstant_common_factor\n",id) end if:
    end if:
end do:
fclose(lf):
printf("SUMMARY squarefree_same=%d shared=%d different=%d phase_deficient=%d\n",same,shared,different,deficient):
quit:
