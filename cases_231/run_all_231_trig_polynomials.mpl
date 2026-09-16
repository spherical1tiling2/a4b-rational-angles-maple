restart:
with(LinearAlgebra):
read "a4b_rational_core.mpl":
read "trig_polynomial_engine.mpl":
read "appendix_all_3vertex_cases.txt":
read "appendix_reference_polys_by_heading.mpl":

logname := "appendix_231_trig_generation.log":
lf := fopen(logname,WRITE):
fprintf(lf,"Equation used: paper (2.6) only; Appendix OR branches are not repeated.\n"):
fprintf(lf,"case,status,kind,terms,reference_loaded\n"):
count := 0: bad := 0:
for row in appendix_all_3vertex_cases do
    id := row[1]:
    if StringTools:-Search("-OR",id) <> 0 then next end if:
    count := count+1:
    T := `a4b_solver/MakeModel`(row[2]):
    if T["status"] <> "ok" then
        bad := bad+1:
        fprintf(lf,"%s,%a,,,,\n",id,T["status"]):
    else
        M := Matrix([[add(row[3][j]*T["affine_expressions"][j],j=1..6),
                      add(row[4][j]*T["affine_expressions"][j],j=1..6)]]):
        # Rebuild the two-by-two phase coefficient matrix in (_u1,_u2).
        PM := Matrix([[coeff(M[1,1],'_u1'),coeff(M[1,2],'_u1')],
                      [coeff(M[1,1],'_u2'),coeff(M[1,2],'_u2')]]):
        if Determinant(PM)=0 then
            F := `a4b_solver/PhasePolynomialFreeCoordinates`(T):
            fprintf(lf,"%s,free_coordinate_fallback,%a,%d,D1=%a;D2=%a\n",id,T["kind"],nops([op(expand(F["polynomial"]))]),F["D1"],F["D2"]):
        else
            P := `a4b_solver/PhasePolynomial`(T,row[3],row[4]):
            fprintf(lf,"%s,ok,%a,%d,%a\n",id,T["kind"],nops([op(expand(P))]),assigned(appendix_reference_polys_by_heading[count])):
        end if:
    end if:
end do:
fprintf(lf,"SUMMARY generated=%d bad=%d\n",count,bad):
fclose(lf):
printf("SUMMARY generated=%d bad=%d\n",count,bad):
quit:
