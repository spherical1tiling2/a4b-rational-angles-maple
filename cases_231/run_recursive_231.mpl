restart:
with(LinearAlgebra):
read "a4b_rational_core.mpl":
read "canonical_phase_engine.mpl":
read "recursive_torsion_engine.mpl":
read "phase_backsubstitution.mpl":
read "appendix_all_3vertex_cases.txt":

# Recursive exact-elimination audit for the 230 primary cases; legacy filename.
# Each factor is handled separately.  The coefficient-field norm is taken
# only after each resultant, avoiding degree blow-up before elimination.
first_case := 1:
last_case := nops(appendix_all_3vertex_cases):

logname := "recursive_231_tree.log":
lf := fopen(logname,WRITE):
mapname := "recursive_231_phase_maps.log":
mf := fopen(mapname,WRITE):
fprintf(lf,"Recursive torsion-elimination tree\n"):
fprintf(lf,"No bounded f or denominator grid is used in this stage.\n"):
fprintf(lf,"case,factor,N,D1,D2,degreeX,degreeY,branches,nonzero_branches,cyclotomic_branch_count\n"):

processed := 0: failed := 0:
for idx from first_case to last_case do
  row := appendix_all_3vertex_cases[idx]:
  if StringTools:-Search("-OR",row[1])<>0 then next end if:
  processed := processed+1:
  id := row[1]:
  try
    X:='X': Y:='Y': zeta:='zeta':
    T := `a4b_solver/MakeModel`(row[2]):
    if T["status"]<>"ok" then
      failed := failed+1:
      fprintf(lf,"%s,MODEL_FAILURE,,,,,,,%s\n",id,T["status"]):
      next
    end if:
    Q := `a4b_solver/CanonicalPhasePolynomial`(T,row[3],row[4]):
    fprintf(mf,"%s kind=%a affine=%a phase_rows=%a canonical=(D1=%a,D2=%a,N=%a)\n",
      id,T["kind"],T["affine_expressions"],Q["phase_coordinates"],Q["D1"],Q["D2"],Q["N"]):
    P := `a4b_solver/PrimitivePolynomial`(Q["polynomial"],Q["variables"]):
    fac := factors(P)[2]:
    for idxf from 1 to nops(fac) do
      F := fac[idxf][1]:
      active := []:
      for v in Q["variables"] do
        if has(F,v) then active := [op(active),v] end if:
      end do:
      if nops(active)=1 then
        diag1 := `a4b_solver/UnitRootDiagnostics`(F,active[1]):
        fprintf(lf,"%s,%d,%a,%a,%a,%a,%a,1,1,%d\n",id,idxf,Q["N"],Q["D1"],Q["D2"],
          degree(F,active[1]),0,nops(diag1["cyclotomic_orders"])):
      elif nops(active)=2 then
        tree := `a4b_solver/RecursiveTorsion2DNormed`(F,active,Q["N"],zeta):
        nz := 0: cyc := 0:
        for b in tree["univariate_branches"] do
          if b["resultant"]<>0 then nz:=nz+1 end if:
          if nops(b["diagnostics"]["cyclotomic_orders"])>0 then cyc:=cyc+1 end if
        end do:
        fprintf(lf,"%s,%d,%a,%a,%a,%a,%a,%a,%a,%a\n",id,idxf,Q["N"],Q["D1"],Q["D2"],
          degree(F,active[1]),degree(F,active[2]),tree["step"]["branch_count"],nz,cyc):
      else
        fprintf(lf,"%s,%d,%a,%a,%a,%a,%a,SKIP,SKIP,SKIP\n",id,idxf,Q["N"],Q["D1"],Q["D2"],0,0):
      end if:
    end do:
  catch:
    failed := failed+1:
    fprintf(lf,"%s,ERROR,,,,,,,,\n",id):
  end try:
end do:

fprintf(lf,"SUMMARY processed=%d failed=%d\n",processed,failed):
fclose(lf):
fclose(mf):
printf("SUMMARY processed=%d failed=%d\n",processed,failed):
quit:
