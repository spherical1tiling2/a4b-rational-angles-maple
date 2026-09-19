restart:
with(LinearAlgebra):
read "a4b_rational_core.mpl":
read "canonical_phase_engine.mpl":
read "recursive_torsion_engine.mpl":
read "root_backsubstitution.mpl":
read "support_lattice.mpl":
read "phase_backsubstitution.mpl":
read "appendix_all_3vertex_cases.txt":

first_case := 1:
last_case := nops(appendix_all_3vertex_cases):
tol := 10^(-28):
if FileTools:-Exists("solution_audit_range.mpl") then
  read "solution_audit_range.mpl":
end if:
if last_case>nops(appendix_all_3vertex_cases) then
  last_case:=nops(appendix_all_3vertex_cases):
end if:

R26 := proc(aa,bb,cc,dd,ee)
  ((1-cos(Pi*bb))*sin(Pi*(dd-aa/2))
   -(1-cos(Pi*cc))*sin(Pi*(ee-aa/2)))*sin(Pi*(dd-ee)/2)
   -(1-cos(Pi*(bb-cc)))*sin(Pi*aa/2)*sin(Pi*(dd+ee)/2)
end proc:

caseByID := table():
for row in appendix_all_3vertex_cases do caseByID[row[1]]:=row end do:
appendixExact := []: appendixUnique := []:
candidateText := FileTools:-Text:-ReadFile("appendix_3vertex_solution_candidates.txt"):
for line in StringTools:-Split(candidateText,"\n") do
  if type(line,string) then
    candidateLine := StringTools:-Trim(line):
    if candidateLine<>"" then
      if substring(candidateLine,1..1)<>"#" then
        parts := StringTools:-Split(candidateLine,"|"):
        nums := map(parse,StringTools:-Split(parts[3],",")):
        den := parse(parts[4]): f0 := parse(parts[2]):
        av0 := [seq(nums[j]/den,j=1..5)]:
        key0 := cat(parts[1],"|",f0,"|",av0):
        keyu := cat(f0,"|",av0):
        if not member(key0,appendixExact) then appendixExact:=[op(appendixExact),key0] end if:
        if not member(keyu,appendixUnique) then appendixUnique:=[op(appendixUnique),keyu] end if:
      end if:
    end if:
  end if:
end do:

foundExact := []: foundUnique := []:
checkpointFile := "recursive_231_solution_checkpoint.txt":
if FileTools:-Exists(checkpointFile) then
  checkpointText := FileTools:-Text:-ReadFile(checkpointFile):
  for line in StringTools:-Split(checkpointText,"\n") do
    if type(line,string) then
      if StringTools:-Trim(line)<>"" then
      checkpointParts := StringTools:-Split(line,"\t"):
      if nops(checkpointParts)>=2 then
        if not member(checkpointParts[1],foundExact) then foundExact:=[op(foundExact),checkpointParts[1]] end if:
        if not member(checkpointParts[2],foundUnique) then foundUnique:=[op(foundUnique),checkpointParts[2]] end if:
      end if:
      end if:
    end if:
  end do:
end if:
if FileTools:-Exists("recursive_231_solution_audit.csv") then
  csvText := FileTools:-Text:-ReadFile("recursive_231_solution_audit.csv"):
  if type(csvText,string) then
    if StringTools:-Trim(csvText)<>"" then
      lf := fopen("recursive_231_solution_audit.csv",APPEND):
    else
      lf := fopen("recursive_231_solution_audit.csv",WRITE):
      fprintf(lf,"case,f,angles_pi,residual,status,branch\n"):
    end if:
  else
    lf := fopen("recursive_231_solution_audit.csv",WRITE):
    fprintf(lf,"case,f,angles_pi,residual,status,branch\n"):
  end if:
else
  lf := fopen("recursive_231_solution_audit.csv",WRITE):
  fprintf(lf,"case,f,angles_pi,residual,status,branch\n"):
end if:

Progress := proc(msg)
  local lp;
  lp := fopen("recursive_231_solution_progress.log",APPEND):
  fprintf(lp,"%s\n",msg):
  fclose(lp):
end proc:
Progress("START full audit"):
normCaseIDs := ["APP-019","APP-021","APP-048","APP-083","APP-084","APP-133",
                "APP-184","APP-197","APP-200","APP-202","APP-210","APP-223"]:
reduceBeforeNorm := false:

`a4b_solver/AcceptAffineCandidate` := proc(id,av,branch,lf)
  local qv,fv,res,keye,keyu,cf;
  global foundExact, foundUnique;
  if not andmap(v->evalb(v>0 and v<2),av[1..5]) then return NULL end if:
  if av[2]=av[3] or av[4]=av[5] then return NULL end if:
  qv := normal(av[6]):
  if qv<=0 then return NULL end if:
  fv := normal(1/qv):
  if denom(fv)<>1 or not type(fv,integer) or fv<12 or irem(fv,2)<>0 then return NULL end if:
  res := abs(evalf(R26(op(av)),50)):
  if res>=tol then return NULL end if:
  keye := cat(id,"|",fv,"|",av[1..5]):
  keyu := cat(fv,"|",av[1..5]):
  if not member(keye,foundExact) then
    foundExact := [op(foundExact),keye]:
    fprintf(lf,"%s,%a,%a,%a,PASS,%s\n",id,fv,av[1..5],res,branch):
    cf := fopen(checkpointFile,APPEND):
    fprintf(cf,"%s\t%s\n",keye,keyu):
    fclose(cf):
  end if:
  if not member(keyu,foundUnique) then foundUnique := [op(foundUnique),keyu] end if:
  return [fv,av[1..5]]:
end proc:

# Derive a finite, exact integer-lift box from two linearly independent
# angle coordinates.  The old fixed [-20,20] box was only a heuristic.
# If 0<a_i<2 and 0<a_j<2, solving the two affine equations for k1,k2
# bounds the lift variables by the extrema on the square [0,2]^2.
`a4b_solver/AffineLiftBounds` := proc(back)
  local b1,b2,aff,i,j,a11,a12,a21,a22,det,c1,c2,s1,s2,v1,v2,
        e1,e2,vals1,vals2,lo1,hi1,lo2,hi2,found;
  b1:=back["branch_variables"][1]: b2:=back["branch_variables"][2]:
  aff:=back["affine_angles"]: found:=false:
  for i from 1 to 5 do
    for j from i+1 to 5 do
      a11:=normal(coeff(aff[i],b1)): a12:=normal(coeff(aff[i],b2)):
      a21:=normal(coeff(aff[j],b1)): a22:=normal(coeff(aff[j],b2)):
      det:=normal(a11*a22-a12*a21):
      if det<>0 and not found then
        c1:=normal(subs(b1=0,b2=0,aff[i])):
        c2:=normal(subs(b1=0,b2=0,aff[j])):
        vals1:=[]: vals2:=[]:
        for s1 in [0,2] do
          for s2 in [0,2] do
            v1:=s1: v2:=s2:
            e1:=normal(((v1-c1)*a22-a12*(v2-c2))/det):
            e2:=normal((a11*(v2-c2)-(v1-c1)*a21)/det):
            vals1:=[op(vals1),e1]: vals2:=[op(vals2),e2]:
          end do:
        end do:
        lo1:=floor(min(op(vals1)))-1: hi1:=ceil(max(op(vals1)))+1:
        lo2:=floor(min(op(vals2)))-1: hi2:=ceil(max(op(vals2)))+1:
        found:=true:
      end if:
    end do:
  end do:
  if not found then error "the five angle coordinates do not bound two lift variables" end if:
  return table(["k1_min"=lo1,"k1_max"=hi1,"k2_min"=lo2,"k2_max"=hi2,
                "row_pair"=[i,j]]):
end proc:

`a4b_solver/EnumerateRecoveredBranches` := proc(id,T,Q,tree,N,z,lf,red)
  local cand,rawcandidates,lifted,nx,rx,ny,ry,back,bounds,k1,k2,av,branch,got,j;
  got := 0:
  rawcandidates:=`a4b_solver/Recover2DTorsionCandidates`(tree,N,z):
  # If the factor was reduced to a full exponent lattice, this is the exact
  # finite-index lift from (U,V)-torsion to the original (X,Y)-torsion.
  if assigned(red["composed_base"]) then
    lifted:=`a4b_solver/LiftReducedTorsionCandidatesComposite`(
               rawcandidates,red,red["composed_base"]):
  else
    lifted:=`a4b_solver/LiftReducedTorsionCandidates`(rawcandidates,red):
  end if:
  for cand in lifted do
    nx:=cand["x_order"]: rx:=cand["x_residue"]:
    ny:=cand["y_order"]: ry:=cand["y_residue"]:
    back := `a4b_solver/CanonicalRootBackSubstitute`(T,Q,nx,rx,ny,ry):
    bounds := `a4b_solver/AffineLiftBounds`(back):
    for k1 from bounds["k1_min"] to bounds["k1_max"] do
      for k2 from bounds["k2_min"] to bounds["k2_max"] do
        av := [seq(normal(subs(back["branch_variables"][1]=k1,
                              back["branch_variables"][2]=k2,
                              back["affine_angles"][j])),j=1..6)]:
        branch := cat(cand["branch"],";k=",k1,",",k2):
        if `a4b_solver/AcceptAffineCandidate`(id,av,branch,lf)<>NULL then got:=got+1 end if:
      end do:
    end do:
  end do:
  return got:
end proc:

treeLog := fopen("recursive_231_solution_tree.log",WRITE):
fprintf(treeLog,"case,factor,status,candidate_count,accepted_count\n"):
processed:=0: failed:=0: accepted:=0:
coveredPrimary:=0:
totalRawIndex:=0:
for totalRow in appendix_all_3vertex_cases do
  totalRawIndex:=totalRawIndex+1:
  if totalRawIndex>last_case then next end if:
  if type(totalRow,list) then
    if nops(totalRow)>0 then
      totalID:=totalRow[1]:
      if type(totalID,string) then
        if StringTools:-Search("-OR",totalID)=0 then
          coveredPrimary:=coveredPrimary+1
        end if:
      end if:
    end if:
  end if:
end do:

for idx from first_case to last_case do
  row:=appendix_all_3vertex_cases[idx]:
  if not type(row,list) then next end if:
  if nops(row)=0 then next end if:
  if StringTools:-Search("-OR",row[1])<>0 then next end if:
  processed:=processed+1: id:=row[1]:
  Progress(sprintf("CASE_BEGIN %s raw_index=%d primary=%d",id,idx,processed)):
  try
    X:='X':Y:='Y':zeta:='zeta':
    T:=`a4b_solver/MakeModel`(row[2]):
    Progress(sprintf("MODEL_DONE %s status=%s",id,T["status"])):
    if T["status"]<>"ok" then failed:=failed+1: next end if:
    Q:=`a4b_solver/CanonicalPhasePolynomial`(T,row[3],row[4]):
    Progress(sprintf("PHASE_DONE %s N=%a D1=%a D2=%a",id,Q["N"],Q["D1"],Q["D2"])):
    P:=`a4b_solver/PrimitivePolynomial`(Q["polynomial"],Q["variables"]):
    fac:=factors(P)[2]:
    Progress(sprintf("FACTORS_DONE %s count=%d",id,nops(fac))):
    factor_no:=0:
    for idxf from 1 to nops(fac) do
      Fbase:=fac[idxf][1]:
      try
        fieldfac:=`a4b_solver/CyclotomicFactorList`(Fbase,Q["variables"],Q["N"],zeta):
      catch:
        fieldfac:=[[Fbase,1]]:
        Progress(sprintf("FIELD_FACTORS_FALLBACK %s factor=%d detail=%a",id,idxf,lasterror)):
      end try:
      Progress(sprintf("FIELD_FACTORS_DONE %s factor=%d count=%d",id,idxf,nops(fieldfac))):
      for idxg from 1 to nops(fieldfac) do
        factor_no:=factor_no+1:
        F:=fieldfac[idxg][1]: active:=[]:
        for v in Q["variables"] do if has(F,v) then active:=[op(active),v] end if end do:
        Progress(sprintf("FACTOR_BEGIN %s factor=%d variables=%a",id,factor_no,active)):
        if nops(active)=2 then
          # For a factor over Q(zeta_N), take its exact coefficient norm
          # before applying the rational seven-transform theorem.  This is a
          # deliberate superset step: every torsion point of F is a torsion
          # point of Norm(F), and the final real residual check removes
          # conjugate-factor candidates that do not solve the original model.
          baseRed:=`a4b_solver/FullLatticeReduction2D`(F,active):
          reduced_before_norm:=false:
          if reduceBeforeNorm and Q["N"]<>1 and member(id,normCaseIDs) and baseRed["status"]="reduced" then
            F:=baseRed["polynomial"]: active:=baseRed["variables"]:
            reduced_before_norm:=true:
            Progress(sprintf("SUPPORT_REDUCED_BEFORE_NORM %s factor=%d index=%a",id,factor_no,baseRed["index"])):
          end if:
          solverF:=F: solverN:=Q["N"]:
          if Q["N"]<>1 and member(id,normCaseIDs) then
            try
              solverF:=`a4b_solver/PrimitivePolynomial`(
                           `a4b_solver/CyclotomicNorm`(F,Q["N"],zeta),active):
              solverN:=1:
              Progress(sprintf("NORM_DONE %s factor=%d degree=%a,%a",id,factor_no,
                               degree(solverF,active[1]),degree(solverF,active[2]))):
            catch:
              Progress(sprintf("NORM_FALLBACK %s factor=%d detail=%a",id,factor_no,lasterror)):
              solverF:=F: solverN:=Q["N"]:
            end try:
          end if:
          red:=`a4b_solver/FullLatticeReduction2D`(solverF,active):
          if reduced_before_norm and solverN=1 then red["composed_base"]:=baseRed end if:
          Progress(sprintf("LATTICE_DONE %s factor=%d status=%s index=%a rank=%a",id,factor_no,red["status"],red["index"],
                           `a4b_solver/AnalyzeSupportLattice2D`(solverF,active)["rank"])):
          if red["status"]="rank1" then
            r1red:=`a4b_solver/Rank1Reduction2D`(solverF,active):
            r1diag:=`a4b_solver/UnitRootDiagnostics`(r1red["polynomial"],r1red["variable"]):
            Progress(sprintf("RANK1_FAMILY %s factor=%d direction=%a univariate=%a root_orders=%a",id,factor_no,
                             r1red["direction"],r1red["polynomial"],r1diag["cyclotomic_orders"])):
            fprintf(treeLog,"%s,%d,rank1_family,0,0\n",id,factor_no):
            next
          end if:
          tree:=`a4b_solver/RecursiveTorsion2DGCDSplitNorm`(red["polynomial"],red["variables"],solverN,zeta):
          Progress(sprintf("TREE_DONE %s factor=%d branches=%d",id,factor_no,nops(tree["univariate_branches"]))):
          before:=nops(foundExact):
          got:=`a4b_solver/EnumerateRecoveredBranches`(id,T,Q,tree,solverN,zeta,lf,red):
          for sub in tree["common_subtrees"] do
            got:=got+`a4b_solver/EnumerateRecoveredBranches`(id,T,Q,sub,solverN,zeta,lf,red):
          end do:
          Progress(sprintf("BACKSUB_DONE %s factor=%d candidates=%d accepted=%d",id,factor_no,got,nops(foundExact)-before)):
          accepted:=accepted+(nops(foundExact)-before):
          fprintf(treeLog,"%s,%d,ok,%d,%d\n",id,factor_no,got,nops(foundExact)-before):
        elif nops(active)=1 then
          fprintf(treeLog,"%s,%d,one_variable_factor,0,0\n",id,factor_no):
        else
          fprintf(treeLog,"%s,%d,skipped,0,0\n",id,factor_no):
        end if:
      end do:
    end do:
    Progress(sprintf("CASE_DONE %s accepted_total=%d",id,nops(foundExact))):
  catch:
    failed:=failed+1:
    fprintf(treeLog,"%s,ERROR,error,0,0\n",id):
    Progress(sprintf("CASE_ERROR %s detail=%a",id,lasterror)):
  end try:
end do:
fclose(lf): fclose(treeLog):
cmp:=fopen("recursive_231_solution_comparison.md",WRITE):
fprintf(cmp,"# Independent solution comparison\n\n"):
fprintf(cmp,"- Processed primary cases: %d\n",processed):
fprintf(cmp,"- Primary cases covered through raw index %d: %d\n",last_case,coveredPrimary):
fprintf(cmp,"- Failed cases: %d\n",failed):
fprintf(cmp,"- Appendix exact rows: %d\n",nops(appendixExact)):
fprintf(cmp,"- Appendix unique solutions: %d\n",nops(appendixUnique)):
fprintf(cmp,"- Program exact solutions: %d\n",nops(foundExact)):
fprintf(cmp,"- Program unique solutions: %d\n\n",nops(foundUnique)):
fprintf(cmp,"## Appendix rows not independently recovered\n\n"):
for key in appendixExact do if not member(key,foundExact) then fprintf(cmp,"- `%s`\n",key) end if end do:
fprintf(cmp,"\n## Program rows absent from Appendix\n\n"):
for key in foundExact do if not member(key,appendixExact) then fprintf(cmp,"- `%s`\n",key) end if end do:
fprintf(cmp,"\n## Unique solutions absent from Appendix\n\n"):
for key in foundUnique do if not member(key,appendixUnique) then fprintf(cmp,"- `%s`\n",key) end if end do:
fclose(cmp):
printf("SUMMARY processed=%d failed=%d appendix=%d found=%d unique_found=%d\n",processed,failed,nops(appendixExact),nops(foundExact),nops(foundUnique)):
quit:
