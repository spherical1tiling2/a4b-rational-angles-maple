restart:
kernelopts(numcpus=1):
with(LinearAlgebra):
read "recursive_torsion_engine.mpl":
read "unique_2vertex_case.txt":
try
V := Matrix(unique_2v_case[2]):
A := Matrix(3,7,0):
for i from 1 to 2 do
    for j from 1 to 5 do A[i,j]:=V[i,j] end do:
    A[i,7]:=2:
end do:
for j from 1 to 5 do A[3,j]:=1 end do:
A[3,6]:=-4: A[3,7]:=3:
R := ReducedRowEchelonForm(A):
printf("Affine system RREF = %a\n",R):
x:='x': y:='y': z:='z':
L := x^4*y^2*z-2*x^4*y*z-x^3*y^2*z+x^2*y^2*z^2+x^4*z-x^3*y^2
     +4*x^3*y*z-x^3*z^2-x^3*z+2*x^2*y^2-6*x^2*y*z+2*x^2*z^2
     -x*y^2*z-x*y^2+4*x*y*z-x*z^2+y^2*z+x^2-x*z-2*y*z+z:
step := `a4b_solver/ResultantStep`(L,[x,y,z],z):
if step["branch_count"]<>15 then error "expected 15 resultant branches" end if:
rf:=fopen("two_vertex_resultants.mpl",WRITE):
fprintf(rf,"# Resultants in z of the manuscript polynomial and its sign/square transforms.\n"):
fprintf(rf,"two_vertex_resultants := [\n"):
for i from 1 to nops(step["branches"]) do
    b:=step["branches"][i]:
    fprintf(rf,"  [\"%s\", %a]",b["id"],factor(b["resultant"])):
    if i<nops(step["branches"]) then fprintf(rf,",\n") else fprintf(rf,"\n") end if:
end do:
fprintf(rf,"]:\n"): fclose(rf):
out:=fopen("two_vertex_candidates.csv",WRITE):
fprintf(out,"case,f,angles_pi,polynomial_exact_zero,equation_26_exact_zero,selection\n"):
checked:=0: retained:=0:
for entry in unique_2v_case[4] do
    ang:=[seq(entry[1][j]/entry[2],j=1..5)]:
    aa,bb,cc,dd,ee:=op(ang):
    qval:=normal((add(ang)-3)/4): fval:=normal(1/qval):
    if fval<>20 then error "manuscript two-vertex tile count" end if:
    for v in unique_2v_case[2] do
        if add(v[j]*ang[j],j=1..5)<>2 then error "vertex equation" end if
    end do:
    polynomialValue:=simplify(evalc(subs(x=exp(4*I*Pi/fval),y=exp(I*Pi*(dd+ee)),z=exp(I*Pi*(dd-ee)),L))):
    r26:=simplify(((1-cos(Pi*bb))*sin(Pi*(dd-aa/2))-(1-cos(Pi*cc))*sin(Pi*(ee-aa/2)))*sin(Pi*(dd-ee)/2)
                 -(1-cos(Pi*(bb-cc)))*sin(Pi*aa/2)*sin(Pi*(dd+ee)/2)):
    if polynomialValue<>0 or r26<>0 then error "candidate does not satisfy the exact equations" end if:
    if (bb-cc)*(dd-ee)>=0 then selection:="excluded_angle_order_lemma"
    elif aa>=3/2 then selection:="excluded_self_intersection_case"
    else selection:="retained_manuscript_candidate": retained:=retained+1 end if:
    checked:=checked+1:
    fprintf(out,"B2V-1,%a,\"%a\",true,true,%s\n",fval,ang,selection):
end do:
fclose(out):
if checked<>4 or retained<>1 then error "two-vertex candidate count" end if:
printf("PASS: resultant_branches=15; exact_candidates=4; retained_manuscript_candidates=1\n"):
catch:
    printf("FAIL: %a\n",lastexception):
end try:
quit:

