restart:
with(LinearAlgebra):
read "unique_2vertex_case.txt":

printf("=== UNIQUE TWO-VERTEX CASE ===\n"):
printf("case = %s\n",unique_2v_case[1]):
printf("vertices = %a\n",unique_2v_case[2]):
printf("phase basis = %s\n",unique_2v_case[3]):
printf("equation used: paper (2.6) only; no extra factor from (2.3).\n\n"):

# Generic affine solve for two vertex rows plus the pentagon total-angle row.
V := Matrix(unique_2v_case[2]):
n := RowDimension(V):
A := Matrix(n+1,6,0): Aug := Matrix(n+1,7,0):
for i from 1 to n do
  for j from 1 to 5 do A[i,j] := V[i,j] end do:
  Aug[i,7] := 2:
end do:
for j from 1 to 5 do A[n+1,j] := 1 end do:
A[n+1,6] := -4: Aug[n+1,7] := 3:
for i from 1 to n+1 do
  for j from 1 to 6 do Aug[i,j] := A[i,j] end do:
end do:

R := ReducedRowEchelonForm(Aug):
pivotcols := []: prow := table():
for i from 1 to n+1 do
  p := 0:
  for j from 1 to 6 do
    if R[i,j] <> 0 then p := j: break end if:
  end do:
  if p <> 0 then pivotcols := [op(pivotcols),p]: prow[p] := i end if:
end do:
freecols := []:
for j from 1 to 6 do
  if not member(j,pivotcols) then freecols := [op(freecols),j] end if:
end do:

params := [seq(cat(_u,k),k=1..nops(freecols))]:
exprs := Vector(6,0):
for k from 1 to nops(freecols) do exprs[freecols[k]] := params[k] end do:
for k from 1 to nops(pivotcols) do
  p := pivotcols[k]: i := prow[p]:
  exprs[p] := normal(R[i,7]-add(R[i,freecols[j]]*params[j],j=1..nops(freecols))):
end do:

printf("RREF = %a\n",R):
printf("free columns = %a\n",freecols):
printf("affine expressions [a,b,c,d,e,q] = %a\n",[seq(exprs[j],j=1..6)]):

# The second factor of the paper's three-variable polynomial after u=yz,w=y/z.
x := 'x': u := 'u': w := 'w':
L := u^2*w*x^4+u^2*w^2*x^2-u^2*w*x^3-2*u*w*x^4-u^2*x^3
     +4*u*w*x^3-w^2*x^3+w*x^4-u^2*w*x+2*u^2*x^2-6*u*w*x^2
     +2*w^2*x^2-w*x^3+u^2*w-u^2*x+4*u*w*x-w^2*x-2*u*w-w*x+x^2+w:

printf("\nCandidate verification (the four rational candidates listed in the paper):\n"):
valid := []:
for entry in unique_2v_case[4] do
  nums := entry[1]: den := entry[2]:
  ang := [seq(nums[j]/den,j=1..5)]:
  qval := normal((add(ang)-3)/4): fval := normal(1/qval):
  xv := exp(4*I*Pi/fval):
  uv := exp(I*Pi*(ang[4]+ang[5])):
  wv := exp(I*Pi*(ang[4]-ang[5])):
  residual := evalf[30](evalc(subs({x=xv,u=uv,w=wv},L))):
  # The sign relation from Lemma 2 is recorded diagnostically only.  It is
  # not imposed as an input filter in this algebraic batch.
  sign_relation := evalb((ang[2]-ang[3])*(ang[4]-ang[5]) < 0):
  non_symmetric := evalb(ang[2]<>ang[3] and ang[4]<>ang[5]):
  alpha_ok := evalb(ang[1] < 3/2):
  simple_ok := non_symmetric and alpha_ok:
  printf("angles = %a*pi, f = %a, residual(L) = %a\n",ang,fval,residual):
  printf("sign_relation_diagnostic = %a, non_symmetric = %a, alpha_bound_ok = %a, simple_candidate = %a\n",
         sign_relation,non_symmetric,alpha_ok,simple_ok):
  if simple_ok then valid := [op(valid),[ang,fval]] end if:
end do:

if nops(valid)=0 then
  printf("FINAL RESULT = none\n"):
else
  printf("FINAL VALID FIVE-ANGLE SOLUTIONS = %a\n",valid):
end if:

quit:
