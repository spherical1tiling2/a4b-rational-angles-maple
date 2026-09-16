restart:
with(LinearAlgebra):

# Four-variable computational experiment for formula (2.6).
#
# Linear constraints:
#   alpha+beta+gamma = 2,
#   alpha+beta+gamma+delta+epsilon = 3+4/f.
#
# Half-angle torus coordinates:
#   x = exp(I*Pi*beta/2),
#   y = exp(I*Pi*gamma/2),
#   z = exp(I*Pi*delta/2),
#   t = exp(I*Pi*(1/2+2/f)).
# Hence
#   exp(I*Pi*alpha/2)   = -1/(x*y),
#   exp(I*Pi*epsilon/2) = t/z.

read "recursive_torsion_engine.mpl":

x := 'x': y := 'y': z := 'z': t := 't':

# If P=exp(I*Pi*u), then C0(P)=2(1-cos(Pi*u)) and
# S0(P)=2*I*sin(Pi*u).  The common nonzero scalar in (2.6) is omitted.
C0 := P -> expand(2-P-1/P):
S0 := P -> expand(P-1/P):

Pbeta       := x^2:
Pgamma      := y^2:
PdeltaAhalf := -x*y*z^2:
PepsilonAhalf := -x*y*t^2/z^2:
PhalfDeltaMinusEpsilon := z^2/t:
PbetaMinusGamma := x^2/y^2:
PalphaHalf := -1/(x*y):
PhalfDeltaPlusEpsilon := t:

Flaurent := expand(
    (C0(Pbeta)*S0(PdeltaAhalf)
     -C0(Pgamma)*S0(PepsilonAhalf))
       *S0(PhalfDeltaMinusEpsilon)
    -C0(PbetaMinusGamma)*S0(PalphaHalf)*S0(PhalfDeltaPlusEpsilon)
):

F := `a4b_solver/PrimitivePolynomial`(Flaurent,[x,y,z,t]):
Sqf := sqrfree(F,[x,y,z,t]):
Fac := factors(F):

# All four exponent coordinates are even.  Compress the support to the
# full-size variables X=x^2, Y=y^2, Z=z^2, T=t^2.
X := 'X': Y := 'Y': Z := 'Z': T := 'T':
terms := [op(expand(F))]:
G := 0: support := []:
for term in terms do
    ex := degree(term,x): ey := degree(term,y):
    ez := degree(term,z): et := degree(term,t):
    if irem(ex,2)<>0 or irem(ey,2)<>0 or
       irem(ez,2)<>0 or irem(et,2)<>0 then
        error "the expected even-exponent reduction failed"
    end if:
    coeff0 := normal(term/(x^ex*y^ey*z^ez*t^et)):
    G := G+coeff0*X^(ex/2)*Y^(ey/2)*Z^(ez/2)*T^(et/2):
    support := [op(support),[ex/2,ey/2,ez/2,et/2]]:
end do:
G := expand(G):
Gfac := factors(G):

# The compressed support still has index two because only even powers of Z
# occur.  Put W=Z^2=exp(2*I*Pi*delta); this is the full support lattice.
W := 'W': Hpoly := 0:
for term in [op(G)] do
    eX := degree(term,X): eY := degree(term,Y):
    eZ := degree(term,Z): eT := degree(term,T):
    if irem(eZ,2)<>0 then error "the residual Z lattice was not even" end if:
    coeff1 := normal(term/(X^eX*Y^eY*Z^eZ*T^eT)):
    Hpoly := Hpoly+coeff1*X^eX*Y^eY*W^(eZ/2)*T^eT:
end do:
Hpoly := expand(Hpoly):
Hsupport := []:
for term in [op(Hpoly)] do
    Hsupport := [op(Hsupport),[degree(term,X),degree(term,Y),
                               degree(term,W),degree(term,T)]]
end do:
he0 := Hsupport[1]: hdiffs := []:
for j from 2 to nops(Hsupport) do
    hdiffs := [op(hdiffs),[seq(Hsupport[j][k]-he0[k],k=1..4)]]
end do:
HH := HermiteForm(Matrix(hdiffs)):
HHb := []:
for j from 1 to RowDimension(HH) do
    if [seq(HH[j,k],k=1..4)]<>[0,0,0,0] then
        HHb := [op(HHb),[seq(HH[j,k],k=1..4)]]
    end if:
end do:

# Lattice rank/index of the compressed support.
e0 := support[1]: diffs := []:
for j from 2 to nops(support) do
    diffs := [op(diffs),[seq(support[j][k]-e0[k],k=1..4)]]:
end do:
H := HermiteForm(Matrix(diffs)):
Hb := []:
for j from 1 to RowDimension(H) do
    if [seq(H[j,k],k=1..4)]<>[0,0,0,0] then
        Hb := [op(Hb),[seq(H[j,k],k=1..4)]]
    end if:
end do:

printf("=== FOUR-VARIABLE FORMULA (2.6) EXPERIMENT ===\n"):
printf("coordinates: x=exp(I*Pi*beta/2), y=exp(I*Pi*gamma/2), z=exp(I*Pi*delta/2)\n"):
printf("             t=exp(I*Pi*(1/2+2/f)), f even and f>=12\n"):
printf("alpha phase=-1/(x*y), epsilon phase=t/z\n"):
printf("expanded term count=%d\n",nops([op(expand(F))])):
printf("degrees [x,y,z,t]=%a\n",[degree(F,x),degree(F,y),degree(F,z),degree(F,t)]):
printf("squarefree decomposition=%a\n",Sqf):
printf("factorization=%a\n",Fac):
printf("reduced variables: X=x^2, Y=y^2, Z=z^2, T=t^2=-exp(4*I*Pi/f)\n"):
printf("reduced polynomial=%a\n",G):
printf("reduced term count=%d\n",nops([op(G)])):
printf("reduced degrees [X,Y,Z,T]=%a\n",[degree(G,X),degree(G,Y),degree(G,Z),degree(G,T)]):
printf("reduced factorization=%a\n",Gfac):
printf("support rank=%d, Hermite basis=%a\n",Rank(Matrix(diffs)),Hb):
if nops(Hb)=4 then printf("support-lattice index=%a\n",abs(Determinant(Matrix(Hb)))) end if:
printf("full-lattice variables: X=exp(I*Pi*beta), Y=exp(I*Pi*gamma), W=exp(2*I*Pi*delta), T=-exp(4*I*Pi/f)\n"):
printf("full-lattice polynomial=%a\n",Hpoly):
printf("full-lattice degrees [X,Y,W,T]=%a\n",[degree(Hpoly,X),degree(Hpoly,Y),degree(Hpoly,W),degree(Hpoly,T)]):
printf("full-lattice rank=%d, Hermite basis=%a\n",Rank(Matrix(hdiffs)),HHb):
if nops(HHb)=4 then printf("full-lattice index=%a\n",abs(Determinant(Matrix(HHb)))) end if:

quit:
