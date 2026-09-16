restart:
with(PolynomialTools):

read "recursive_torsion_engine.mpl":
read "formula26_four_variable_core.mpl":

X := 'X': Y := 'Y': W := 'W': T := 'T':
F := `a4b_solver/Formula26FourVariablePolynomial`(X,Y,W,T):
vars := [T,X,Y,W]:
transforms := `a4b_solver/SignSquareTransforms`(F,vars):

printf("=== FORMULA (2.6): FOUR TO THREE VARIABLES ===\n"):
printf("input terms=%d, degrees [T,X,Y,W]=%a\n",
       nops([op(F)]),[degree(F,T),degree(F,X),degree(F,Y),degree(F,W)]):
printf("transform count=%d (expected 31)\n",nops(transforms)):

nontrivial_gcd := 0: zero_resultants := 0: nonzero_resultants := 0:
branch_no := 0:
for tr in transforms do
    branch_no := branch_no+1:
    H := tr["expression"]:
    GCDpoly := gcd(F,H):
    if degree(GCDpoly,{T,X,Y,W})>0 then
        nontrivial_gcd := nontrivial_gcd+1:
        A := `a4b_solver/PrimitivePolynomial`(quo(F,GCDpoly,W),vars):
        B := `a4b_solver/PrimitivePolynomial`(quo(H,GCDpoly,W),vars):
    else
        A := F: B := H:
    end if:
    R := resultant(A,B,W):
    if R=0 then
        zero_resultants := zero_resultants+1:
        printf("%02d,%s,gcddeg=%d,resultant=ZERO\n",
               branch_no,tr["id"],degree(GCDpoly,{T,X,Y,W})):
    else
        nonzero_resultants := nonzero_resultants+1:
        R := `a4b_solver/PrimitivePolynomial`(R,[T,X,Y]):
        Rfac := factors(R):
        printf("%02d,%s,gcddeg=%d,deg=%a,terms=%d,factors=%d\n",
               branch_no,tr["id"],degree(GCDpoly,{T,X,Y,W}),
               [degree(R,T),degree(R,X),degree(R,Y)],
               nops([op(expand(R))]),nops(Rfac[2])):
    end if:
end do:

printf("SUMMARY branches=%d nontrivial_gcd=%d zero_resultants=%d nonzero_resultants=%d\n",
       nops(transforms),nontrivial_gcd,zero_resultants,nonzero_resultants):
quit:
