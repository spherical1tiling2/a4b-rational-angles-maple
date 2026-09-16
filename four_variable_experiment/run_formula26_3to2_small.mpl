restart:
with(PolynomialTools):

read "recursive_torsion_engine.mpl":
read "experimental_recursive_utils.mpl":
read "four_to_three_unique.m":

T:='T': X:='X': Y:='Y':
vars3 := [T,X,Y]: vars2 := [T,X]:
twoFactors := []: twoOrigins := []: twoMultiplicity := []:
processed3 := []: skippedLarge3 := []:

for caseIndex from 1 to nops(uniqueFactors) do
    P := uniqueFactors[caseIndex]:
    termCount := nops([op(expand(P))]):
    activeCount := 0:
    for v in vars3 do if has(P,v) then activeCount:=activeCount+1 end if end do:

    if activeCount<3 then
        printf("U4_%03d direct active=%d terms=%d\n",caseIndex,activeCount,termCount):
        next
    elif termCount>200 then
        skippedLarge3 := [op(skippedLarge3),caseIndex]:
        printf("U4_%03d deferred-large terms=%d deg=%a\n",caseIndex,termCount,
               [degree(P,T),degree(P,X),degree(P,Y)]):
        next
    end if:

    processed3 := [op(processed3),caseIndex]:
    transforms := `a4b_solver/SignSquareTransforms`(P,vars3):
    printf("U4_%03d start terms=%d deg=%a branches=%d\n",caseIndex,termCount,
           [degree(P,T),degree(P,X),degree(P,Y)],nops(transforms)):

    localRaw := 0: localNew := 0: localGCD := 0:
    for tr in transforms do
        Q := tr["expression"]:
        GCDpoly := gcd(P,Q):
        if degree(GCDpoly,{T,X,Y})>0 then
            localGCD := localGCD+1:
            A := `a4b_solver/PrimitivePolynomial`(quo(P,GCDpoly,Y),vars3):
            B := `a4b_solver/PrimitivePolynomial`(quo(Q,GCDpoly,Y),vars3):
        else
            A := P: B := Q:
        end if:
        R := resultant(A,B,Y):
        if R<>0 then
            R := `a4b_solver/PrimitivePolynomial`(R,vars2):
            for entry in factors(R)[2] do
                localRaw := localRaw+1:
                C := `a4b_solver/CanonicalScalarPolynomial`(entry[1],vars2):
                found := `a4b_solver/FindProjectivePolynomial`(C,twoFactors,vars2):
                if found=0 then
                    twoFactors := [op(twoFactors),C]:
                    twoOrigins := [op(twoOrigins),[sprintf("U4_%03d/%s",caseIndex,tr["id"])]]:
                    twoMultiplicity := [op(twoMultiplicity),1]:
                    localNew := localNew+1:
                else
                    twoOrigins[found] := [op(twoOrigins[found]),sprintf("U4_%03d/%s",caseIndex,tr["id"])]:
                    twoMultiplicity[found] := twoMultiplicity[found]+1:
                end if:
            end do:
        end if:
    end do:
    printf("U4_%03d done raw2=%d newUnique2=%d gcdBranches=%d totalUnique2=%d\n",
           caseIndex,localRaw,localNew,localGCD,nops(twoFactors)):
    save twoFactors,twoOrigins,twoMultiplicity,processed3,skippedLarge3,
         "three_to_two_small_checkpoint.m":
end do:

printf("SUMMARY processed3=%a deferredLarge=%a unique2=%d\n",
       processed3,skippedLarge3,nops(twoFactors)):
for j from 1 to nops(twoFactors) do
    P := twoFactors[j]:
    printf("U3_%04d deg=%a terms=%d occurrences=%d\n",j,
           [degree(P,T),degree(P,X)],nops([op(expand(P))]),twoMultiplicity[j]):
end do:
quit:
