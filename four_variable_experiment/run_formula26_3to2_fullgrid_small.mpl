restart:
with(PolynomialTools):

read "recursive_torsion_engine.mpl":
read "experimental_recursive_utils.mpl":
read "support_lattice_nd.mpl":
read "four_to_three_unique.m":

T:='T': X:='X': Y:='Y': originalVars:=[T,X,Y]:
checkpointFile := "three_to_two_fullgrid_small.m":
if FileTools:-Exists(checkpointFile) then
    read checkpointFile:
    printf("RESUME processed=%a rankDeficient=%a deferred=%a\n",
           processedCases,rankDeficientCases,deferredCases):
else
    caseResults := table(): familyProjections := table():
    processedCases := []: deferredCases := []: rankDeficientCases := []:
end if:

for caseIndex from 1 to nops(uniqueFactors) do
    if member(caseIndex,processedCases) or member(caseIndex,rankDeficientCases) or
       member(caseIndex,deferredCases) then next end if:
    originalP := uniqueFactors[caseIndex]:
    L3 := `a4b_solver/SupportLatticeReductionND`(originalP,originalVars):
    if L3["rank"]<3 then
        rankDeficientCases := [op(rankDeficientCases),caseIndex]:
        familyProjections[caseIndex] := L3:
        printf("U4_%03d family rank=%d free=%d reduced=%a basis=%a\n",
               caseIndex,L3["rank"],3-L3["rank"],L3["polynomial"],L3["basis"]):
        next
    end if:

    P := L3["polynomial"]: vars3 := L3["variables"]:
    termCount := nops([op(expand(P))]):
    if termCount>200 then
        deferredCases := [op(deferredCases),caseIndex]:
        printf("U4_%03d deferred-large fullgrid_terms=%d deg=%a\n",
               caseIndex,termCount,[seq(degree(P,v),v in vars3)]):
        next
    end if:

    processedCases := [op(processedCases),caseIndex]:
    vars2 := [vars3[1],vars3[2]]:
    elim := vars3[3]:
    transforms := `a4b_solver/SignSquareTransforms`(P,vars3):
    localFactors := []: localMultiplicity := []:
    rawCount := 0: gcdCount := 0: rankDropCount := 0:
    printf("U4_%03d start lattice=%s index=%a terms=%d deg=%a\n",
           caseIndex,L3["status"],L3["index"],termCount,
           [seq(degree(P,v),v in vars3)]):

    for tr in transforms do
        Q := tr["expression"]:
        GCDpoly := gcd(P,Q):
        if degree(GCDpoly,{op(vars3)})>0 then
            gcdCount := gcdCount+1:
            # A full-grid irreducible factor should not be invariant.  Keep
            # the branch visible and use only its residual intersection.
            A := `a4b_solver/PrimitivePolynomial`(quo(P,GCDpoly,elim),vars3):
            B := `a4b_solver/PrimitivePolynomial`(quo(Q,GCDpoly,elim),vars3):
        else
            A := P: B := Q:
        end if:
        R := resultant(A,B,elim):
        if R=0 then next end if:
        R := `a4b_solver/PrimitivePolynomial`(R,vars2):
        for entry in factors(R)[2] do
            if degree(entry[1],{op(vars2)})=0 then next end if:
            rawCount := rawCount+1:
            L2 := `a4b_solver/SupportLatticeReductionND`(entry[1],vars2):
            if L2["rank"]<2 then rankDropCount:=rankDropCount+1 end if:
            C := `a4b_solver/CanonicalScalarPolynomial`(L2["polynomial"],L2["variables"]):
            found := `a4b_solver/FindProjectivePolynomial`(C,localFactors,L2["variables"]):
            if found=0 then
                localFactors := [op(localFactors),C]:
                localMultiplicity := [op(localMultiplicity),1]:
            else
                localMultiplicity[found] := localMultiplicity[found]+1:
            end if:
        end do:
    end do:
    caseResults[caseIndex] := table(["source_lattice"=L3,
                                     "factors"=localFactors,
                                     "multiplicity"=localMultiplicity,
                                     "raw_factor_count"=rawCount,
                                     "gcd_branch_count"=gcdCount,
                                     "rank_drop_occurrences"=rankDropCount]):
    printf("U4_%03d done raw2=%d uniqueFullgrid2=%d gcd=%d rankDrops=%d\n",
           caseIndex,rawCount,nops(localFactors),gcdCount,rankDropCount):
    save caseResults,familyProjections,processedCases,deferredCases,
         rankDeficientCases,
         checkpointFile:
end do:

printf("SUMMARY processed=%a rankDeficient=%a deferred=%a\n",
       processedCases,rankDeficientCases,deferredCases):
quit:
