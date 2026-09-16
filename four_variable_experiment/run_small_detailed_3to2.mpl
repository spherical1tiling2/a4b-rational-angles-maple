restart:
with(PolynomialTools):

read "recursive_torsion_engine.mpl":
read "experimental_recursive_utils.mpl":
read "support_lattice_nd.mpl":
read "four_to_three_unique.m":

# Recompute U4_001--U4_026 with the same detailed, recoverable checkpoint
# format used for the large cases.  The earlier aggregate checkpoint kept
# the reduced factors but not their finite-index lattice maps, so it was not
# sufficient for exact back-substitution.
representatives:=[seq(j,j=1..26)]:
printf("SMALL_PROCESS_GROUP=ALL representatives=%a\n",representatives):

baseRoot:="full_run/":
T:='T': X:='X': Y:='Y': originalVars:=[T,X,Y]:

for caseIndex in representatives do
    baseDir:=cat(baseRoot,sprintf("case_%03d/",caseIndex)):
    if not FileTools:-Exists(baseDir) then FileTools:-MakeDirectory(baseDir) end if:
    completeFile:=cat(baseDir,"case_complete.m"):
    familyFile:=cat(baseDir,"case_rank_deficient.m"):
    if FileTools:-Exists(completeFile) or FileTools:-Exists(familyFile) then
        printf("CASE=%d SKIP checkpoint\n",caseIndex):
        next
    end if:

    sourceP:=uniqueFactors[caseIndex]:
    L3:=`a4b_solver/SupportLatticeReductionND`(sourceP,originalVars):
    sourceLatticeMap:=[L3["rank"],L3["status"],L3["index"],L3["basis"],
                      L3["original_variables"],L3["variables"],
                      L3["base_exponent"],L3["coordinate_shift"]]:
    if L3["rank"]<>3 then
        familyPolynomial:=L3["polynomial"]:
        save sourceLatticeMap,familyPolynomial,familyFile:
        printf("CASE=%d RANK_DEFICIENT rank=%d status=%s\n",
               caseIndex,L3["rank"],L3["status"]):
        next
    end if:

    P:=L3["polynomial"]: vars3:=L3["variables"]:
    elim:=vars3[3]: vars2:=[vars3[1],vars3[2]]:
    transforms:=`a4b_solver/SignSquareTransforms`(P,vars3):
    printf("CASE=%d terms=%d deg=%a lattice=%s index=%a branches=%d\n",
           caseIndex,nops([op(expand(P))]),[seq(degree(P,v),v in vars3)],
           L3["status"],L3["index"],nops(transforms)):

    for branchIndex from 1 to nops(transforms) do
        branchFile:=cat(baseDir,sprintf("branch_%02d.m",branchIndex)):
        if FileTools:-Exists(branchFile) then
            printf("CASE=%d BRANCH=%02d SKIP checkpoint\n",caseIndex,branchIndex):
            next
        end if:
        tr:=transforms[branchIndex]: startCPU:=time():
        Q:=tr["expression"]: GCDpoly:=gcd(P,Q):
        gcdDegree:=degree(GCDpoly,{op(vars3)}):
        if gcdDegree>0 then
            A:=`a4b_solver/PrimitivePolynomial`(quo(P,GCDpoly,elim),vars3):
            B:=`a4b_solver/PrimitivePolynomial`(quo(Q,GCDpoly,elim),vars3):
        else
            A:=P: B:=Q:
        end if:
        R:=resultant(A,B,elim):
        resultantCPU:=time()-startCPU:
        if R<>0 then R:=`a4b_solver/PrimitivePolynomial`(R,vars2) end if:

        branchFactors:=[]: branchMaps:=[]: branchMultiplicity:=[]:
        rawFactorCount:=0: factorStart:=time():
        if R<>0 then
            SqLocal:=sqrfree(R,vars2): rawFactors:=[]:
            for sqEntry in SqLocal[2] do
                componentFactors:=factors(sqEntry[1],method="Wang")[2]:
                for componentEntry in componentFactors do
                    rawFactors:=[op(rawFactors),
                        [componentEntry[1],sqEntry[2]*componentEntry[2]]]
                end do:
            end do:
            for entry in rawFactors do
                if degree(entry[1],{op(vars2)})=0 then next end if:
                rawFactorCount:=rawFactorCount+1:
                L2:=`a4b_solver/SupportLatticeReductionND`(entry[1],vars2):
                C:=`a4b_solver/CanonicalScalarPolynomial`(L2["polynomial"],L2["variables"]):
                maprec:=[L2["rank"],L2["status"],L2["index"],L2["basis"],
                         L2["original_variables"],L2["variables"],
                         L2["base_exponent"],L2["coordinate_shift"]]:
                found:=0:
                for j from 1 to nops(branchFactors) do
                    if maprec=branchMaps[j] and
                       evalb(expand(C-branchFactors[j])=0) then found:=j: break end if:
                end do:
                if found=0 then
                    branchFactors:=[op(branchFactors),C]:
                    branchMaps:=[op(branchMaps),maprec]:
                    branchMultiplicity:=[op(branchMultiplicity),entry[2]]:
                else
                    branchMultiplicity[found]:=branchMultiplicity[found]+entry[2]:
                end if:
            end do:
        end if:
        factorCPU:=time()-factorStart:
        branchMeta:=[caseIndex,branchIndex,tr["id"],gcdDegree,
                     resultantCPU,factorCPU,rawFactorCount,nops(branchFactors),
                     if R=0 then 0 else nops([op(expand(R))]) end if,
                     if R=0 then [] else [seq(degree(R,v),v in vars2)] end if]:
        save branchFactors,branchMaps,branchMultiplicity,branchMeta,sourceLatticeMap,
             branchFile:
        printf("CASE=%d BRANCH=%02d COMPLETE id=%a gcd=%d raw=%d stored=%d\n",
               caseIndex,branchIndex,tr["id"],gcdDegree,
               rawFactorCount,nops(branchFactors)):
    end do:
    caseComplete:=caseIndex: save caseComplete,completeFile:
    printf("CASE=%d COMPLETE_3TO2\n",caseIndex):
end do:
printf("ALL_SMALL_DETAILED COMPLETE_3TO2\n"):
quit:
