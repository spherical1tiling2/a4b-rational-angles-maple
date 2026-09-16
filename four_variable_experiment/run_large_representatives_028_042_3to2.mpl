restart:
with(PolynomialTools):

read "recursive_torsion_engine.mpl":
read "experimental_recursive_utils.mpl":
read "support_lattice_nd.mpl":
read "four_to_three_unique.m":

# 31,32,39,40 are obtained from 29,30,37,38 by X<->Y.  The remaining
# representatives can be split into three disjoint process groups.
processGroup:=getenv("A4B_REP_GROUP"):
if processGroup="A" then
    representatives:=[28]
elif processGroup="B" then
    representatives:=[29]
elif processGroup="C" then
    representatives:=[30]
elif processGroup="D" then
    representatives:=[33,36]
elif processGroup="E" then
    representatives:=[34,37]
elif processGroup="F" then
    representatives:=[35,38]
elif processGroup="G" then
    representatives:=[41]
elif processGroup="H" then
    representatives:=[42]
else
    representatives:=[28,29,30,33,34,35,36,37,38,41,42]
end if:
printf("PROCESS_GROUP=%a representatives=%a\n",processGroup,representatives):
baseRoot:="full_run/":
T:='T': X:='X': Y:='Y': originalVars:=[T,X,Y]:

for caseIndex in representatives do
    baseDir:=cat(baseRoot,sprintf("case_%03d/",caseIndex)):
    if not FileTools:-Exists(baseDir) then FileTools:-MakeDirectory(baseDir) end if:
    completeFile:=cat(baseDir,"case_complete.m"):
    if FileTools:-Exists(completeFile) then
        printf("CASE=%d SKIP complete checkpoint\n",caseIndex):
        next
    end if:

    sourceP:=uniqueFactors[caseIndex]:
    L3:=`a4b_solver/SupportLatticeReductionND`(sourceP,originalVars):
    if L3["rank"]<>3 then error "large representative unexpectedly dropped rank" end if:
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
        tr:=transforms[branchIndex]:
        resultantFile:=cat(baseDir,sprintf("branch_%02d_resultant.m",branchIndex)):
        if FileTools:-Exists(resultantFile) then
            read resultantFile:
            printf("CASE=%d BRANCH=%02d RESULTANT checkpoint terms=%d deg=%a\n",
                   caseIndex,branchIndex,
                   if R=0 then 0 else nops([op(expand(R))]) end if,
                   if R=0 then [] else [seq(degree(R,v),v in vars2)] end if):
        else
            startCPU:=time():
            Q:=tr["expression"]:
            GCDpoly:=gcd(P,Q):
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
            save R,gcdDegree,resultantCPU,resultantFile:
            printf("CASE=%d BRANCH=%02d RESULTANT done cpu=%.2f terms=%d deg=%a\n",
                   caseIndex,branchIndex,resultantCPU,
                   if R=0 then 0 else nops([op(expand(R))]) end if,
                   if R=0 then [] else [seq(degree(R,v),v in vars2)] end if):
        end if:

        branchFactors:=[]: branchMaps:=[]: branchMultiplicity:=[]:
        rawFactorCount:=0:
        if R<>0 then
            factoredFile:=cat(baseDir,sprintf("branch_%02d_factored.m",branchIndex)):
            if FileTools:-Exists(factoredFile) then
                read factoredFile:
                printf("CASE=%d BRANCH=%02d FACTOR checkpoint count=%d\n",
                       caseIndex,branchIndex,nops(rawFactors)):
            else
                factorStart:=time():
                SqLocal:=sqrfree(R,vars2):
                rawFactors:=[]:
                for sqEntry in SqLocal[2] do
                    component:=sqEntry[1]: componentMultiplicity:=sqEntry[2]:
                    componentFactors:=factors(component,method="Wang")[2]:
                    for componentEntry in componentFactors do
                        rawFactors:=[op(rawFactors),
                            [componentEntry[1],
                             componentMultiplicity*componentEntry[2]]]
                    end do:
                end do:
                factorCPU:=time()-factorStart:
                save rawFactors,factorCPU,factoredFile:
                printf("CASE=%d BRANCH=%02d FACTOR done cpu=%.2f count=%d\n",
                       caseIndex,branchIndex,factorCPU,nops(rawFactors)):
            end if:

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
                       evalb(expand(C-branchFactors[j])=0) then
                        found:=j: break
                    end if:
                end do:
                if found=0 then
                    branchFactors:=[op(branchFactors),C]:
                    branchMaps:=[op(branchMaps),maprec]:
                    branchMultiplicity:=[op(branchMultiplicity),entry[2]]:
                else
                    branchMultiplicity[found]:=branchMultiplicity[found]+entry[2]
                end if:
            end do:
        else
            factorCPU:=0:
        end if:

        branchMeta:=[caseIndex,branchIndex,tr["id"],gcdDegree,
                     resultantCPU,factorCPU,rawFactorCount,nops(branchFactors),
                     if R=0 then 0 else nops([op(expand(R))]) end if,
                     if R=0 then [] else [seq(degree(R,v),v in vars2)] end if]:
        sourceLatticeMap:=[L3["rank"],L3["status"],L3["index"],L3["basis"],
                          L3["original_variables"],L3["variables"],
                          L3["base_exponent"],L3["coordinate_shift"]]:
        save branchFactors,branchMaps,branchMultiplicity,branchMeta,sourceLatticeMap,
             branchFile:
        printf("CASE=%d BRANCH=%02d COMPLETE id=%a gcd=%d raw=%d stored=%d\n",
               caseIndex,branchIndex,tr["id"],gcdDegree,
               rawFactorCount,nops(branchFactors)):
    end do:
    caseComplete:=caseIndex:
    save caseComplete,completeFile:
    printf("CASE=%d COMPLETE_3TO2\n",caseIndex):
end do:
printf("ALL_LARGE_REPRESENTATIVES COMPLETE_3TO2\n"):
quit:
