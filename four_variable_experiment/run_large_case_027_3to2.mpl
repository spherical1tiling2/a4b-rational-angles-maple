restart:
with(PolynomialTools):

read "recursive_torsion_engine.mpl":
read "experimental_recursive_utils.mpl":
read "support_lattice_nd.mpl":
read "four_to_three_unique.m":

caseIndex:=27:
baseDir:="full_run/case_027/":
T:='T': X:='X': Y:='Y': originalVars:=[T,X,Y]:
sourceP:=uniqueFactors[caseIndex]:
L3:=`a4b_solver/SupportLatticeReductionND`(sourceP,originalVars):
if L3["rank"]<>3 then error "case 027 unexpectedly dropped rank" end if:
P:=L3["polynomial"]: vars3:=L3["variables"]:
elim:=vars3[3]: vars2:=[vars3[1],vars3[2]]:
transforms:=`a4b_solver/SignSquareTransforms`(P,vars3):

printf("CASE=%d terms=%d deg=%a lattice=%s index=%a branches=%d\n",
       caseIndex,nops([op(expand(P))]),[seq(degree(P,v),v in vars3)],
       L3["status"],L3["index"],nops(transforms)):

for branchIndex from 1 to nops(transforms) do
    branchFile:=cat(baseDir,sprintf("branch_%02d.m",branchIndex)):
    if FileTools:-Exists(branchFile) then
        printf("BRANCH=%02d SKIP checkpoint\n",branchIndex):
        next
    end if:
    tr:=transforms[branchIndex]:
    resultantFile:=cat(baseDir,sprintf("branch_%02d_resultant.m",branchIndex)):
    if FileTools:-Exists(resultantFile) then
        read resultantFile:
        printf("BRANCH=%02d RESULTANT checkpoint terms=%d deg=%a\n",
               branchIndex,if R=0 then 0 else nops([op(expand(R))]) end if,
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
        printf("BRANCH=%02d RESULTANT done cpu=%.2f terms=%d deg=%a\n",
               branchIndex,resultantCPU,
               if R=0 then 0 else nops([op(expand(R))]) end if,
               if R=0 then [] else [seq(degree(R,v),v in vars2)] end if):
    end if:
    branchFactors:=[]: branchMaps:=[]: branchMultiplicity:=[]:
    rawFactorCount:=0:
    if R<>0 then
        factoredFile:=cat(baseDir,sprintf("branch_%02d_factored.m",branchIndex)):
        if FileTools:-Exists(factoredFile) then
            read factoredFile:
            printf("BRANCH=%02d FACTOR checkpoint count=%d\n",branchIndex,nops(rawFactors)):
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
            printf("BRANCH=%02d FACTOR done cpu=%.2f count=%d\n",
                   branchIndex,factorCPU,nops(rawFactors)):
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
                if maprec=branchMaps[j] and evalb(expand(C-branchFactors[j])=0) then
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
                 resultantCPU,factorCPU,rawFactorCount,
                 nops(branchFactors),
                 if R=0 then 0 else nops([op(expand(R))]) end if,
                 if R=0 then [] else [seq(degree(R,v),v in vars2)] end if]:
    sourceLatticeMap:=[L3["rank"],L3["status"],L3["index"],L3["basis"],
                      L3["original_variables"],L3["variables"],
                      L3["base_exponent"],L3["coordinate_shift"]]:
    save branchFactors,branchMaps,branchMultiplicity,branchMeta,sourceLatticeMap,
         branchFile:
    printf("BRANCH=%02d id=%a gcd=%d resCPU=%.2f factorCPU=%.2f raw=%d stored=%d resTerms=%d resDeg=%a\n",
           branchIndex,tr["id"],gcdDegree,resultantCPU,factorCPU,
           rawFactorCount,nops(branchFactors),branchMeta[9],branchMeta[10]):
end do:
printf("CASE=%d COMPLETE_3TO2\n",caseIndex):
quit:
