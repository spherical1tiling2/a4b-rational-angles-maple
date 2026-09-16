restart:
with(PolynomialTools):

read "recursive_torsion_engine.mpl":
read "formula26_four_variable_core.mpl":

CanonicalScalar := proc(P,vars::list)
    local G,terms,term,best,bestexp,exps,j,better,c,mon;
    G := `a4b_solver/PrimitivePolynomial`(P,vars):
    terms := [op(expand(G))]:
    best := terms[1]:
    bestexp := [seq(degree(best,vars[j]),j=1..nops(vars))]:
    for term in terms[2..-1] do
        exps := [seq(degree(term,vars[j]),j=1..nops(vars))]:
        better := false:
        for j from 1 to nops(vars) do
            if exps[j]>bestexp[j] then better:=true: break
            elif exps[j]<bestexp[j] then break
            end if:
        end do:
        if better then best:=term: bestexp:=exps end if:
    end do:
    mon := mul(vars[j]^bestexp[j],j=1..nops(vars)):
    c := normal(best/mon):
    return expand(G/c):
end proc:

AddUnique := proc(P,origin,vars::list,$)
    global uniqueFactors,uniqueOrigins,uniqueMultiplicity;
    local C,j;
    C := CanonicalScalar(P,vars):
    for j from 1 to nops(uniqueFactors) do
        if evalb(expand(C-uniqueFactors[j])=0) then
            uniqueOrigins[j] := [op(uniqueOrigins[j]),origin]:
            uniqueMultiplicity[j] := uniqueMultiplicity[j]+1:
            return j
        end if:
    end do:
    uniqueFactors := [op(uniqueFactors),C]:
    uniqueOrigins := [op(uniqueOrigins),[origin]]:
    uniqueMultiplicity := [op(uniqueMultiplicity),1]:
    return nops(uniqueFactors):
end proc:

X := 'X': Y := 'Y': W := 'W': T := 'T':
F := `a4b_solver/Formula26FourVariablePolynomial`(X,Y,W,T):
vars4 := [T,X,Y,W]: vars3 := [T,X,Y]:
transforms := `a4b_solver/SignSquareTransforms`(F,vars4):
uniqueFactors := []: uniqueOrigins := []: uniqueMultiplicity := []:
rawFactorCount := 0:

for tr in transforms do
    R := resultant(F,tr["expression"],W):
    if R<>0 then
        R := `a4b_solver/PrimitivePolynomial`(R,vars3):
        Rfac := factors(R)[2]:
        for entry in Rfac do
            rawFactorCount := rawFactorCount+1:
            AddUnique(entry[1],tr["id"],vars3):
        end do:
    end if:
end do:

printf("=== UNIQUE FACTORS AFTER 4 -> 3 ===\n"):
printf("raw factor occurrences=%d, unique factors=%d\n",
       rawFactorCount,nops(uniqueFactors)):
for j from 1 to nops(uniqueFactors) do
    P := uniqueFactors[j]:
    active := []:
    for v in vars3 do if has(P,v) then active:=[op(active),v] end if end do:
    printf("U4_%03d active=%a deg=%a terms=%d occurrences=%d origins=%a\n",
           j,active,[seq(degree(P,v),v in vars3)],
           nops([op(expand(P))]),uniqueMultiplicity[j],uniqueOrigins[j]):
end do:

save uniqueFactors,uniqueOrigins,uniqueMultiplicity,
     "four_to_three_unique.m":
quit:
