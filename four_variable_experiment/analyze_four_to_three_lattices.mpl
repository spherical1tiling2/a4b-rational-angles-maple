restart:
read "recursive_torsion_engine.mpl":
read "support_lattice_nd.mpl":
read "four_to_three_unique.m":

T:='T': X:='X': Y:='Y': vars:=[T,X,Y]:
printf("=== SUPPORT LATTICES OF UNIQUE 3-VARIABLE PROJECTIONS ===\n"):
for j from 1 to nops(uniqueFactors) do
    P:=uniqueFactors[j]:
    L:=`a4b_solver/SupportLatticeReductionND`(P,vars):
    printf("U4_%03d status=%s rank=%d index=%a olddeg=%a oldterms=%d",
           j,L["status"],L["rank"],L["index"],
           [seq(degree(P,v),v in vars)],nops([op(expand(P))])):
    if L["status"]="reduced" or L["status"]="rank_deficient" then
        printf(" newdeg=%a newterms=%d basis=%a",
               [seq(degree(L["polynomial"],v),v in L["variables"])],
               nops([op(expand(L["polynomial"]))]),L["basis"])
    end if:
    printf("\n"):
end do:
quit:
