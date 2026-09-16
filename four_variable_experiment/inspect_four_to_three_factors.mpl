restart:
read "four_to_three_unique.m":
T:='T': X:='X': Y:='Y':
printf("=== SMALL UNIQUE FACTORS AFTER 4 -> 3 ===\n"):
for j from 1 to nops(uniqueFactors) do
    P:=uniqueFactors[j]: nterms:=nops([op(expand(P))]):
    if nterms<=100 then
        printf("U4_%03d terms=%d deg=%a occurrences=%d\n%a\n",
               j,nterms,[degree(P,T),degree(P,X),degree(P,Y)],
               uniqueMultiplicity[j],P):
    end if:
end do:
quit:
