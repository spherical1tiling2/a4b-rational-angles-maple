restart:
read "full_run/case_027/branch_08_sqrfree.m":
T:='T': X:='X': P:=Sq[2][3][1]:
printf("core terms=%d deg=%a\n",nops([op(expand(P))]),[degree(P,T),degree(P,X)]):
startCPU:=time():
FM:=Factors(P) mod 3:
printf("mod3 factor cpu=%.2f count=%d unit=%a\n",time()-startCPU,nops(FM[2]),FM[1]):
for j from 1 to nops(FM[2]) do
    Q:=FM[2][j][1]:
    printf("factor=%d multiplicity=%d terms=%d deg=%a\n",j,FM[2][j][2],
           nops([op(expand(Q))]),[degree(Q,T),degree(Q,X)]):
end do:
save FM,"full_run/case_027/branch_08_core_mod3.m":
quit:
