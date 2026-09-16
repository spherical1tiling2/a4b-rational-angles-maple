restart:
read "full_run/case_027/branch_08_sqrfree.m":
T:='T': X:='X': P:=Sq[2][3][1]:
printf("core terms=%d deg=%a\n",nops([op(expand(P))]),[degree(P,T),degree(P,X)]):
startCPU:=time():
FW:=factors(P,method="Wang"):
printf("Wang cpu=%.2f count=%d\n",time()-startCPU,nops(FW[2])):
for j from 1 to nops(FW[2]) do
    Q:=FW[2][j][1]:
    printf("factor=%d multiplicity=%d terms=%d deg=%a\n",j,FW[2][j][2],
           nops([op(expand(Q))]),[degree(Q,T),degree(Q,X)]):
end do:
save FW,"full_run/case_027/branch_08_core_wang.m":
quit:
