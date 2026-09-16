restart:
read "full_run/case_027/branch_08_resultant.m":
T:='T': X:='X': vars:=[T,X]:
printf("input terms=%d deg=%a\n",nops([op(expand(R))]),[degree(R,T),degree(R,X)]):
startCPU:=time():
Sq:=sqrfree(R,vars):
printf("sqrfree cpu=%.2f content=%a components=%d\n",time()-startCPU,Sq[1],nops(Sq[2])):
for j from 1 to nops(Sq[2]) do
    P:=Sq[2][j][1]:
    printf("component=%d multiplicity=%d terms=%d deg=%a\n",j,Sq[2][j][2],
           nops([op(expand(P))]),[degree(P,T),degree(P,X)]):
end do:
save Sq,"full_run/case_027/branch_08_sqrfree.m":
quit:
