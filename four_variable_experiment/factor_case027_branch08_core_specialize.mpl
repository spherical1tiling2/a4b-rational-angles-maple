restart:
read "full_run/case_027/branch_08_sqrfree.m":
T:='T': X:='X': P:=Sq[2][3][1]:
cX:=content(P,X):
printf("content_in_X=%a generic_degX=%d\n",cX,degree(P,X)):
for t0 in [2,3,-1,5] do
    p0:=expand(subs(T=t0,P)):
    if degree(p0,X)<>degree(P,X) then
        printf("T=%d skipped degree=%d\n",t0,degree(p0,X)):
        next
    end if:
    startCPU:=time(): f0:=factors(p0)[2]: elapsed:=time()-startCPU:
    printf("T=%d cpu=%.2f factors=%d degrees=%a\n",t0,elapsed,nops(f0),
           [seq([degree(f0[j][1],X),f0[j][2]],j=1..nops(f0))]):
    if nops(f0)=1 and degree(f0[1][1],X)=degree(P,X) and f0[1][2]=1 then
        printf("IRREDUCIBILITY_CERTIFICATE specialization_T=%d\n",t0):
        certT:=t0: certFactorization:=f0:
        save certT,certFactorization,cX,
             "full_run/case_027/branch_08_core_irreducible_certificate.m":
        break
    end if:
end do:
quit:
