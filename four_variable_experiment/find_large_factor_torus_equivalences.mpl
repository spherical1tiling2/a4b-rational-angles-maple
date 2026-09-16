restart:
read "recursive_torsion_engine.mpl":
read "experimental_recursive_utils.mpl":
read "four_to_three_unique.m":

T:='T': X:='X': Y:='Y': vars:=[T,X,Y]:
largeIndices:=[seq(j,j=27..42)]:
canonicalLarge:=table():
for j in largeIndices do
    canonicalLarge[j]:=`a4b_solver/CanonicalScalarPolynomial`(uniqueFactors[j],vars)
end do:

equivalences:=[]:
for source in largeIndices do
  P0:=uniqueFactors[source]:
  for doSwap in [false,true] do
    if doSwap then P1:=subs({X=Y,Y=X},P0) else P1:=P0 end if:
    for signMask from 0 to 7 do
      sT:=if irem(signMask,2)=1 then -1 else 1 end if:
      sX:=if irem(iquo(signMask,2),2)=1 then -1 else 1 end if:
      sY:=if irem(iquo(signMask,4),2)=1 then -1 else 1 end if:
      P2:=subs({T=sT*T,X=sX*X,Y=sY*Y},P1):
      for inverseMask from 0 to 7 do
        iT:=irem(inverseMask,2):
        iX:=irem(iquo(inverseMask,2),2):
        iY:=irem(iquo(inverseMask,4),2):
        P3:=subs({T=`if`(iT=1,1/T,T),
                  X=`if`(iX=1,1/X,X),
                  Y=`if`(iY=1,1/Y,Y)},P2):
        C:=`a4b_solver/CanonicalScalarPolynomial`(P3,vars):
        for target from source+1 to 42 do
          if evalb(expand(C-canonicalLarge[target])=0) then
            rec:=[source,target,doSwap,[sT,sX,sY],[iT,iX,iY]]:
            if not member(rec,equivalences) then
              equivalences:=[op(equivalences),rec]:
              printf("source=%d target=%d swap=%a signs=%a inversions=%a\n",
                     source,target,doSwap,[sT,sX,sY],[iT,iX,iY])
            end if:
          end if:
        end do:
      end do:
    end do:
  end do:
end do:
printf("equivalence_count=%d\n",nops(equivalences)):
save equivalences,
     "full_run/large_factor_equivalences.m":
quit:
