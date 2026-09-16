restart:
read "recursive_torsion_engine.mpl":
read "experimental_recursive_utils.mpl":
read "support_lattice_nd.mpl":

X:='X': Y:='Y': W:='W':
Q:=W^2*X^2*Y^2+W*(Y-X)*(X*Y+1)-1:
vars3:=[X,Y,W]: vars2:=[X,Y]:
L3:=`a4b_solver/SupportLatticeReductionND`(Q,vars3):
printf("=== alpha=4/f: residual quadratic torsion branch ===\n"):
printf("Q=%a\n",Q):
printf("lattice status=%s rank=%d index=%a basis=%a\n",
       L3["status"],L3["rank"],L3["index"],L3["basis"]):
printf("reduced Q=%a in variables %a\n",L3["polynomial"],L3["variables"]):

P:=L3["polynomial"]: active:=L3["variables"]:
transforms:=`a4b_solver/SignSquareTransforms`(P,active):
unique2:=[]: unique2Multiplicity:=[]: unique2Maps:=[]:
gcdCount:=0: rawCount:=0:
for tr in transforms do
    H:=tr["expression"]: GCDpoly:=gcd(P,H):
    if degree(GCDpoly,{op(active)})>0 then gcdCount:=gcdCount+1 end if:
    R:=resultant(P,H,active[3]):
    if R=0 then next end if:
    R:=`a4b_solver/PrimitivePolynomial`(R,[active[1],active[2]]):
    for entry in factors(R)[2] do
        if degree(entry[1],{active[1],active[2]})=0 then next end if:
        rawCount:=rawCount+1:
        L2:=`a4b_solver/SupportLatticeReductionND`(entry[1],[active[1],active[2]]):
        C:=`a4b_solver/CanonicalScalarPolynomial`(L2["polynomial"],L2["variables"]):
        maprec:=[L2["rank"],L2["status"],L2["index"],L2["basis"],
                 L2["original_variables"],L2["variables"],
                 L2["base_exponent"],L2["coordinate_shift"]]:
        found:=0:
        for j from 1 to nops(unique2) do
            if maprec=unique2Maps[j] and
               evalb(expand(C-unique2[j])=0) then found:=j: break end if
        end do:
        if found=0 then
            unique2:=[op(unique2),C]:
            unique2Multiplicity:=[op(unique2Multiplicity),1]:
            unique2Maps:=[op(unique2Maps),maprec]:
        else
            unique2Multiplicity[found]:=unique2Multiplicity[found]+1
        end if:
    end do:
end do:
printf("3->2 raw=%d unique-with-map=%d gcd=%d\n",rawCount,nops(unique2),gcdCount):
for j from 1 to nops(unique2) do
    printf("AQ_%03d rank=%d index=%a basis=%a poly=%a multiplicity=%d\n",
           j,unique2Maps[j][1],unique2Maps[j][3],
           unique2Maps[j][4],unique2[j],unique2Multiplicity[j]):
end do:
save Q,L3,unique2,unique2Multiplicity,unique2Maps,
     "alpha4_quadratic_3to2.m":
quit:
