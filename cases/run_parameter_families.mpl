restart:
kernelopts(numcpus=1):
read "exact_angle_checks.mpl":
read "appendix_all_3vertex_cases.txt":
read "appendix_reference_polys_by_heading.mpl":
Check:=proc(ok,msg) if not evalb(ok) then error "%1",msg end if end proc:
try
byID:=table(): for row in appendix_all_3vertex_cases do byID[row[1]]:=row end do:
refs:=table(): for row in appendix_reference_polys_by_heading do refs[row[1]]:=row[2] end do:
families:=[
 ["APP-011",[4*q,4*q,1-2*q,1-2*q,1],12],
 ["APP-011-OR1",[4*q,4*q,1-2*q,1-2*q,1],12],
 ["APP-202",[8*q,1-12*q,1-4*q,1/2+2*q,1/2+10*q],14]
]:
out:=fopen("parameter_families.csv",WRITE):
fprintf(out,"case,angles_pi_q_equals_1_over_f,even_f_min,vertex_exact,angle_sum_exact,equation_26_exact\n"):
for item in families do
    id:=item[1]: ang:=item[2]:
    for v in byID[id][2] do Check(normal(add(v[j]*ang[j],j=1..5)-2)=0,"family vertex equation") end do:
    Check(normal(add(ang)-3-4*q)=0,"family angle sum"):
    Check(`a4b_solver/FamilyEquation26`(ang,q),"family equation (2.6)"):
    fprintf(out,"%s,\"%a\",%d,true,true,true\n",id,ang,item[3]):
end do:
fclose(out):
# APP-202 has x=exp(2*I*Pi*delta), y=exp(4*I*Pi/f).
# The factor x+y gives delta=1/2+2/f+k for an integer k.
Check(rem(refs["APP-202"],x+y,x)=0,"APP-202 phase factor"):
vars:=[aa,bb,cc,dd,ee]:
eqs:={seq(add(v[j]*vars[j],j=1..5)=2,v in byID["APP-202"][2]),
       add(vars)=3+4*q,dd=1/2+2*q+k}:
sol:=solve(eqs,{op(vars)}):
lift:=map(normal,subs(sol,vars)):
Check(lift=[8*q,1-12*q,1-4*q,1/2+2*q+k,1/2+10*q-k],"phase lift derivation"):
# With even f>=12 and positive angles, beta>0 gives f>=14.
# Then 0<delta<2 and epsilon>0 allow k=0, or k=1 with f<20.
# All k<=-1 have delta<0; all k>=2 have delta>2.
Check(`a4b_solver/FamilyEquation26`(subs(k=1,lift),q),"second lift equation (2.6)"):
derived:=[]:
for fv from 14 to 18 by 2 do
    ang:=subs(q=1/fv,k=1,lift):
    Check(andmap(a->evalb(0<a and a<2),ang),"finite lift angle bounds"):
    Check(`a4b_solver/ExactEquation26`(ang),"finite lift exact equation"):
    derived:=[op(derived),["APP-202",fv,ang]]:
end do:
out:=fopen("family_lift_candidates.mpl",WRITE):
fprintf(out,"# Finite angle lifts derived from the APP-202 factor x+y.\nfamily_lift_candidates := %a:\n",derived): fclose(out):
printf("PASS: parameter_family_records=3; exact_symbolic_identities=3; independently_derived_finite_lifts=3\n"):
catch:
    printf("FAIL: %a\n",lastexception):
end try:
quit:
