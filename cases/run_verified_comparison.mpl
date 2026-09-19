restart:
kernelopts(numcpus=1):
read "exact_angle_checks.mpl":
read "appendix_all_3vertex_cases.txt":
read "search_candidates.mpl":
read "family_lift_candidates.mpl":
Check:=proc(ok,msg) if not evalb(ok) then error "%1",msg end if end proc:
try
byID:=table(): for row in appendix_all_3vertex_cases do byID[row[1]]:=row end do:
expected:={}:
for line in StringTools:-Split(FileTools:-Text:-ReadFile("appendix_3vertex_solution_candidates.txt"),"\n") do
    if StringTools:-Search("APP-",line)=1 then
        fields:=StringTools:-Split(line,"|"):
        nums:=map(parse,StringTools:-Split(fields[3],",")):
        expected:=expected union {[fields[1],parse(fields[2]),nums/parse(fields[4])]}:
    end if:
end do:
observed:={op(search_candidates),op(degree_excluded_candidates),op(family_lift_candidates)}:
accepted:={}: rejected:={}: exactCount:=0:
out:=fopen("finite_candidate_checks.csv",WRITE):
fprintf(out,"case,f,angles_pi,vertex_exact,angle_sum_exact,equation_26_exact,minimum_f_from_vertex_degrees,status,source\n"):
for item in sort([op(observed)]) do
    id,fv,ang:=op(item):
    Check(assigned(byID[id]),"unknown candidate case"):
    Check(type(fv,even) and fv>=12,"tile count"):
    Check(andmap(a->evalb(0<a and a<2),ang),"angle bounds"):
    Check(ang[2]<>ang[3] and ang[4]<>ang[5],"nonsymmetric domain"):
    for v in byID[id][2] do Check(add(v[j]*ang[j],j=1..5)=2,"candidate vertex equation") end do:
    Check(add(ang)=3+4/fv,"candidate angle sum"):
    Check(`a4b_solver/ExactEquation26`(ang),cat("exact equation (2.6): ",id)):
    exactCount:=exactCount+1:
    minF:=12+2*add(add(v)-3,v in byID[id][2]):
    if fv<minF then status:="excluded_vertex_degree_bound": rejected:=rejected union {item}
    else status:="matches_manuscript": accepted:=accepted union {item} end if:
    if member(item,family_lift_candidates) then origin:="phase_family_lift"
    else origin:="recursive_search" end if:
    fprintf(out,"%s,%a,\"%a\",true,true,true,%a,%s,%s\n",id,fv,ang,minF,status,origin):
end do:
fclose(out):
missing:=expected minus accepted: extra:=accepted minus expected:
printf("COMPARISON expected=%d accepted=%d excluded=%d missing=%d extra=%d exact_checked=%d\n",nops(expected),nops(accepted),nops(rejected),nops(missing),nops(extra),exactCount):
if nops(missing)>0 then printf("MISSING: %a\n",missing) end if:
if nops(extra)>0 then printf("EXTRA: %a\n",extra) end if:
Check(nops(missing)=0 and nops(extra)=0,"manuscript candidate comparison"):
out:=fopen("manuscript_finite_candidates.csv",WRITE):
fprintf(out,"case,f,angles_pi\n"):
for item in sort([op(accepted)]) do fprintf(out,"%s,%a,\"%a\"\n",op(item)) end do:
fclose(out):
printf("PASS: all finite manuscript records agree after exact checks and the vertex degree bound.\n"):
catch:
    printf("FAIL: %a\n",lastexception):
end try:
quit:
