restart:
read "a4b_rational_core.mpl":
read "trig_polynomial_engine.mpl":
read "appendix_all_3vertex_cases.txt":
read "appendix_reference_polys_by_heading.mpl":
read "unique_2vertex_case.txt":

Check := proc(ok,msg)
    if not evalb(ok) then error "%1",msg end if
end proc:
try
byID := table(): seen := table(): models := table():
primary := 0: alternatives := 0:
Check(nops(appendix_all_3vertex_cases)=257,"expected 257 three-vertex records"):
for row in appendix_all_3vertex_cases do
    id := row[1]:
    Check(not assigned(byID[id]),cat("duplicate identifier: ",id)):
    key := {op(row[2])}:
    Check(not assigned(seen[key]),cat("duplicate vertex combination: ",id)):
    byID[id] := row: seen[key] := id:
    if StringTools:-Search("-OR",id)=0 then primary:=primary+1
    else alternatives:=alternatives+1 end if:
    T := `a4b_solver/MakeModel`(row[2]):
    Check(T["status"]="ok",cat("invalid affine model: ",id)):
    av := T["affine_expressions"]:
    for v in row[2] do
        Check(normal(add(v[j]*av[j],j=1..5)-2)=0,cat("vertex equation: ",id))
    end do:
    Check(normal(add(av[j],j=1..5)-3-4*av[6])=0,cat("angle sum: ",id)):
    models[id] := av:
end do:
Check(primary=230 and alternatives=27,"expected 230 primary and 27 OR records"):
Check(not assigned(byID["APP-158"]),"APP-158 must be retired"):
refs := table():
for row in appendix_reference_polys_by_heading do
    Check(not assigned(refs[row[1]]),"duplicate reference ID"):
    Check(assigned(byID[row[1]]),cat("orphan reference: ",row[1])):
    Check(StringTools:-Search("-OR",row[1])=0,"reference must have primary ID"):
    refs[row[1]] := row[2]:
end do:
Check(nops(appendix_reference_polys_by_heading)=primary,"reference coverage"):
Check(unique_2v_case[2]=[[1,0,0,1,1],[0,2,1,0,0]],"two-vertex identity"):
Check(models["APP-177"]=models["APP-177-OR1"],"OR affine equivalence"):
base := byID["APP-177"]: alt := byID["APP-177-OR1"]:
Check(base[3..4]=alt[3..4],"OR phase basis"):
P := `a4b_solver/PhasePolynomial`(`a4b_solver/MakeModel`(base[2]),base[3],base[4]):
Q := `a4b_solver/PhasePolynomial`(`a4b_solver/MakeModel`(alt[2]),alt[3],alt[4]):
Check(normal(P-Q)=0,"OR phase polynomial equality"):
candidate := [6/7,4/7,2/7,3/7,1,1/28]:
for v in byID["APP-133"][2] do
    Check(add(v[j]*candidate[j],j=1..5)=2,"transferred candidate vertex")
end do:
Check(add(candidate[j],j=1..5)=3+4*candidate[6],"transferred candidate angle sum"):
printf("PASS: 257 distinct three-vertex records = 230 primary + 27 OR; one separate two-vertex case.\n"):
printf("PASS: all affine equations, primary reference coverage, transferred candidate, and APP-177 OR equivalence.\n"):
catch:
    printf("FAIL: %a\n",lasterror):
    quit:
end try:
quit:
