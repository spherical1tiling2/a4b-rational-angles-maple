restart:
read "a4b_rational_core.mpl":
read "appendix_all_3vertex_cases.txt":

logname := "appendix_all_3vertex_affine.log":
interface(errorbreak=0):
logfile := fopen(logname,WRITE):
fprintf(logfile,"A4BR appendix batch; equation (2.6) reserved for the trigonometric stage only.\n"):
fprintf(logfile,"primary appendix headings = %d\n",nops(appendix_all_3vertex_cases)):
fprintf(logfile,"case,status,kind,vertex_rank,extended_rank,affine_dimension\n"):

okcount := 0: badcount := 0:
for row in appendix_all_3vertex_cases do
    id := row[1]: V := row[2]:
    T := `a4b_solver/MakeModel`(V):
    if T["status"] = "ok" then
        okcount := okcount+1:
        fprintf(logfile,"%s,%a,%a,%a,%a,%a\n",id,T["status"],T["kind"],
            T["vertex_rank"],T["extended_rank"],nops(T["free_columns"])):
        fprintf(logfile,"%s affine = %a\n",id,T["affine_expressions"]):
    else
        badcount := badcount+1:
        fprintf(logfile,"%s,%a,,,,\n",id,T["status"]):
        if assigned(T["message"]) then
            fprintf(logfile,"%s message = %a\n",id,T["message"]):
        end if:
    end if:
end do:

fprintf(logfile,"SUMMARY ok=%d bad=%d total=%d\n",okcount,badcount,nops(appendix_all_3vertex_cases)):
fclose(logfile):
printf("SUMMARY ok=%d bad=%d total=%d\n",okcount,badcount,nops(appendix_all_3vertex_cases)):
quit:
