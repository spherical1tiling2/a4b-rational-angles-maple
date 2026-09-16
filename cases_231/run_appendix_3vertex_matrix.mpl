restart:
read "a4b_rational_core.mpl":
read "appendix_3vertex_cases.txt":

printf("=== APPENDIX SIX THREE-VERTEX CASES ===\n"):
printf("Equation used for the trigonometric stage: paper (2.6) only.\n"):
printf("The extra cosine factor in (2.3) is not used.\n\n"):

for row in appendix_3v_cases do
    id := row[1]: V := row[2]: phase := row[3]: expected := row[4]:
    T := `a4b_solver/MakeModel`(V):
    printf("CASE %s\n",id):
    printf("vertices = %a\n",V):
    printf("phase basis = %s\n",phase):
    printf("status = %a\n",T["status"]):
    if T["status"] = "ok" then
        printf("kind = %a, vertex_rank = %a, extended_rank = %a\n",
               T["kind"],T["vertex_rank"],T["extended_rank"]):
        printf("affine_expressions = %a\n",T["affine_expressions"]):
        printf("appendix_expected_final = %s\n",expected):
    else
        printf("message = %a\n",T["message"]):
    end if:
    printf("\n"):
end do:

quit:
