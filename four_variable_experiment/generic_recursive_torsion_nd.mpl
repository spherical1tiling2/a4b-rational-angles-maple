# Generic dimension-recursive wrapper for the root-of-unity elimination tree.
# It uses the existing SignSquareTransforms/ResultantStep routines and
# recursively reduces n variables to n-1 variables.  Inactive variables are
# reported as free variables instead of being silently eliminated.

if not assigned(a4b_solver) then a4b_solver := table() end if:

`a4b_solver/GenericRecursiveTorsion` := proc(F, vars::list)
    local G,active,free,elim,kept,step,b,children,n,sub,v,j;
    G := `a4b_solver/PrimitivePolynomial`(F,vars):
    active := []: free := []:
    for v in vars do
        if has(G,v) then active := [op(active),v]
        else free := [op(free),v]
        end if:
    end do:
    n := nops(active):

    if n=0 then
        return table(["dimension"=0,"polynomial"=G,
                      "active_variables"=[],"free_variables"=free,
                      "status"=if G=0 then "identity" else "constant" end if]):
    elif n=1 then
        return table(["dimension"=1,"polynomial"=G,
                      "active_variables"=active,"free_variables"=free,
                      "status"="univariate",
                      "diagnostics"=`a4b_solver/UnitRootDiagnostics`(G,active[1])]):
    end if:

    # Eliminate the last active variable.  Each n-to-(n-1) step has
    # (2^n-1) sign branches plus 2^n square branches.
    elim := active[n]:
    kept := [seq(active[j],j=1..n-1)]:
    step := `a4b_solver/ResultantStep`(G,active,elim):
    children := []:
    for b in step["branches"] do
        if b["resultant"]<>0 then
            sub := `a4b_solver/GenericRecursiveTorsion`(b["resultant"],kept):
            children := [op(children),table(["id"=b["id"],
                         "resultant"=b["resultant"],"subtree"=sub])]:
        end if:
    end do:
    return table(["dimension"=n,"polynomial"=G,
                  "active_variables"=active,"free_variables"=free,
                  "eliminated"=elim,"retained"=kept,
                  "branch_count"=step["branch_count"],
                  "nonzero_branch_count"=nops(children),
                  "children"=children,"status"="projected"]):
end proc:

`a4b_solver/PrintGenericSummary` := proc(T,indent::string)
    local C,c;
    printf("%sdimension=%a active=%a free=%a status=%a\n",indent,
        T["dimension"],T["active_variables"],T["free_variables"],T["status"]):
    if assigned(T["branch_count"]) then
        printf("%sbranches=%a nonzero=%a eliminated=%a retained=%a\n",indent,
            T["branch_count"],T["nonzero_branch_count"],
            T["eliminated"],T["retained"]):
        for C in T["children"] do
            printf("%sbranch %a: ",indent,C["id"]):
            c := C["subtree"]:
            printf("dimension=%a active=%a free=%a status=%a\n",
                c["dimension"],c["active_variables"],
                c["free_variables"],c["status"]):
        end do:
    end if:
    return NULL:
end proc:

NULL:
