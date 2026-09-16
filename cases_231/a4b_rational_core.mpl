with(LinearAlgebra):

`a4b_solver/MakeModel` := proc(vertices)
    local V, Vaug, VU, A, Aug, R, vars, rankV, rankVaug,
          rankVU, pivotcols, freecols, pivotrow, xexpr,
          P, T, i, j, k, p, rhs;

    vars := [a,b,c,d,e,q]:

    if nops(vertices) <> 3 then
        error "vertices must contain exactly three vectors":
    end if:
    for i from 1 to 3 do
        if nops(vertices[i]) <> 5 then
            error "each vertex vector must have length five":
        end if:
    end do:

    V := Matrix(vertices):
    rhs := Vector([2,2,2]):
    Vaug := <V | rhs>:
    rankV := Rank(V):
    rankVaug := Rank(Vaug):

    T := table():
    T["vertices"] := vertices:
    T["variables"] := vars:
    T["vertex_rank"] := rankV:
    T["vertex_augmented_rank"] := rankVaug:

    if rankVaug > rankV then
        T["status"] := "inconsistent_vertex_equations":
        return T
    end if:
    if rankV < 3 then
        T["status"] := "unsupported_affine_dimension":
        T["message"] := "vertex rank < 3; the current two-parameter kernel stops":
        return T
    end if:

    VU := <V; Matrix(1,5,1)>:
    rankVU := Rank(VU):
    T["extended_rank"] := rankVU:
    if rankVU = 4 then
        T["kind"] := "free_f":
    elif rankVU = 3 then
        T["kind"] := "fixed_f":
    else
        T["status"] := "unsupported_extended_rank":
        return T
    end if:

    A := Matrix(4,6,0):
    for i from 1 to 3 do
        for j from 1 to 5 do
            A[i,j] := V[i,j]:
        end do:
    end do:
    for j from 1 to 5 do
        A[4,j] := 1:
    end do:
    A[4,6] := -4:

    Aug := Matrix(4,7,0):
    for i from 1 to 4 do
        for j from 1 to 6 do
            Aug[i,j] := A[i,j]:
        end do:
    end do:
    Aug[1,7] := 2:
    Aug[2,7] := 2:
    Aug[3,7] := 2:
    Aug[4,7] := 3:

    R := ReducedRowEchelonForm(Aug):
    pivotcols := []:
    pivotrow := table():

    for i from 1 to 4 do
        p := 0:
        for j from 1 to 6 do
            if R[i,j] <> 0 then
                p := j:
                break:
            end if:
        end do:
        if p <> 0 then
            pivotcols := [op(pivotcols),p]:
            pivotrow[p] := i:
        elif R[i,7] <> 0 then
            T["status"] := "inconsistent_full_system":
            return T
        end if
    end do:

    freecols := []:
    for j from 1 to 6 do
        if not member(j,pivotcols) then
            freecols := [op(freecols),j]:
        end if
    end do:

    T["coefficient_matrix"] := A:
    T["rref"] := R:
    T["pivot_columns"] := pivotcols:
    T["free_columns"] := freecols:

    if nops(freecols) <> 2 then
        T["status"] := "unexpected_affine_dimension":
        return T
    end if:

    P := [_u1,_u2]:
    xexpr := Vector(6,0):
    for k from 1 to 2 do
        xexpr[freecols[k]] := P[k]:
    end do:
    for k from 1 to nops(pivotcols) do
        p := pivotcols[k]:
        i := pivotrow[p]:
        xexpr[p] := normal(R[i,7] -
            add(R[i,freecols[j]]*P[j],j=1..2)):
    end do:

    T["status"] := "ok":
    T["temporary_parameters"] := P:
    T["affine_expressions"] := [seq(normal(xexpr[j]),j=1..6)]:
    T["normalized_domain"] := [
        seq(0 < T["affine_expressions"][j],j=1..5),
        seq(T["affine_expressions"][j] < 2,j=1..5)
    ]:

    if T["kind"] = "free_f" then
        T["normalized_domain"] := [op(T["normalized_domain"]),
            0 < T["affine_expressions"][6],
            T["affine_expressions"][6] <= 1/12]:
    else
        T["fixed_f"] := normal(1/T["affine_expressions"][6]):
    end if:

    return T
end proc:

`a4b_solver/PrintModel` := proc(T)
    printf("status = %a\n",T["status"]):
    if T["status"] <> "ok" then
        return NULL:
    end if:
    printf("kind = %a\n",T["kind"]):
    printf("vertex_rank = %a, extended_rank = %a\n",
        T["vertex_rank"],T["extended_rank"]):
    printf("coordinates = %a\n",T["temporary_parameters"]):
    printf("affine_expressions = %a\n",T["affine_expressions"]):
    if T["kind"] = "fixed_f" then
        printf("fixed_f = %a\n",T["fixed_f"]):
    end if:
    printf("normalized_domain = %a\n",T["normalized_domain"]):
    return NULL:
end proc:
