import sympy

def linear_solve(sparse_mat_rep, sub_cvector, length, fvars, iofvars, fvargen, newfvars, tempgen, tempvars, len_oldfvars):
    augmatrix = sympy.zeros(len(sub_cvector), length+1)
    for index, value in sparse_mat_rep.items():
        augmatrix[index[0],index[1]] = value
    for i, c in enumerate(sub_cvector):
        augmatrix[i,-1] = c
    if augmatrix.cols == 2 and augmatrix.rows == 1:
        return {fvars[0]: augmatrix[0,1]/augmatrix[0,0]}
    return sympy.solve_linear_system(augmatrix,*fvars)
