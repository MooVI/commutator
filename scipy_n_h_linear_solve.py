from scipy.linalg import null_space
import numpy

def n_h_linear_solve(sparse_mat_rep, sub_cvector, length, fvars, iofvars,
                 fvargen, newfvars, numeric_dict):
    matrix = numpy.zeros((len(sub_cvector), length))
    for index, value in sparse_mat_rep.items():
        matrix[index[0],index[1]] = value.xreplace(numeric_dict)
    ns = null_space(matrix)
    sols_list = [0]*length
    for nullvec in ns.T:
            newfvar = next(fvargen)
            newfvars.append(newfvar)
            sols_list = [sol + newfvar*nullvecel
                         for sol, nullvecel in zip(sols_list, nullvec)]
   #ipdb.set_trace()
    if not sols_list:
        return {}
    sols = dict(zip(fvars, sols_list))
    return {var: sol for var, sol in sols.items() if var is not sol}
