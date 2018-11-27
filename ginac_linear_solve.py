from mathematica_printer import mstr
import tempfile
import os
import sympy
from sympy import symbols, I, pretty, sympify, Matrix, S
import shutil
from subprocess import check_output
from sympy.parsing.sympy_parser import (parse_expr, convert_xor,
                                        auto_number, auto_symbol_add)


transformations = (convert_xor, auto_number, auto_symbol_add)

command = './linsolve_gauss'

def linear_solve(sparse_mat_rep, sub_cvector, length, fvars, iofvars,
                 fvargen, newfvars, vardict):
    matrix = sympy.zeros(len(sub_cvector), length)
    for index, value in sparse_mat_rep.items():
        matrix[index[0],index[1]] = value
    if matrix.cols == 1 and matrix.rows == 1:
        return {fvars[0]: sub_cvector[0]/matrix[0,0]}
    M_str = mstr(matrix)
    b_str = '{'+','.join('{'+mstr(el)+'}'for el in sub_cvector) + '}'
    script = M_str + '\n' + b_str
    script_file = tempfile.NamedTemporaryFile(mode = 'wt',
                                              delete=False,
                                              dir = './',
                                              suffix=".m")
    checkstring = "Not set."
    script_file.write(script)
    script_file.close()
    del script
    try:
        checkstring = check_output([command, script_file.name])[2:-3].decode("utf-8")
        auto_symbol_add.newfvars = newfvars
        auto_symbol_add.fvargen = fvargen
        temp_vardict = vardict.copy()
        sols_list = [parse_expr(sol,
                                transformations = transformations,
                                local_dict = temp_vardict)
                     for sol in checkstring.split('],[')]
    except Exception as e:
        shutil.copy(script_file.name, './failed_ginac_solve.txt')
        os.remove(script_file.name)
        print(str(e))
        print(checkstring)
        return {}
    else:
        os.remove(script_file.name)
    if not sols_list:
        return {}
    sols = dict(zip(fvars, sols_list))
    return {var: sol for var, sol in sols.items() if var is not sol}
