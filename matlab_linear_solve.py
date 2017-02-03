import matlab.engine
from matlab_printer import matlabstr
import tempfile
import os
import sympy
from sympy import symbols, I, pretty, sympify, Matrix, S
import shutil

eng = matlab.engine.start_matlab("-nojvm")


def mathematica_parser(exprstring):
    return sympify(exprstring.replace('^', '**'), mathematica_parser.vardict)

def matlab_parser(exprstring, matrix = False):
    if matrix:
        return sympify('Matri'+exprstring[5:].replace('^','**').replace('i', '*I'), mathematica_parser.vardict)
    return sympify(exprstring.replace('^','**').replace('i', '*I'), mathematica_parser.vardict)

def linear_solve(sparse_mat_rep, sub_cvector, length, fvars, iofvars, fvargen, newfvars, tempgen, tempvars, len_oldfvars):
    global eng
    MAX_TRYS = 2
    matrix = sympy.zeros(len(sub_cvector), length)
    for index, value in sparse_mat_rep.items():
        matrix[index[0],index[1]] = value
    #return sympy.solve_linear_system(augmatrix,*fvars)
    if matrix.cols == 1 and matrix.rows == 1:
        return {fvars[0]: sub_cvector[0]/matrix[0,0]}
    M_str = matlabstr(matrix)
    b_str = matlabstr(sympy.Matrix(sub_cvector))
    var_str = ';'.join(name+"=sym('"+name+ "', 'real')"
                       for name in mathematica_parser.vardict)+';\n'
    script = ('function ret = linearsolve()\n'+
              var_str+
              'M =' + M_str +';\n'
              'b = ' + b_str + ';\n'
              'Z = null(full(M));\n'
              "Zstr = '';\n"
              'for col= 1:size(Z,2)\n'
              "Zstr = strcat(Zstr,';', char(Z(:,col)));\n"
              'end\n'
              "ret = strcat(char(linsolve(M, b)),'TBS', Zstr);")
    script_file = tempfile.NamedTemporaryFile(mode = 'wt',
                                              delete=False,
                                              dir = './',
                                              suffix=".m")
    checkstring = "Not set."
    script_file.write(script)
    script_file.close()
    del script
    linsolve_run = True
    trys = 0
    while linsolve_run:
        try:
            checkstring = getattr(eng, script_file.name.split('/')[-1][:-2])().split('TBS')
            solstring, nullstring = checkstring
            if solstring[0:6]=='matrix':
                sols_list = matlab_parser(solstring, matrix=True)
            else:
                sols_list = [matlab_parser(solstring)]
            if any(x in sols_list for x in [sympy.Symbol('Inf'), sympy.Symbol('Nan')]):
                trys = MAX_TRYS
                raise ValueError("Inconsistent matrix equation.")
        except Exception as e:
            trys+=1
            if trys < MAX_TRYS:
                try:
                    eng.quit()
                except Exception as equiterror:
                    print(str(equiterror))
                try:
                    eng = matlab.engine.start_matlab("-nojvm")
                except Exception as estarterror:
                    shutil.copy(script_file.name, './failed_matlab_script.m')
                    os.remove(script_file.name)
                    print(str(e))
                    print(checkstring)
                    return {}
            else:
                shutil.copy(script_file.name, './failed_matlab_script.m')
                os.remove(script_file.name)
                print(str(e))
                print(checkstring)
                return {}
        else:
            linsolve_run = False
            os.remove(script_file.name)
    try:
        eng.clearvars(nargout=0)
        eng.clear(script_file.name, nargout=0)
    except Exception as eclear:
        print(eclear)
    if len(nullstring) > 2:
        sepvecs = [matlab_parser(nullvecstr, matrix=True)
                           for nullvecstr in nullstring.split(';')[1:]]
        for nullvec in sepvecs:
            if len_oldfvars > 0:
                raise ValueError("Non empy iofvars not supported in Matlab!")
                newfvar = next(tempgen)
                iofvars.append(newfvar)
                tempvars.append(newfvar)
            else:
                newfvar = next(fvargen)
                newfvars.append(newfvar)
            sols_list = [sol+newfvar*nullvecel
                         for sol, nullvecel in zip(sols_list, nullvec)]
   #ipdb.set_trace()
    if not sols_list:
        return {}
    sols = dict(zip(fvars, sols_list))
    return {var: sol for var, sol in sols.items() if var is not sol}
