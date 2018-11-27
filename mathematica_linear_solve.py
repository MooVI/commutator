from sympy import symbols, I, pretty, sympify, Matrix, S
from mathematica_printer import mstr
import tempfile
import os
from subprocess import check_output
import shutil

command='/home/kempj/bin/MathematicaScript'

def mathematica_parser(exprstring):
    return sympify(exprstring.replace('^', '**'), mathematica_parser.vardict)



def linear_solve(sparse_mat_rep, sub_cvector, length, fvars, iofvars,
                 fvargen, newfvars, vardict):
    mathematica_parser.vardict = vardict
    if len(sparse_mat_rep) == 1:
        return {fvars[0]: sub_cvector[0]/sparse_mat_rep[(0,0)]}
    #return sympy.solve_linear_system(augmatrix,*fvars)
    sparse_str = ('SparseArray[{'
                  + ','.join(['{'+str(ind[0]+1)+','+str(ind[1]+1)+'}' for ind in sparse_mat_rep])
                  + '}->{'
                  + ','.join([mstr(value) for ind, value in sparse_mat_rep.items()])
                  + '}]')
    cvector_str = '{'+','.join([mstr(val) for val in sub_cvector]) + '}'
    script = ('Print['
              'M = ' + sparse_str + ';'
              'b = SparseArray[' + cvector_str + '];'
              'ToString['
              '{'
              'Simplify[LinearSolve[M, b]], TBS,'
              'NullSpace[M]'
              '}'
              ', InputForm] ];')
    script_file = tempfile.NamedTemporaryFile(mode = 'wt', delete=False)
    checkstring = "Not set."
    try:
        script_file.write(script)
        script_file.close()
        del script
        checkstring = check_output([command, '-script', script_file.name])[2:-3].decode("utf-8").split(', TBS, ')
        solstring, nullstring = checkstring
        #print(solstring)
        #print(check_output([command,parameter])[2:-3].decode("utf-8").split(', TBS, '))
        #ipdb.set_trace()
        sols_list =  [mathematica_parser(sol)
                      for sol in solstring[:-1].split(',')]
    except Exception as e:
        shutil.copy(script_file.name, './failed_script.m')
        print(str(e))
        print(checkstring)
        return {}
    finally:
        os.remove(script_file.name)
    if len(nullstring) > 2:
        sepvecs = nullstring[1:-1].split('}, ')
        for nullvec in sepvecs:
            nullvec = nullvec[1:].split(',')
            newfvar = next(fvargen)
            newfvars.append(newfvar)
            sols_list = [sol+newfvar*mathematica_parser(nullvecel)
                         for sol, nullvecel in zip(sols_list, nullvec)]
   #ipdb.set_trace()
    if not sols_list:
        return {}
    sols = dict(zip(fvars, sols_list))
    return {var: sol for var, sol in sols.items() if var is not sol}
