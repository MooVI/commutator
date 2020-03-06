import os
import shutil
import tempfile
from subprocess import check_output
from tokenize import NAME, OP
from keyword import iskeyword


from sympy import zeros, Symbol
from sympy.parsing.sympy_parser import (parse_expr, convert_xor,
                                        auto_number)
from sympy.core.basic import Basic

from mathematica_printer import mstr


command = './linsolve_gauss'

def auto_symbol_add(tokens, local_dict, global_dict):
    """Inserts calls to ``Symbol``/``Function`` for undefined variables."""
    result = []
    prevTok = (None, None)

    tokens.append((None, None))  # so zip traverses all tokens
    for tok, nextTok in zip(tokens, tokens[1:]):
        tokNum, tokVal = tok
        nextTokNum, nextTokVal = nextTok
        if tokNum == NAME:
            name = tokVal

            if (name in ['True', 'False', 'None']
                or iskeyword(name)
                # Don't convert attribute access
                or (prevTok[0] == OP and prevTok[1] == '.')
                # Don't convert keyword arguments
                or (prevTok[0] == OP and prevTok[1] in ('(', ',')
                    and nextTokNum == OP and nextTokVal == '=')):
                result.append((NAME, name))
                continue
            elif name in local_dict:
                if isinstance(local_dict[name], Symbol) and nextTokVal == '(':
                    result.extend([(NAME, 'Function'),
                                   (OP, '('),
                                   (NAME, repr(str(local_dict[name]))),
                                   (OP, ')')])
                else:
                    result.append((NAME, name))
                continue
            elif name in global_dict:
                obj = global_dict[name]
                if isinstance(obj, (Basic, type)) or callable(obj):
                    result.append((NAME, name))
                    continue

            print("here")
            auto_symbol_add.newfvar = next(auto_symbol_add.fvargen)
            auto_symbol_add.newfvars.append(auto_symbol_add.newfvar)
            local_dict[name] = auto_symbol_add.newfvar
            result.append((NAME, name))
        else:
            result.append((tokNum, tokVal))

        prevTok = (tokNum, tokVal)

    return result


transformations = (convert_xor, auto_number, auto_symbol_add)

def linear_solve(sparse_mat_rep, sub_cvector, length, fvars, iofvars,
                 fvargen, newfvars, vardict):
    matrix = zeros(len(sub_cvector), length)
    for index, value in sparse_mat_rep.items():
        matrix[index[0], index[1]] = value
    if matrix.cols == 1 and matrix.rows == 1:
        return {fvars[0]: sub_cvector[0]/matrix[0, 0]}
    M_str = mstr(matrix)
    b_str = '{'+','.join('{'+mstr(el)+'}'for el in sub_cvector) + '}'
    script = M_str + '\n' + b_str
    script_file = tempfile.NamedTemporaryFile(mode='wt',
                                              delete=False,
                                              dir='./',
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
                                transformations=transformations,
                                local_dict=temp_vardict)
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
    print(sols)
    return {var: sol for var, sol in sols.items() if var is not sol}
