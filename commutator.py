from sympy import symbols, I, pretty, sympify, Matrix
from sympy.solvers.solveset import linsolve, linear_eq_to_matrix
from bisect import bisect_right
import ipdb
import json
import yaml
from collections import OrderedDict, defaultdict
import sympy
import re
from subprocess import check_output
from sympy.parsing.mathematica import mathematica

command='/home/kempj/bin/runMath'

def multiple_replace(string, rep_dict):
    pattern = re.compile("|".join([re.escape(k) for k in rep_dict.keys()]), re.M)
    return pattern.sub(lambda x: rep_dict[x.group(0)], string)

class Ncproduct:
    def __init__(self, scalar, product):
        self.scalar = scalar
        if isinstance(product, list):
            self.product = product
        elif isinstance(product, str):
            self.product = self.destringify(product)
        else:
            self.product = [product]

    def __getitem__(self, ind):
        return Ncproduct(self.scalar, self.product[ind])

    def get_unit(self, ind):
        return Ncproduct(1, self.product[ind])

    def get_operator(self,ind):
        return self.product[ind]

    def __setitem__(self, ind, value):
        self.product[ind] = value

    def is_product(self):
        return len(self.product) > 1

    def __add__(self, other):
        """Note add just returns a list, as this is addition"""
        if not isinstance(other,list):
            return [self, other]
        else:
            return [self]+other

    def __radd__(self, other):
        """Note add just returns a list, as this is addition"""
        if not isinstance(other,list):
            return [other, self]
        else:
            return other+[self]

    def __sub__(self, other):
        """Note add just returns a list, as this is addition"""
        if not isinstance(other,list):
            return [self, other*(-1)]
        else:
            return [self]+[a*(-1) for a in other]

    def __rsub__(self, other):
        """Note add just returns a list, as this is addition"""
        if not isinstance(other,list):
            return [other, self*(-1)]
        else:
            return other+[self*(-1)]

    def __mul__(self, other):
        if isinstance(other, Ncproduct):
            return Ncproduct(self.scalar*other.scalar, self.product+other.product)
        else:
            return Ncproduct(self.scalar*other, self.product)

    def __rmul__(self, other):
        if isinstance(other, Ncproduct):
            return Ncproduct(self.scalar*other.scalar, other.product+self.product)
        else:
            return Ncproduct(self.scalar*other, self.product)

    def __repr__(self):
        return str(self.scalar) + ' : ' +str(self.product)

    def __eq__(self, other):
        """Note == ignores the scalar field!"""
        try:
            return self.product == other.product
        except:
            return False

    def __len__(self):
        return len(self.product)

    def __str__(self):
        return '\u22c5'.join([str(self.scalar).replace('**','^').replace('*','\u22c5').replace('I','i')]
                       +[self.stringify(a) for a in self.product])

    def stringify(self, a):
        if a % 2 == 0:
            return 'b' + str(a//2)
        else:
            return 'a' + str((a+1)//2)

    def _texify_stringify(self, a):
        if a % 2 == 0:
            return 'b_' + str(a//2)
        else:
            return 'a_' + str((a+1)//2)

    def destringify(self, string):
        result = []
        string = string.replace('\u22c5', ' ')
        for op in string.split(' '):
            if op[0] == 'a':
                result.append(int(op[1:])*2-1)
            elif op[0] == 'b':
                result.append(int(op[1:])*2)
            else:
                print('Unknown operator ' + op)
        return result

    def texify(self):
        tex = sympy.latex(self.scalar)
        if (sympify(self.scalar).func ==sympy.Add):
            tex = '\\left (' + tex + '\\right )'
        return ' '.join([tex]
                       +[self._texify_stringify(a) for a in self.product])

def postmultiply(group,a):
     return [b*a for b in group]

def premultiply(a, group):
    return [a*b for b in group]

def commute(a,b):
    if a.is_product():
        return postmultiply(commute(a.get_unit(0),b), a[1:]) + premultiply(a.get_unit(0),commute(a[1:],b))
    elif b.is_product():
        return postmultiply(commute(a,b.get_unit(0)), b[1:]) + premultiply(b.get_unit(0),commute(a,b[1:]))
    elif a == b:
        return [0*a*b]
    else:
        return [2*a*b]

def commute_group(group_a, group_b):
    result = []
    for a in group_a:
        for b in group_b:
            result += commute(a,b)
    return result

def remove_zeros(group):
    group[:] = (a for a in group if a.scalar != 0)

def sort_anticommuting_product(ncprod):
    i = 0
    nflips = 0
    a = ncprod.product
    while i < len(a)-1:
        if a[i] > a[i+1]:
            a[i], a[i+1] = a[i+1],a[i]
            nflips += 1
            while i>0 and a[i] < a[i-1]:
                a[i], a[i-1] = a[i-1],a[i]
                nflips +=1
                i -= 1
        i+=1
    ncprod.scalar = ncprod.scalar * (-1)**nflips

def set_squares_to_identity(ncprod):
    i = 0
    a = ncprod.product
    while i < len(a)-1:
        if a[i] == a[i+1]:
            del a[i+1]
            del a[i]
        else:
            i+=1

def collect_terms(group):
    from collections import defaultdict
    D = defaultdict(list)
    for i,ncprod in enumerate(group):
        D[tuple(ncprod.product)].append(i)
    return [Ncproduct(sympy.expand(sum([group[i].scalar for i in D[key]])), list(key)) for key in D]

def full_collect_terms(group):
    from collections import defaultdict
    D = defaultdict(list)
    for i,ncprod in enumerate(group):
        D[tuple(ncprod.product)].append(i)
    return [Ncproduct(sympy.simplify(sum([group[i].scalar for i in D[key]])), list(key)) for key in D]

def simplify_group(group):
    remove_zeros(group)
    for ncprod in group:
        sort_anticommuting_product(ncprod)
        set_squares_to_identity(ncprod)
    group = collect_terms(group)
    remove_zeros(group)
    return group

def full_simplify_group(group):
    remove_zeros(group)
    for ncprod in group:
        sort_anticommuting_product(ncprod)
        set_squares_to_identity(ncprod)
    group = full_collect_terms(group)
    remove_zeros(group)
    return group

def multiply_groups(group_a, group_b):
    return simplify_group([a*b for a in group_a for b in group_b])

def square_group_to_order(group, order, split_orders):
    return simplify_group([(a*b)
                           for i in range(order+1)
                           for a in group[split_orders[i]:split_orders[i+1]]
                           for b in group[split_orders[(order-i)]:split_orders[(order-i+1)]]
                           ])

def square_group_to_order(group, order, split_orders):
    return simplify_group([(a*b)
                           for i in range(order+1)
                           for a in group[split_orders[i]:split_orders[i+1]]
                           for b in group[split_orders[(order-i)]:split_orders[(order-i+1)]]
                           ])

def square_to_find_identity(group):
    return simplify_group([a*a for a in collect_terms(group)])

def square_to_find_identity_scalar_up_to_order(group, order, split_orders):
    D = defaultdict(dict)
    for i,ncprod in enumerate(group):
        D[tuple(ncprod.product)].update({bisect_right(split_orders,i)-1: i})
    result = 0
    for torder in range(order+1):
        for product, positions in D.items():
            for iorder, index in positions.items():
                jorder = torder - iorder
                if jorder in positions:
                    result += (group[positions[iorder]].scalar
                               *group[positions[jorder]].scalar
                               *((2-(len(product)%4))))
    return sympy.expand(result)

def calculate_commutator(group_a,group_b):
    if not isinstance(group_a, list):
        group_a = [group_a]
    if not isinstance(group_b, list):
        group_b = [group_b]
    group = commute_group(group_a, group_b)
    return simplify_group(group)

def print_progress(i, length):
    print(str(i+1)+'/'+str(length), end = '\r')

def find_order(expr, orders):
    """Where order is power of small quantity, and orders a dict of
    symbols with their order.
    """
    expr = sympify(expr)
    order = float("inf")
    if expr.func == sympy.Add:
        for arg in expr.args:
            torder = find_order(arg, orders)
            if torder < order:
                order = torder
    elif expr.func == sympy.Mul:
        order = sum([find_order(arg, orders) for arg in expr.args])
    elif expr.func == sympy.Pow:
        order = find_order(expr.args[0], orders)*expr.args[1]
    elif expr.func == sympy.Symbol:
        order = orders.get(expr,0)
    else:
        order = 0
    return order

def neglect_to_order(expr, order, orders):
    expr = sympy.expand(expr)
    if find_order(expr, orders) > order:
        return 0
    if expr.func == sympy.Add:
        expr = sympy.Add(*[neglect_to_order(term, order, orders) for term in expr.args])
    return expr

def order_group(group, orders):
    return sorted(group, key = lambda a: find_order(a.scalar,orders))

def check_group_at_least_order(group, order, orders):
    for ncprod in group:
        if find_order(ncprod.scalar, orders) <= order:
            test = sympy.simplify(ncprod.scalar)
            if test != 0 and find_order(test, orders) <= order:
                print('Error: ' + str(test) + ': ' + str(ncprod.product))
                return False
    return True

def print_group(group, breaks = True):
    if not hasattr(print_group, 'orders'):
        print_group.orders = {}
    if isinstance(group, list):
        if len(group) > 1:
            group = order_group(group, print_group.orders)
            if breaks:
                print('{'+ ('+'+'\n+'.join(str(a) for a in group)).replace('+-', '\u2212')+'}')
            if not breaks:
                print(' + '.join(str(a) for a in group).replace('+ -', '\u2212'))
        elif group:
                print(group[0])
        else:
            print('0')
    else:
        print(group)

def texify_group(group):
    """Uses same orders as print_group"""
    if not hasattr(print_group, 'orders'):
        print_group.orders = {}
    if isinstance(group, list):
        if len(group) > 1:
            group = order_group(group, print_group.orders)
            return('$$'+' + '.join(a.texify() for a in group).replace('+ -', '-')+'$$')
        else:
            return('$$'+group[0].texify()+'$$')
    else:
        return('$$'+group.texify()+'$$')


def ordered_dump(data, stream=None, Dumper=yaml.Dumper, **kwds):
    class OrderedDumper(Dumper):
        pass
    def _dict_representer(dumper, data):
        return dumper.represent_mapping(
            yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG,
            data.items())
    OrderedDumper.add_representer(OrderedDict, _dict_representer)
    return yaml.dump(data, stream, OrderedDumper, **kwds)


def write_yaml(data, filename, **kwargs):
    """Writes yaml to file. Use py:class:`OrderedDict` to preserve key
    order.
    """
    if filename[-5:] != ".yaml":
        filename = filename + ".yaml"
    with open(filename, mode = 'w') as file_obj:
        ordered_dump(data, stream=file_obj, default_flow_style = False,
                     **kwargs)

def save_group(group, filename, iofvars=None, split_orders = None):
  #  ipdb.set_trace()
    if iofvars is None:
        iofvars = []
    if split_orders is None:
        split_orders = []
    data = OrderedDict([('group', group),
                        ('iofvars', iofvars),
                        ('split_orders', split_orders)])
    write_yaml(data, filename)

def load_group(filename, iofvars = None, split_orders = None):
    ext = '.yaml'
    if filename[-5:] != ext:
            filename = filename + ext
    with open(filename) as f:
        parsed = yaml.load(f)
    if iofvars is not None:
        iofvars[:] = parsed['iofvars']
    if split_orders is not None:
        split_orders[:] = parsed['split_orders']
    return parsed['group']

def substitute_group(group, subs_rules, split_orders = None):
    temp = [Ncproduct(sympify(ncprod.scalar).xreplace(subs_rules),
                     ncprod.product) for ncprod in group]
    remove_zeros(temp)
    if split_orders is not None:
        split_orders[-1] = len(temp)
    return temp


def fill_subspace_rows(to_cancel, matrixrows, subspace, Jpart):
    row_to_fill = matrixrows[subspace[tuple(to_cancel.product)]]
    if not row_to_fill:
        #pdb.set_trace()
        comm = calculate_commutator(Jpart, Ncproduct(1,to_cancel.product))
        row_to_fill[:] = [0]*len(subspace)
        for ncprod in comm:
           try:
                ind = subspace[tuple(ncprod.product)]
           except KeyError:
                ind = len(subspace)
                subspace[tuple(ncprod.product)] = ind
                for row in matrixrows:
                    if row:
                        row.append(0)
                matrixrows.append([])
                fill_subspace_rows(ncprod, matrixrows, subspace, Jpart)
           row_to_fill[ind] = ncprod.scalar

def find_subspace(to_cancel, Jpart):
    subspace = OrderedDict()
    matrixrows = []
    for i, ncprod in enumerate(to_cancel):
        if not tuple(ncprod.product) in subspace:
            subspace[tuple(ncprod.product)] = len(subspace)
            for row in matrixrows:
                if row:
                    row.append(0)
            matrixrows.append([])
            fill_subspace_rows(ncprod, matrixrows, subspace, Jpart)
    return subspace, matrixrows

def build_vector_to_cancel(to_cancel, subspace):
    cvector = [0]*len(subspace)
    for ncprod in to_cancel:
        cvector[subspace[tuple(ncprod.product)]] = ncprod.scalar
    return cvector



def print_subspace(subspace):
    for key, item in subspace.items():
        print(str(item)+ ': ' + ' '.join([Ncproduct.stringify(Ncproduct,a) for a in key]))

def sparse_fill_subspace_rows(to_cancel, matrixrows, subspace, Jpart, ind_col):
        #pdb.set_trace()
        comm = calculate_commutator(Jpart, Ncproduct(1,to_cancel.product))
        for ncprod in comm:
           try:
                ind_row = subspace[tuple(ncprod.product)]
           except KeyError:
                ind_row = len(subspace)
                subspace[tuple(ncprod.product)] = ind_row
                matrixrows[ind_row] = []
                sparse_fill_subspace_rows(ncprod, matrixrows, subspace, Jpart, ind_row)
           matrixrows[ind_row].append((ind_col, ncprod.scalar))

def sparse_find_subspace(to_cancel, Jpart):
    subspace = OrderedDict()
    matrixrows = {}
    length = len(to_cancel)
    for i, ncprod in enumerate(to_cancel):
        if not tuple(ncprod.product) in subspace:
            ind = len(subspace)
            subspace[tuple(ncprod.product)] = ind
            matrixrows[ind] = []
            sparse_fill_subspace_rows(ncprod, matrixrows, subspace, Jpart, ind)
        print_progress(i, length)
    return subspace, matrixrows

def sparse_normalise(psi_total, order, orders, coeffs, cvector, matrixrows, split_orders, start_ind = 0):
    norm = square_group_to_order(psi_total, order, split_orders)
    to_cancel = [ncprod.scalar for ncprod in norm if ncprod.product]
    ind = start_ind
    for term in to_cancel:
        term = sympy.expand(term)
        matrixrows[ind] = []
        row = matrixrows[ind]
        for ind_col, coeff in enumerate(coeffs):
            product = term.coeff(coeff)
            if product != 0:
                row.append((ind_col, product))
        const_term = term.as_coeff_add(*coeffs)[0]
        cvector.append(-const_term)
        ind+=1

def merge(lsts):
    """From stackoverflow: http://stackoverflow.com/questions/9110837/"""
    sets = [set(lst) for lst in lsts if lst]
    merged = 1
    while merged:
        merged = 0
        results = []
        while sets:
            common, rest = sets[0], sets[1:]
            sets = []
            for x in rest:
                if x.isdisjoint(common):
                    sets.append(x)
                else:
                    merged = 1
                    common |= x
            results.append(common)
        sets = results
    return sets

def find_sub_subspaces(matrixrows):
    return [list(space) for space in merge([[el[0] for el in row] for rownum, row in matrixrows.items()])]


def linear_solve(augmatrix, fvars, fvargen, newfvars, tempgen, tempvars, len_oldfvars):
    #return sympy.solve_linear_system(augmatrix,*fvars)
    mstr = multiple_replace(str(augmatrix)[7:-1],{ '[':'{',']':'}', '**':'^'})
    parameter = ('M = ' + mstr + ';'
                 'ToString['
                 '{'
                 'Simplify[LinearSolve[M[[1 ;; -1, 1 ;; -2]], M[[All, -1]]]], TBS,'
                 'NullSpace[M[[1 ;; -1, 1 ;; -2]]]'
                 '}'
                 ', InputForm]')
    solstring, nullstring = check_output([command,parameter])[2:-3].decode("utf-8").split(', TBS, ')
    #print(solstring)
    #print(check_output([command,parameter])[2:-3].decode("utf-8").split(', TBS, '))
    #ipdb.set_trace()
    try:
        sols_list =  [sympify(sol.replace('^', '**'))
                      for sol in solstring[:-1].split(',')]
    except Exception as e:
        print(str(e))
        return {}
    if len(nullstring) > 2:
        sepvecs = nullstring[1:-1].split('}, ')
        for nullvec in sepvecs:
            nullvec = nullvec[1:].split(',')
            if any(nullvecel != 0
                   for nullvecel in nullvec[len(nullvec)-len_oldfvars:]):
                newfvar = next(tempgen)
                iofvars.append(newfvar)
                tempvars.append(newfvar)
            else:
                newfvar = next(fvargen)
                newfvars.append(newfvar)
            sols_list = [sol+newfvar*sympify(nullvecel)
                         for sol, nullvecel in zip(sols_list, nullvec)]
    #ipdb.set_trace()
    if not sols_list:
        return {}
    sols = dict(zip(fvars, sols_list))
    return {var: sol for var, sol in sols.items() if var is not sol}

def solve_for_sub_subspace(matrixrows, sub_sub_space,
                           coeffs, cvector,
                           iofvars, subs_rules,
                           fvargen, newfvars, tempgen, tempvars):
    sspacedict = dict(zip(sub_sub_space, range(len(sub_sub_space))))
    length = len(sub_sub_space)
    augmatrixrows = []
    rownumstore = [] #for debugging
    for rownum, row in matrixrows.items():
        if row and row[0][0] in sub_sub_space:
            augmatrixrows.append(length*[0]+[cvector[rownum]])
            rownumstore.append(rownum)
            for el in row:
                augmatrixrows[-1][sspacedict[el[0]]] = el[1]
    fvars = [coeffs[ind] for ind in sub_sub_space]
    augmatrix = Matrix(augmatrixrows)
    oldfvars = []
    if iofvars:
        #ipdb.set_trace()
        atoms = augmatrix.atoms(sympy.Symbol)
        for iofvar in subs_rules:
            if iofvar in atoms:
                augmatrix = augmatrix.xreplace({iofvar: subs_rules[iofvar]})
        atoms = augmatrix.atoms(sympy.Symbol)
        for iofvar in iofvars:
            if iofvar not in subs_rules and iofvar in atoms:
                fvars.append(iofvar)
                oldfvars.append(iofvar)
                augmatrix = augmatrix.col_insert(-1, sympy.zeros(augmatrix.rows,1))
                for row_ind in range(len(augmatrix[:,0])):
                    coeff_val = -sympy.expand(augmatrix[row_ind,-1]).coeff(iofvar)
                    augmatrix[row_ind,-2] = coeff_val
                    augmatrix[row_ind,-1] += coeff_val*iofvar
    sols = linear_solve(augmatrix, fvars, fvargen, newfvars, tempgen, tempvars, len(oldfvars))
    if not sols:
        print(repr(augmatrix))
        print(fvars)
        print(rownumstore)
        print(iofvars)
        print(subs_rules)
        raise ValueError("Failure. No solutions.")
    for oldfvar in oldfvars:
        if oldfvar in sols:
            subs_rules.update({var: rule.xreplace({oldfvar: sols[oldfvar]})
                          for var, rule in subs_rules.items()})
            subs_rules[oldfvar] = sympy.simplify(sols[oldfvar])
    return sols


def sparse_solve_for_commuting_term(cvector, psi_lower, order, orders,
                                    matrixrows, subspace, norm = False,
                                    fvarname = 'A', iofvars = None, subs_rules = None, split_orders = None):
    fvar_gen = sympy.numbered_symbols('fvar')
    fvars = [next(fvar_gen) for i in range(len(subspace))]
    psi_order = [Ncproduct(-fvars[subspace[key]], list(key))
                 for i,key in enumerate(subspace)]
    if norm:
        psi_total = psi_lower + psi_order
        new_orders = orders.copy()
        new_orders.update(dict(zip(fvars, [order]*len(fvars))))
        sparse_normalise(psi_total, order, new_orders, fvars, cvector, matrixrows, split_orders, start_ind = len(fvars))
    sub_sub_spaces = find_sub_subspaces(matrixrows)
    print(sub_sub_spaces)
    solutions = {}
    if subs_rules is None:
        subs_rules = {}
    #deal with empty rows
    # for i, row in matrixrows.items():
    #     if not row:
    #         if sympy.simplify(cvector[i]) != 0:
    #             poss = False
    #             for iofvar in iofvars:
    #                 if iofvar in cvector[i].atoms(sympy.Symbol):
    #                     sub_sub_spaces.append([i])
    #                     print('Warning, new sspace: ' + str(i))
    #                     poss = True
    #             if not poss:
    #                 print(matrixrows)
    #                 print(cvector)
    #                 raise ValueError('Error term to cancel in null')
    length_ss = len(sub_sub_spaces)
    fvargen = sympy.numbered_symbols(fvarname)
    tempgen = sympy.numbered_symbols('temp')
    tempvars = []
    newfvars = []
    for i, ss_space in enumerate(sub_sub_spaces):
        #if i == 4:
        #    ipdb.set_trace()
        solutions.update(solve_for_sub_subspace(matrixrows, ss_space,
                                                fvars, cvector, iofvars,
                                                subs_rules, fvargen, newfvars, tempgen, tempvars))
        print_progress(i, length_ss)
    solvector = []
    for fvar in fvars:
        try:
            solvector.append(solutions[fvar])
        except KeyError:
            newfvars.append(next(fvargen))
            solvector.append(newfvars[-1])
    for tempvar in tempvars:
        if tempvar not in subs_rules:
            newfvars.append(next(fvargen))
            subs_rules[tempvar] = newfvars[-1]
    if newfvars:
        if iofvars is not None:
            iofvars[:] = newfvars
    if not subs_rules:
        return simplify_group([Ncproduct(-solvector[i], list(key))
                           for i,key in enumerate(subspace)])
    else:
        ret = simplify_group([Ncproduct(-solvector[i].xreplace(subs_rules), list(key))
                               for i,key in enumerate(subspace)])
        for tempvar in tempvars:
            subs_rules.pop(tempvar, None)
        return ret



def check_normalisable(psi, fvars, order, orders, split_orders, update_splits = True):
    matrixrows = {}
    cvector = []
    solutions = {}
    if update_splits:
        split_orders.append(len(psi))
    if not fvars:
        norm = square_group_to_order(psi, order, split_orders)
        to_cancel = [ncprod for ncprod in norm if ncprod.product]
        if to_cancel:
           raise ValueError('Non-normalisable: '+ str(to_cancel))
    sparse_normalise(psi, order, orders, fvars, cvector, matrixrows, split_orders)
    for i, row in matrixrows.items():
        if not row:
            if sympy.simplify(cvector[i]) != 0:
                raise ValueError('Non-normalisable: '+ str(cvector[i]))
    sub_sub_spaces = find_sub_subspaces(matrixrows)
    length_ss = len(sub_sub_spaces)
    for i, ss_space in enumerate(sub_sub_spaces):
        solutions.update(solve_for_sub_subspace(matrixrows, ss_space,
                                                fvars, cvector, None,
                                                None, None, None))
        print_progress(i, length_ss)
    return solutions
