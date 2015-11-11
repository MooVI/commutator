from sympy import symbols, I, pretty, sympify, Matrix
from sympy.solvers.solveset import linsolve, linear_eq_to_matrix
import ipdb
from collections import OrderedDict
import sympy

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
        return self.product == other.product

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

def square_to_find_identity(group):
    return simplify_group([a*a for a in group])

def calculate_commutator(group_a,group_b):
    if not isinstance(group_a, list):
        group_a = [group_a]
    if not isinstance(group_b, list):
        group_b = [group_b]
    group = commute_group(group_a, group_b)
    return simplify_group(group)

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
        if find_order(ncprod.scalar, orders) < order:
            if sympy.simplify(ncprod.scalar) != 0:
                print('Error: ' + str(ncprod))
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

def substitute_group(group, subs_rules):
    temp = [Ncproduct(sympify(ncprod.scalar).subs((var, rule) for var, rule in subs_rules.items()),
		     ncprod.product) for ncprod in group]
    return full_simplify_group(temp)

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
    for ncprod in to_cancel:
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
    for ncprod in to_cancel:
        if not tuple(ncprod.product) in subspace:
            ind = len(subspace)
            subspace[tuple(ncprod.product)] = ind
            matrixrows[ind] = []
            sparse_fill_subspace_rows(ncprod, matrixrows, subspace, Jpart, ind)
    return subspace, matrixrows

def sparse_normalise(psi_total, order, orders, coeffs, cvector, matrixrows, start_ind = 0):
    norm = multiply_groups(psi_total, psi_total)
    to_zero = [ncprod.scalar for ncprod in norm if (find_order(ncprod.scalar, orders) <= order
                                                    and ncprod.product)]
    ind = start_ind
    for term in to_zero:
        term = sympy.expand(term)
        matrixrows[ind] = []
        row = matrixrows[ind]
        for ind_col, coeff in enumerate(coeffs):
            product = neglect_to_order(term.coeff(coeff), 0, orders)
            if product != 0:
                row.append((ind_col, product))
        const_term = term.as_coeff_add(*coeffs)[0]
        const_term = neglect_to_order(const_term, order, orders)
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


def linear_solve(augmatrix, fvars):
    #return sympy.solve_linear_system(augmatrix,*fvars)
    sols_set = linsolve(augmatrix, *fvars)
    if not sols_set:
        return {}
    sols = dict(zip(fvars, list(sols_set)[0]))
    return {var: sympy.simplify(sol) for var, sol in sols.items() if var is not sol}

def solve_for_sub_subspace(matrixrows, sub_sub_space, coeffs, cvector, iofvars, subs_rules, debug = False):
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
                augmatrix = augmatrix.subs(iofvar,subs_rules[iofvar])
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
    sols = linear_solve(augmatrix, fvars)
    if not sols:
        print(repr(augmatrix))
        print(fvars)
        print(rownumstore)
        print(iofvars)
        print(subs_rules)
        raise ValueError("Failure. No solutions.")
    for oldfvar in oldfvars:
        if oldfvar in sols:
            subs_rules.update({var: rule.subs(oldfvar, sols[oldfvar])
                          for var, rule in subs_rules.items()})
            subs_rules[oldfvar] = sympy.simplify(sols[oldfvar])
    return sols
    
                    
def sparse_solve_for_commuting_term(cvector, psi_lower, order, orders,
                                    matrixrows, subspace, norm = True,
                                    fvarname = 'A', iofvars = None, subs_rules = None):
    fvar_gen = sympy.numbered_symbols('fvar')
    fvars = [next(fvar_gen) for i in range(len(subspace))]
    psi_order = [Ncproduct(-fvars[subspace[key]], list(key))
                 for i,key in enumerate(subspace)]
    if norm:
        psi_total = psi_lower + psi_order
        new_orders = orders.copy()
        new_orders.update(dict(zip(fvars, [order]*len(fvars))))
        sparse_normalise(psi_total, order, new_orders, fvars, cvector, matrixrows, start_ind = len(fvars))
    sub_sub_spaces = find_sub_subspaces(matrixrows)
    print(sub_sub_spaces)
    solutions = {}
    if subs_rules is None:
        subs_rules = {}
    length_ss = len(sub_sub_spaces)
    for i, ss_space in enumerate(sub_sub_spaces):
        #if i == 4:
        #    ipdb.set_trace()
        solutions.update(solve_for_sub_subspace(matrixrows, ss_space,
                                                fvars, cvector, iofvars,
                                                subs_rules))
        print(str(i+1)+'/'+str(length_ss), end = '\r')
    solvector = []
    newfvars = []
    oldfvars  = []
    freevars_gen = sympy.numbered_symbols(fvarname)
    for fvar in fvars:
        try:
            solvector.append(solutions[fvar])
        except KeyError:
            newfvar = next(freevars_gen)
            solvector.append(newfvar)
            newfvars.append(newfvar)
            oldfvars.append(fvar)
            print(str(fvar)+': ' + str(newfvar))
    if newfvars:
        rules = [sub for sub in zip(oldfvars, newfvars)]
        solvector = [sol.subs(rules) for sol in solvector]
        if iofvars is not None:
            iofvars[:] = newfvars
    if not subs_rules:
        return simplify_group([Ncproduct(-solvector[i], list(key))
                           for i,key in enumerate(subspace)])
    else:
        return simplify_group([Ncproduct(-solvector[i].subs(
            [(var, rule) for var, rule in subs_rules.items()]), list(key))
                           for i,key in enumerate(subspace)])


def check_normalisable(psi, fvars, order, orders):
    matrixrows = {}
    cvector = []
    solutions = {}
    sparse_normalise(psi, order, orders, fvars, cvector, matrixrows)
    sub_sub_spaces = find_sub_subspaces(matrixrows)
    for ss_space in sub_sub_spaces:
        solutions.update(solve_for_sub_subspace(matrixrows, ss_space,
                                                fvars, cvector, None,
                                                None))
    return solutions
