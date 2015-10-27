from sympy import symbols, I, pretty, sympify, Matrix
from sympy.solvers.solveset import linsolve
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
        if (self.scalar.func ==sympy.Add):
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

def simplify_group(group):
    remove_zeros(group)
    for ncprod in group:
        sort_anticommuting_product(ncprod)
        set_squares_to_identity(ncprod)
    group = collect_terms(group)
    remove_zeros(group)
    return group

def multiply_groups(group_a, group_b):
    return simplify_group([a*b for a in group_a for b in group_b])

def calculate_commutator(group_a,group_b):
    if not isinstance(group_a, list):
        group_a = [group_a]
    if not isinstance(group_b, list):
        group_b = [group_b]
    group = commute_group(group_a, group_b)
    return simplify_group(group)

def find_order(expr,orders):
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

def order_group(group, orders):
    return sorted(group, key = lambda a: find_order(a.scalar,orders))

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
        else:
            print(group[0])
    else:
        print(group)

def texify_group(group):
    """Uses same orders as print_group"""
    if not hasattr(print_group, 'orders'):
        print_group.orders = {}
    if isinstance(group, list):
        if len(group) > 1:
            group = order_group(group, print_group.orders)
            print(' + '.join(a.texify() for a in group).replace('+ -', '-'))
        else:
            print(group[0].texify())
    else:
        print(group.texify())

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

def _normalise_fill_coeff(product, rows, cvector, coeffs):
    if product.func != sympy.Mul:
        raise ValueError
    for i,coeff in enumerate(coeffs):
        if coeff in product.args:
            product = sympy.Mul(*[prod for prod in product.args if prod is not coeff])
            rows[-1][i] = product
            break
    else:
        cvector[-1] = -product

def normalise(psi_total, order, orders, coeffs, cvector):
    norm = multiply_groups(psi_total, psi_total)
    to_zero = [ncprod.scalar for ncprod in norm if (find_order(ncprod.scalar, orders) <= order
                                                    and ncprod.product)]
    rows = []
    for term in to_zero:
        rows.append([0]*len(coeffs))
        cvector.append(0)
        if term.func == sympy.Add:
            for arg in term.args:
                if find_order(arg, orders) <= order:
                    _normalise_fill_coeff(arg, rows, cvector, coeffs)
        else:
            _normalise_fill_coeff(term, rows, cvector, coeffs)
    return rows
            
            
def solve_for_commuting_term(cvector, psi_lower, order, orders, matrixrows, subspace):
    matrix = Matrix(matrixrows)
    augmatrix = matrix.col_insert(len(subspace), Matrix(cvector))
    if matrix.rank() <= augmatrix.rank():
        if matrix.rank() == len(subspace):
            solvector = matrix.LUsolve(Matrix(cvector))
        else:
            rref = augmatrix.rref(simplify=True)
            iter = 0
            while iter < len(rref[0][:,0]):
                if all(a==0 for a in rref[0][iter,0:-1]):
                    rref[0].row_del(iter)
                else:
                    iter+=1
            cvector = [n for sn in rref[0][:,-1].tolist() for n in sn] #flatten list
            a = sympy.numbered_symbols('b')
            b = [next(a) for i in range(len(subspace))]
            psi_order = [Ncproduct(b[subspace[key]], list(key))for i,key in enumerate(subspace)]
            psi_total = psi_lower + psi_order
            new_orders = orders.copy()
            new_orders.update(dict(zip(b, [order]*len(b))))
            rows = normalise(psi_total, order, new_orders, b, cvector)
            matrix = Matrix([row for row
                             in [rref[0][i,0:-1] for i in range(len(rref[0][:,0]))]]
                             +rows)
            augmatrix = matrix.col_insert(len(matrix[0,:]), Matrix(cvector))
            solutions = sympy.solve_linear_system(augmatrix,*b)
            if not solutions:
                print('Failed. Inconsistency.')
                print(augmatrix)
                return None
            solvector = [solutions[b[i]] for i in range(len(b))]
    else:
        print('Matrix not invertible. Something has gone wrong.')
        return None
    return simplify_group([Ncproduct(solvector[i], list(key)) for i,key in enumerate(subspace)])

        
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

def _sparse_normalise_fill_coeff(product, matrixrows, cvector, coeffs, ind_row):
    if product.func != sympy.Mul:
        raise ValueError
    for ind_col,coeff in enumerate(coeffs):
        if coeff in product.args:
            product = sympy.Mul(*[prod for prod in product.args if prod is not coeff])
            matrixrows[ind_row].append((ind_col, product))
            break
    else:
        cvector[-1] = product

def sparse_normalise(psi_total, order, orders, coeffs, cvector, matrixrows):
    norm = multiply_groups(psi_total, psi_total)
    to_zero = [ncprod.scalar for ncprod in norm if (find_order(ncprod.scalar, orders) <= order
                                                    and ncprod.product)]
    ind = len(coeffs)
    for term in to_zero:
        matrixrows[ind] = []
        cvector.append(0)
        if term.func == sympy.Add:
            for arg in term.args:
                if find_order(arg, orders) <= order:
                    _sparse_normalise_fill_coeff(arg, matrixrows, cvector, coeffs, ind)
        else:
            _sparse_normalise_fill_coeff(term, matrixrows, cvector, coeffs, ind)
        ind += 1

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

def solve_for_sub_subspace(matrixrows, sub_sub_space, coeffs, cvector):
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
    sols = linsolve(Matrix(augmatrixrows),fvars)
    if not sols:
        print(repr(Matrix(augmatrixrows)))
        print(fvars)
        print(rownum)
        raise ValueError("Failure. No solutions.")
    return dict(zip(fvars, list(sols)[0]))
    
                    
def sparse_solve_for_commuting_term(cvector, psi_lower, order, orders, matrixrows, subspace):
    fvar_gen = sympy.numbered_symbols('fvar')
    fvars = [next(fvar_gen) for i in range(len(subspace))]
    psi_order = [Ncproduct(fvars[subspace[key]], list(key))for i,key in enumerate(subspace)]
    psi_total = psi_lower + psi_order
    new_orders = orders.copy()
    new_orders.update(dict(zip(fvars, [order]*len(fvars))))
    sparse_normalise(psi_total, order, new_orders, fvars, cvector, matrixrows)
    sub_sub_spaces = find_sub_subspaces(matrixrows)
    solutions = {}
    for ss_space in sub_sub_spaces:
        solutions.update(solve_for_sub_subspace(matrixrows, ss_space, fvars, cvector))
    solvector = [solutions[fvars[i]] for i in range(len(fvars))]
    return simplify_group([Ncproduct(-solvector[i], list(key)) for i,key in enumerate(subspace)])
