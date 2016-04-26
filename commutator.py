from sympy import symbols, I, pretty, sympify, Matrix, S
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
import tempfile
import os
from mathematica_printer import mstr
from matlab_printer import matlabstr
import shutil
import matlab.engine


eng = matlab.engine.start_matlab("-nojvm")

command='/home/kempj/bin/MathematicaScript'
qcommand = '/home/kempj/bin/runMath'

def multiple_replace(string, rep_dict):
    pattern = re.compile("|".join([re.escape(k) for k in rep_dict.keys()]), re.M)
    return pattern.sub(lambda x: rep_dict[x.group(0)], string)

def mathematica_parser(exprstring):
    return sympify(exprstring.replace('^', '**'), mathematica_parser.vardict)

def matlab_parser(exprstring, matrix = False):
    if matrix:
        return sympify('Matri'+exprstring[5:].replace('^','**').replace('i', '*I'), mathematica_parser.vardict)
    return sympify(exprstring.replace('^','**').replace('i', '*I'), mathematica_parser.vardict)


mathematica_parser.vardict = {}

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

    def is_identity(self):
        return len(self.product) == 0

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

    def conjugate(self):
        return Ncproduct(sympy.conjugate(self.scalar)*(2-len(self.product)%4), self.product)

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
        if len(string) > 0:
            string = string.replace('\u22c5', ' ')
            for op in string.split(' '):
                if op[0] == 'a':
                    result.append(int(op[1:])*2-1)
                elif op[0] == 'b':
                    result.append(int(op[1:])*2)
                else:
                    print('Unknown operator ' + op)
            return result
        else:
            return []

    def texify(self):
        tex = sympy.latex(self.scalar)
        if (sympify(self.scalar).func ==sympy.Add):
            tex = '\\left (' + tex + '\\right )'
        return ' '.join([tex]
                       +[self._texify_stringify(a) for a in self.product])


class SigmaProduct(Ncproduct):
    def stringify(self, a):
        if a % 2 == 0:
            return 'x' + str(a//2)
        else:
            return 'z' + str((a+1)//2)

    def _texify_stringify(self, a):
        if a % 2 == 0:
            return '\\sigma^x_{' + str(a//2)+'}'
        else:
            return '\\sigma^z_{' + str((a+1)//2)+'}'

    def destringify(self, string):
        result = []
        string = string.replace('\u22c5', ' ')
        for op in string.split(' '):
            if op[0] == 'z':
                result.append(int(op[1:])*2-1)
            elif op[0] == 'x':
                result.append(int(op[1:])*2)
            else:
                print('Unknown operator ' + op)
        return result




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
            if not (a.is_identity() or b.is_identity()):
                result += commute(a,b)
    return result


def commute_up_to_order(group_a, group_b, order, split_orders_a, split_orders_b):
    result = []
    a_order_max = len(split_orders_a)-1
    b_order_max = len(split_orders_b)-1
    if a_order_max + b_order_max <= order:
        return commute_group(group_a, group_b)
    for aorder in range(a_order_max):
        for border in range(min([order-aorder, b_order_max])):
            result += commute_group(group_a[split_orders[aorder]:split_orders[aorder+1]],
                                    group_b[split_orders[border]:split_orders[border+1]])
    return result

def remove_zeros(group):
    group[:] = (a for a in group if a.scalar != 0)

def conjugate_group(group):
    return [a.conjugate() for a in group]


def sort_pauli_product(ncprod):
    i = 0
    nflips = 0
    a = ncprod.product
    while i < len(a)-1:
        diff = a[i] - a[i+1]
        if diff > 0:
            a[i], a[i+1] = a[i+1], a[i]
            if diff == 1:
                if a[i] % 2 == 1:
                    nflips += 1
            while i>0 and a[i] < a[i-1]:
                a[i], a[i-1] = a[i-1], a[i]
                if a[i] - a[i-1] == 1:
                    if a[i-1] % 2 == 1:
                        nflips +=1
                i -= 1
        i+=1
    ncprod.scalar = ncprod.scalar * (-1)**nflips

def sort_anticommuting_product(ncprod):
    i = 0
    nflips = 0
    a = ncprod.product
    while i < len(a)-1:
        if a[i] > a[i+1]:
            a[i], a[i+1] = a[i+1],a[i]
            nflips += 1
            while i>0 and a[i] < a[i-1]:
                a[i], a[i-1] = a[i-1], a[i]
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

def convert_to_sigma(ncprod):
    #ipdb.set_trace()
    ret = SigmaProduct(1, [])
    ret.scalar = ncprod.scalar
    for el in ncprod.product:
        if el % 2 == 0:
            ret.scalar *= I
            ret.product += [i for i in range(2,el,2)]+ [el-1, el]
        else:
            ret.product += [i for i in range(2,el,2)]+[el]
    sort_pauli_product(ret)
    set_squares_to_identity(ret)
    return ret


def convert_from_sigma(sigma):
    ret = Ncproduct(1, [])
    ret.scalar = sigma.scalar
    for el in sigma.product:
        if el % 2 == 0:
            ret.scalar *= -I
            ret.product += [el-1, el]
        else:
            ret.product += [i for i in range(1,el,1)]+[el]
            ret.scalar *= (-I)**((el-1)//2)
    sort_anticommuting_product(ret)
    set_squares_to_identity(ret)
    return ret

def convert_group(group):
    if not isinstance(group, list):
        group = [group]
    if isinstance(group[0], SigmaProduct):
        return [convert_from_sigma(el) for el in group]
    elif isinstance(group[0], Ncproduct):
        return [convert_to_sigma(el) for el in group]
    else:
        raise ValueError('Unrecognised conversion asked for!')

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


def mathematica_simplify(scalar, use_tempfile = False):
    sstr = mstr(scalar)
    parameter = ('ToString['
                 'Simplify[' + sstr +']'
                 ', InputForm]')
    simpstring = "Not set."
    if use_tempfile:
        script_file = tempfile.NamedTemporaryFile(mode = 'wt', delete=False)
        script_file.write("Print["+parameter+"]")
        script_file.close()
        simpstring = check_output([command, '-script', script_file.name])[0:-1].decode("utf-8")
        os.remove(script_file.name)
    else:
        simpstring = check_output([qcommand,parameter])[0:-1].decode("utf-8")
    try:
        return mathematica_parser(simpstring)
    except Exception as e:
        print(str(e))
        print(simpstring)

def mathematica_series(scalar, varlist, order):
    sstr = mstr(scalar)
    varsstr = ','.join(['{'+str(var) +', 0, ' + str(order) +'}' for var in varlist])
    parameter = ('ToString[' +
                 'Series['+
                 sstr+
                 ',' +varsstr+']'
                 ', InputForm]')
    simpstring = check_output([qcommand,parameter])[0:-1].decode("utf-8")
    return mathematica_parser(simpstring)

def math_collect_terms(group, use_tempfile = False):
    from collections import defaultdict
    D = defaultdict(list)
    for i,ncprod in enumerate(group):
        D[tuple(ncprod.product)].append(i)
    return [Ncproduct(mathematica_simplify(sum([group[i].scalar for i in D[key]]), use_tempfile), list(key)) for key in D]

def mathematica_simplify_group(group, use_tempfile = False):
    remove_zeros(group)
    for ncprod in group:
        sort_anticommuting_product(ncprod)
        set_squares_to_identity(ncprod)
    group = math_collect_terms(group, use_tempfile)
    remove_zeros(group)
    return group

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

def multiply_conjugate_groups(group_a, group_b):
    return simplify_group([a.conjugate()*b for a in group_a for b in group_b])

def square_group_to_order(group, order, split_orders):
    """Assumes Hermitian"""
    return simplify_group([(a*b)
                           for i in range(order+1)
                           for a in group[split_orders[i]:split_orders[i+1]]
                           for b in group[split_orders[(order-i)]:split_orders[(order-i+1)]]
                           ])

def square_to_find_identity(group):
    """Assumes Hermitian"""
    return simplify_group([a*a for a in collect_terms(group)])

def square_to_find_identity_scalar_up_to_order(group, order, split_orders):
    """Assumes Hermitian and operator strings odd in length"""
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
                               *(2-len(product)%4)*group[positions[jorder]].scalar)
    return sympy.expand(result)

def square_to_find_identity_scalar_up_to_order_poss_even(group, order, split_orders):
    """Assumes Hermitian and operator strings odd in length"""
    D = defaultdict(dict)
    for i,ncprod in enumerate(group):
        D[tuple(ncprod.product)].update({bisect_right(split_orders,i)-1: i})
    result = 0
    for torder in range(order+1):
        for product, positions in D.items():
            for iorder, index in positions.items():
                jorder = torder - iorder
                if jorder in positions:
                    length = len(product)
                    if length % 2 == 0:
                        length +=1
                    result += (group[positions[iorder]].scalar
                               *(2-length%4)*group[positions[jorder]].scalar)
    return sympy.expand(result)

def calculate_commutator(group_a,group_b):
    if not isinstance(group_a, list):
        group_a = [group_a]
    if not isinstance(group_b, list):
        group_b = [group_b]
    group = commute_group(group_a, group_b)
    return simplify_group(group)

def commute_with_perturbing_H(Hpert, group, split_orders):
    if not isinstance(group, list):
        group = [group]
    result = []
    for a in Hpert:
        for b in group[split_orders[-2]:split_orders[-1]]:
            if not (a.is_identity() or b.is_identity()):
                result += commute(a,b)
    return simplify_group(result)

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
        return S.Zero
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

def print_group(group, breaks = True, neglect_order = None):
    if not hasattr(print_group, 'orders'):
        print_group.orders = {}
    if isinstance(group, list):
        if len(group) > 1:
            group = order_group(group, print_group.orders)
            if neglect_order is not None:
                for ncprod in group:
                    ncprod.scalar  = neglect_to_order(ncprod.scalar, neglect_order, print_group.orders)
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

def texify_group(group, newlines = False):
    """Uses same orders as print_group"""
    if not hasattr(print_group, 'orders'):
        print_group.orders = {}
    if isinstance(group, list):
        if len(group) > 1:
            if not newlines:
                group = order_group(group, print_group.orders)
                return('$$'+' + '.join(a.texify() for a in group).replace('+ -', '-').replace('\\\\','\\')+'$$')
            else:
                group = order_group(group, print_group.orders)
                return('\\begin{align*}\n'+'\\\\\n'.join('&' + a.texify().replace('\\\\','\\') for a in group)+'\n\\end{align*}')
        else:
            return('$$'+group[0].texify().replace('\\\\','\\')+'$$')
    else:
        return('$$'+group.texify().replace('\\\\','\\')+'$$')


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

def save_group(group, filename, iofvars=None, split_orders = None, normdict = None):
  #  ipdb.set_trace()
    if iofvars is None:
        iofvars = []
    if split_orders is None:
        split_orders = []
    if normdict is None:
        normdict = []
    data = OrderedDict([('group', group),
                        ('iofvars', iofvars),
                        ('split_orders', split_orders),
                        ('normdict', normdict)])
    write_yaml(data, filename)

def load_group(filename, iofvars = None, split_orders = None, normdict = None):
    ext = '.yaml'
    if filename[-5:] != ext:
            filename = filename + ext
    with open(filename) as f:
        parsed = yaml.load(f)
    if iofvars is not None:
        iofvars[:] = parsed['iofvars']
    if split_orders is not None:
        split_orders[:] = parsed['split_orders']
    if normdict is not None:
        normdict.update(parsed['normdict'])
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
        ind_cols = [ind_col]
        to_cancels = [to_cancel]
        count = 0
        while count < len(ind_cols):
            comm = calculate_commutator(Jpart, Ncproduct(1,to_cancels[count].product))
            for ncprod in comm:
                try:
                    ind_row = subspace[tuple(ncprod.product)]
                except KeyError:
                    ind_row = len(subspace)
                    subspace[tuple(ncprod.product)] = ind_row
                    matrixrows[ind_row] = []
                    ind_cols.append(ind_row)
                    to_cancels.append(ncprod)
                matrixrows[ind_row].append((ind_cols[count], ncprod.scalar))
            count+=1

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

def sparse_normalise(psi_total, order, orders, coeffs, cvector, matrixrows, split_orders, start_ind = 0, subspace = None):
    if subspace is None:
        subspace = []
    norm = square_group_to_order(psi_total, order, split_orders)
    to_cancel = [ncprod for ncprod in norm if ncprod.product]
    ind = start_ind
    for ncprod in to_cancel:
        term = sympy.expand(ncprod.scalar)
        matrixrows[ind] = []
        row = matrixrows[ind]
        for ind_col, coeff in enumerate(coeffs):
            product = term.coeff(coeff)
            if product != 0:
                row.append((ind_col, product))
        const_term = term.as_coeff_add(*coeffs)[0]
        cvector.append(-const_term)
        subspace.append(Ncproduct(-const_term, ncprod.product))
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


def linear_solve(augmatrix, fvars, iofvars, fvargen, newfvars, tempgen, tempvars, len_oldfvars):
    global eng
    MAX_TRYS = 2
    #return sympy.solve_linear_system(augmatrix,*fvars)
    if augmatrix.cols == 2 and augmatrix.rows == 1:
        return {fvars[0]: augmatrix[0,1]/augmatrix[0,0]}
    M_str = matlabstr(augmatrix[:,:-1])
    b_str = matlabstr(augmatrix[:,-1])
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
    sols = linear_solve(augmatrix, fvars, iofvars, fvargen, newfvars, tempgen, tempvars, len(oldfvars))
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
        #ipdb.set_trace()
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



def _not_edge_normalise(psi0, to_cancel):
    #ipdb.set_trace()
    szero = set(psi0.product)
    scancel  = set(to_cancel.product)
    intersect = (scancel & szero)
    newproduct = sorted(list((scancel | szero) - intersect))
    length = len(szero)-len(intersect)
    if length % 2 == 0:
        length += 1
        print('Warning: Non-normalisable with odd number of operators: ' + str(to_cancel))
    return Ncproduct((length%4-2)*to_cancel.scalar/(2*psi0.scalar), newproduct)


def check_normalisable(psi,
                       fvars,
                       order,
                       orders,
                       split_orders,
                       zero_not_needed = False,
                       update_splits = True,
                       make_norm = True,
                       not_edge = False,
                       simplify = sympy.simplify,
                       Jpart = None):
    matrixrows = {}
    cvector = []
    solutions = {}
    if update_splits:
        split_orders.append(len(psi))
    if not fvars:
        norm = square_group_to_order(psi, order, split_orders)
        to_cancel = [ncprod for ncprod in norm if ncprod.product]
        for ncprod in to_cancel:
            if simplify(ncprod.scalar) != 0:
                if make_norm and ncprod.product[0] != 1 and not not_edge:
                    psi += [Ncproduct(-ncprod.scalar/2,[1]+ncprod.product)]
                    print("Adding " + str(psi[-1]) )
                    split_orders[-1] += 1
                elif (make_norm
                      and not_edge):
                      #and ncprod.product[:len(psi[0].product)] != psi[0].product):
                    psi += [_not_edge_normalise(psi[0], ncprod)]
                    print("Adding " + str(psi[-1]) )
                    split_orders[-1] += 1
                else:
                    raise ValueError('Non-normalisable: '+ str(to_cancel))
                if Jpart and calculate_commutator(Jpart, Ncproduct(1, psi[-1].product)):
                    raise ValueError('Non-normalisable: '+ str(to_cancel)
                                     + ".\n Canceling term " + str(psi[-1])
                                     + " does not commute.")
    subspace = []
    sparse_normalise(psi, order, orders, fvars, cvector, matrixrows, split_orders, subspace=subspace)
    for i, row in matrixrows.items():
        if not row:
            if simplify(cvector[i]) != 0:
                if make_norm and subspace[i].product[0] != 1:
                    psi += [Ncproduct(subspace[i].scalar/2,[1]+subspace[i].product)]
                    print("Adding " + str(psi[-1]) )
                    split_orders[-1] += 1
                elif (make_norm
                      and not_edge
                      and subspace[i].product[:len(psi[0].product)] != psi[0].product):
                    psi += [_not_edge_normalise(psi[0], -subspace[i])]
                    print("Adding " + str(psi[-1]) )
                    split_orders[-1] += 1
                else:
                    raise ValueError('Non-normalisable: '+ str(subspace[i]))
                if Jpart and calculate_commutator(Jpart, Ncproduct(1, psi[-1].product)):
                    raise ValueError('Non-normalisable: '+ str(subspace[i])
                                     + ".\n Canceling term " + str(psi[-1])
                                     + " does not commute.")
    sub_sub_spaces = find_sub_subspaces(matrixrows)
    length_ss = len(sub_sub_spaces)
    for i, ss_space in enumerate(sub_sub_spaces):
        solutions.update(solve_for_sub_subspace(matrixrows, ss_space,
                                                fvars, cvector, None,
                                                None, None, None, None, None))
        print_progress(i, length_ss)
    if zero_not_needed:
        for fvar in fvars:
            if fvar not in solutions:
                solutions[fvar] = 0
    return solutions

def truncate_fill(to_cancel, coeffs, matrixrows, cvector, start_ind = 0):
    ind = start_ind
    to_cancel = [el.scalar for el in to_cancel]
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

def check_truncate(to_cancel, fvars):
    matrixrows = {}
    cvector = []
    solutions = {}
    if not fvars:
        if to_cancel:
           return False
    truncate_fill(to_cancel, fvars, matrixrows, cvector)
    for i, row in matrixrows.items():
        if not row:
            if sympy.simplify(cvector[i]) != 0:
                return False
    sub_sub_spaces = find_sub_subspaces(matrixrows)
    length_ss = len(sub_sub_spaces)
    for i, ss_space in enumerate(sub_sub_spaces):
        try:
            solutions.update(solve_for_sub_subspace(matrixrows, ss_space,
                                                fvars, cvector, None,
                                                None, None, None, None, None))
        except Exception as e:
            print(str(e))
            return False
        print_progress(i, length_ss)
    return solutions

def _build_entire_subspace(result, former, start, end):
    for i in range(start, end):
        term = former + [i]
        result[:] += [term]
        _build_entire_subspace(result, term, i+1, end)

def build_entire_subspace(L):
    start = 1
    end = L*2
    result = []
    former = []
    _build_entire_subspace(result, former, start, end)
    return [Ncproduct(1, el) for el in result]

def build_norm_subspace(L):
    start = 2
    end = L*2
    result = []
    former = [1]
    _build_entire_subspace(result, former, start, end)
    return [Ncproduct(1, el) for el in result]

def build_odd_subspace(L):
    return [ncprod for ncprod in build_entire_subspace(L) if len(ncprod.product) % 2 == 1]

def build_odd_norm_subspace(L):
    return [ncprod for ncprod in build_norm_subspace(L) if len(ncprod.product) % 2 == 1]

def solve_at_once(H, L, iofvars = None):
    subspace_ops = build_odd_subspace(L)
    len_subspace = len(subspace_ops)
    subspace = OrderedDict(zip([tuple(el.product) for el in subspace_ops],
                                [i for i in range(len_subspace)]))
    matrixrows = {}
    for i in range(len_subspace):
        matrixrows[i] = []
    for ind_col in range(len_subspace):
        comm = calculate_commutator(H, subspace_ops[ind_col])
        for ncprod in comm:
            ind_row = subspace[tuple(ncprod.product)]
            matrixrows[ind_row].append((ind_col, ncprod.scalar))
    cvector = [0]*len_subspace
    if iofvars is None:
        iofvars = []
    return sparse_solve_for_commuting_term(cvector,
                                       None,
                                       None,
                                       None,
                                       matrixrows,
                                       subspace,
                                       norm = False,
                                       fvarname = 'F',
                                       iofvars=iofvars)

def fill_subspace(Jpart, order):
    L = 2*order+1
    subspace_ops = build_odd_norm_subspace(L)
    len_subspace = len(subspace_ops)
    matrixrows = defaultdict(list)
    subspace = OrderedDict(zip([tuple(el.product) for el in subspace_ops],
                                [i for i in range(len_subspace)]))
    for product, ind_col in subspace.items():
        comm = calculate_commutator(Jpart, Ncproduct(1, list(product)))
        for ncprod in comm:
           ind_row = subspace[tuple(ncprod.product)]
           matrixrows[ind_row].append((ind_col, ncprod.scalar))
    return subspace, matrixrows

def fill_subspace_norm(Jpart, to_cancel, order):
    L = 2*order+1
    subspace_ops = build_odd_norm_subspace(L)
    return sparse_find_subspace(to_cancel+subspace_ops, Jpart)

def get_yaml_poly(poly, gens):
    s = sympy.Poly(poly, gens)
    pows = [[float(e) for e in el] for el in s.monoms()]
    coeffs = [float(sympy.N(coeff)) for coeff in s.coeffs()]
    return {'pows':pows,'coeffs':coeffs}


def save_to_c_poly(group,
                   couplings,
                   filename,
                   norm_file = None,
                   split_orders = None,
                   order = None):
    if not norm_file:
        if split_orders:
            if not order:
                order = len(split_orders)-1
            norm = square_to_find_identity_scalar_up_to_order(group, order, split_orders)
        else:
            norm = square_to_find_identity(group)[0].scalar
    data = OrderedDict([('Psi',[]), ('Norm', {'poly': get_yaml_poly(norm, couplings)})])
    for ncprod in reversed(group):
        xpos = []
        zpos = []
        sigprod = convert_to_sigma(ncprod)
        for a in sigprod.product:
            if a % 2 == 0:
                xpos.append((a//2)-1)
            else:
                zpos.append(((a+1)//2)-1)
        data['Psi'].append({'xpos': xpos,
                              'zpos': zpos,
                              'poly': get_yaml_poly(sigprod.scalar, couplings)})
    write_yaml(data, filename)
