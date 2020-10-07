import tempfile
import os
from collections import OrderedDict, defaultdict
from itertools import chain, combinations
import re
from sympy import  I, sympify, S, factorial
from sympy.printing.ccode import ccode
from bisect import bisect_right
#import yaml
import pickle
import sympy
from subprocess import check_output

from mathematica_printer import mstr

from choose_linear_solve import linear_solve
from choose_n_h_linear_solve import n_h_linear_solve
h_linear_solve = linear_solve

# Bash commands to run mathematica script or server.
# Only need to be set if using Mathematica to simplify.
command='/home/kempj/bin/MathematicaScript'
qcommand = '/home/kempj/bin/runMath'

# Some list and permutation utility functions

def powerset(iterable):
    """From https://stackoverflow.com/questions/464864/how-to-get-all-possible-combinations-of-a-list-s-elements"""
    s = list(iterable)  # allows duplicate elements
    return chain.from_iterable(combinations(s, r) for r in range(len(s)+1))

def multiple_replace(string, rep_dict):
    pattern = re.compile("|".join([re.escape(k) for k in rep_dict.keys()]), re.M)
    return pattern.sub(lambda x: rep_dict[x.group(0)], string)


def next_permutationS(l):
    '''Changes a list to its next permutation, in place. From
    http://stackoverflow.com/questions/6534430/why-does-pythons-iterto
    ols-permutations-contain-duplicates-when-the-original
    Returns true unless wrapped around so result is lexicographically smaller.
    '''
    n = len(l)
    #Step 1: Find tail
    last = n-1 #tail is from `last` to end
    while last>0:
        if l[last-1] < l[last]: break
        last -= 1
    #Step 2: Increase the number just before tail
    if last>0:
        small = l[last-1]
        big = n-1
        while l[big] <= small: big -= 1
        l[last-1], l[big] = l[big], small
    #Step 3: Reverse tail
    i = last
    j = n-1
    while i < j:
        l[i], l[j] = l[j], l[i]
        i += 1
        j -= 1
    return last>0

def accel_asc(n):
    """From http://jeromekelleher.net/generating-integer-partitions.html
    Generates integer partions
    """
    a = [0 for i in range(n + 1)]
    k = 1
    y = n - 1
    while k != 0:
        x = a[k - 1] + 1
        k -= 1
        while 2 * x <= y:
            a[k] = x
            y -= x
            k += 1
        l = k + 1
        while x <= y:
            a[k] = x
            a[l] = y
            yield a[:k + 2]
            x += 1
            y -= 1
        a[k] = x + y
        y = x + y - 1
        yield a[:k + 1]


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

#YAML and other utility functions

def print_progress(i, length):
    print(str(i+1)+'/'+str(length), end = '\r')

# def ordered_dump(data, stream=None, Dumper=yaml.Dumper, **kwds):
#     class OrderedDumper(Dumper):
#         pass
#     def _dict_representer(dumper, data):
#         return dumper.represent_mapping(
#             yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG,
#             data.items())
#     OrderedDumper.add_representer(OrderedDict, _dict_representer)
#     return yaml.dump(data, stream, OrderedDumper, **kwds)


# def write_yaml(data, filename, **kwargs):
#     """Writes yaml to file. Use py:class:`OrderedDict` to preserve key
#     order.
#     """
#     if filename[-5:] != ".yaml":
#         filename = filename + ".yaml"
#     with open(filename, mode = 'w') as file_obj:
#         ordered_dump(data, stream=file_obj, default_flow_style = False,
#                      **kwargs)

# Mathematica Parsing Functions

def mathematica_parser(exprstring):
    return sympify(exprstring.replace('^', '**'), mathematica_parser.vardict)

mathematica_parser.vardict = {}

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
    definition = """multiTaylor[f_, {vars_?VectorQ, pt_?VectorQ, n_Integer?NonNegative}] :=
Sum[Nest[(vars - pt).# &, (D[f, {vars, \[FormalK]}] /.
Thread[vars -> pt]), \[FormalK]]/\[FormalK]!, {\[FormalK], 0,
n}, Method -> "Procedural"]""".replace('\n','')
    sstr = mstr(scalar)
    varsstr = ('{'
               +'{'+','.join([str(var) for var in varlist])+'},'
               +'{'+','.join([str(0)]*len(varlist))+'},'
               + str(order+1)+'}')
    parameter = (definition+';ToString[' +
                 'Expand[multiTaylor['+
                 sstr+
                 ',' +varsstr+']]'
                 ', InputForm]')
    simpstring = check_output([qcommand,parameter])[0:-1].decode("utf-8")
    return mathematica_parser(simpstring)



# Classes for Majorana and Sigma Products and Utility functions

class NCProduct:
    """Class to represents eiter a Pauli or Majorana string.
    product = the string itself, represented as a list of integers.
    scalar = a scalar multipying the string, a numeric or Sympy expression.

    For construction either:
        - give product in list of integers directly,
        - give space-separated string e.g. Paulis "x1 z2 y6", Majoranas
          "a1 b4 a5."
        - give tuple of string and list of lattice site numbers, e.g.
          ("xyz", [1, 2, 6]), which will be expanded to "x1 z2 y6".

    Assume 1-indexing.
    """

    def __init__(self, scalar, product):
        self.scalar = scalar
        if isinstance(product, list):
            self.product = product
        elif isinstance(product, tuple):
            prodstr = ' '.join(c +str(n) for c, n in zip(product[0], product[1]))
            self.product, scalar = self.destringify(prodstr)
            self.scalar = self.scalar * scalar
        elif isinstance(product, str):
            self.product, scalar = self.destringify(product)
            self.scalar = self.scalar * scalar
        else:
            self.product = [product]

    def is_product(self):
        return len(self.product) > 1

    def is_identity(self):
        return len(self.product) == 0

    def get_operator(self, ind):
        return type(self)(self.scalar, self.product[ind])

    def get_unit(self):
        return type(self)(1, self.product)

    def hermitian_conjugate(self):
        return type(self)(sympy.conjugate(self.scalar)*self.reverse_sign(), self.product)

    def conjugate(self):
        return type(self)(sympy.conjugate(self.scalar), self.product)

    def __getitem__(self, ind):
        return self.product[ind]

    def __setitem__(self, ind, value):
        self.product[ind] = value

    def __len__(self):
        return len(self.product)

    def __add__(self, other):
        """Note add just returns a NCSum, as this is addition"""
        if isinstance(other,NCProduct):
            return simplify_group([self, other])
        elif isinstance(other,NCSum):
            return simplify_group([self]+other.group)
        raise TypeError('Not an NCSum or NCProduct added to NCSum')


    def __radd__(self, other):
        """Note add just returns a NCSum, as this is addition"""
        if isinstance(other,NCProduct):
            return simplify_group([other, self])
        elif isinstance(other,NCSum):
            return simplify_group(other.group+[self])
        raise TypeError('Not an NCSum or NCProduct added to NCSum')

    def __sub__(self, other):
        """Note add just returns a NCSum, as this is addition"""
        if isinstance(other,NCProduct):
            return simplify_group([self, -other])
        elif isinstance(other,NCSum):
            return simplify_group([self]+[-a for a in other])
        raise TypeError('Not an NCSum or NCProduct added to NCSum')


    def __rsub__(self, other):
        """Note add just returns a NCSum, as this is addition"""
        if isinstance(other,NCProduct):
            return simplify_group([other, -self])
        elif isinstance(other,NCSum):
            return simplify_group(other.group+[-self])
        raise TypeError('Not an NCSum or NCProduct added to NCSum')

    def __mul__(self, other):
        typ = type(self)
        if isinstance(other, typ):
            return typ(self.scalar*other.scalar, self.product+other.product)
        elif isinstance(other, NCSum):
            return multiply_groups(NCSum(self), other)
        else:
            return typ(self.scalar*other, self.product)

    def __rmul__(self, other):
        typ = type(self)
        if isinstance(other, typ):
            return typ(self.scalar*other.scalar, other.product+self.product)
        elif isinstance(other, NCSum):
            return multiply_groups(other, NCSum(self))
        else:
            return typ(self.scalar*other, self.product)

    def __neg__(self):
        return type(self)(self.scalar*(-1), self.product)

    def __repr__(self):
        return str(self)
        #return str(self.scalar) + ' : ' +str(self.product)

    def __eq__(self, other):
        """Is scalar and product equal?"""
        try:
            return self.product == other.product and self.scalar == other.scalar
        except:
            return False

    def __str__(self):
        strings, scalar = self.self_stringify()
        return '\u202f'.join([str(self.scalar*scalar).replace('**','^').replace('*','\u202f').replace('I','i')]
                             + strings)

    def translate(self, nSites):
        db = nSites*2
        return type(self)(self.scalar, [el+db for el in self.product])

def sort_anticommuting_list(a):
    i = 0
    nflips = 0
    l = len(a)-1
    while i < l:
        if a[i] > a[i+1]:
            a[i], a[i+1] = a[i+1],a[i]
            nflips += 1
            while i>0 and a[i] < a[i-1]:
                a[i], a[i-1] = a[i-1], a[i]
                nflips +=1
                i -= 1
        i+=1
    return (-1)**nflips


class MajoranaProduct(NCProduct):
    """ Majoranas on site n are labelled a_n and b_n and have integer
    representations 2*n -1, 2*n
    """
    @staticmethod
    def stringify(a):
        if a % 2 == 0:
            return 'b' + str(a//2)
        else:
            return 'a' + str((a+1)//2)

    def self_stringify(self):
        return [self.stringify(a) for a in self.product], 1

    @staticmethod
    def _texify_stringify(a, init_string = None):
        if init_string is None:
            if a % 2 == 0:
                return 'b_' + str(a//2)
            else:
                return 'a_' + str((a+1)//2)
        else:
            if a % 2 == 0:
                return 'b_{'+init_string + ('+'+str(a//2 - 1) if a != 2 else '') +'}'
            else:
                return 'a_{'+init_string+('+'+ str((a+1)//2 - 1) if a != 1 else '')+'}'

    @staticmethod
    def destringify(string):
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
            return result, 1
        else:
            return [], 1

    @staticmethod
    def reverse_sign_product(product):
        """ The sign in the equation: product = +- (reverse-ordered product)
        for a pre-SORTED, NON-DUPLICATE product.
        """
        return 2*( len(product)%4 < 2) - 1

    def reverse_sign(self):
        return  2*( len(self.product)%4 < 2) - 1

    def texify(self, init_string = None):
        tex = sympy.latex(self.scalar)
        if (sympify(self.scalar).func ==sympy.Add):
            tex = '\\left (' + tex + '\\right )'
        return ' '.join([tex]
                        +[self._texify_stringify(a, init_string) for a in self.product])

    def sort(self):
        self.scalar *= sort_anticommuting_list(self.product)

    def commute(self, right):
        """Commutes so that self.commute(right) = [self, right]"""
        total = self.product + right.product
        rev =  right.product + self.product
        sign = sort_anticommuting_list(total)- sort_anticommuting_list(rev)
        return MajoranaProduct(sign*self.scalar*right.scalar, total)

#Legacy compatibility.
Ncproduct = MajoranaProduct

def sort_pauli_list(a):
    i = 0
    nflips = 0
    l = len(a) - 1
    while i < l:
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
    return (-1)**nflips

def count_ys_in_pauli_list(a):
    l = len(a) - 1
    n = 0
    i = 0
    while i < l:
        n += (a[i+1] == a[i] + 1)
        i += 2
    return n

class SigmaProduct(NCProduct):
    """ Paulis on site n are labelled x_n, y_n, z_n and have integer
    representations z_n =  2*n-1, x_n = 2*n and y_n = i*[2n*-1, 2*n].
    i.e. in the product y_n is given by z_n x_n and the i goes into the
    scalar, which keeps the product effectively real.
    """
    @staticmethod
    def stringify(a):
        if a % 2 == 0:
            return 'x' + str(a//2)
        else:
            return 'z' + str((a+1)//2)

    def self_stringify(self):
        "If product is not sorted, may not perform all y conversions."
        strings = []
        scalar = 1
        i = 0
        while i  < len(self):
            a = self.product[i]
            a2 = self.product[i+1] if i != len(self) -1 else None
            if a % 2 == 0:
                strings.append('x' + str(a//2))
            elif a2 == a + 1:
                strings.append('y' + str((a+1)//2))
                scalar = scalar * I
                i = i + 1
            else:
                strings.append('z' + str((a+1)//2))
            i = i +1
        return strings, scalar

    @staticmethod
    def _texify_stringify(a, init_string = None):
        if init_string is None:
            if a % 2 == 0:
                return '\\sigma^x_{' + str(a//2)+'}'
            else:
                return '\\sigma^z_{' + str((a+1)//2)+'}'
        else:
            if a % 2 == 0:
                return '\\sigma^x_{' +init_string + ('+'+str(a//2 -1) if a != 2 else '')+'}'
            else:
                return '\\sigma^z_{' +init_string + ('+'+str((a+1)//2 -1) if a != 1 else '')+'}'

    @staticmethod
    def destringify(string):
        result = []
        string = string.replace('\u22c5', ' ')
        scalar = 1
        for op in string.split(' '):
            if op[0] == 'z':
                result.append(int(op[1:])*2-1)
            elif op[0] == 'x':
                result.append(int(op[1:])*2)
            elif op[0] == 'y':
                result.append(int(op[1:])*2-1)
                result.append(int(op[1:])*2)
                scalar = scalar * (-I)
            else:
                print('Unknown operator ' + op)
        return result, scalar

    def sort(self):
        self.scalar *= sort_pauli_list(self.product)


    @staticmethod
    def reverse_sign_product(product):
        """ The sign in the equation: product = +- (reverse-ordered product)
        for a pre-SORTED, NON-DUPLICATE product.
        """
        return 1-2*(count_ys_in_pauli_list(product)%2)

    def reverse_sign(self):
        return 1-2*(count_ys_in_pauli_list(self.product)%2)

    def commute(self, right):
        """Commutes so that self.commute(right) = [self, right]"""
        total = self.product + right.product
        rev =  right.product + self.product
        sign = sort_pauli_list(total)- sort_pauli_list(rev)
        return SigmaProduct(sign*self.scalar*right.scalar, total)


def set_squares_to_identity(ncprod):
    i = 0
    a = ncprod.product
    while i < len(a)-1:
        if a[i] == a[i+1]:
            del a[i+1]
            del a[i]
        else:
            i+=1

def convert_to_sigma(majprod):
    ret = SigmaProduct(1, [])
    ret.scalar = majprod.scalar
    for el in majprod.product:
        if el % 2 == 0:
            ret.scalar *= I
            ret.product += [i for i in range(2,el,2)]+ [el-1, el]
        else:
            ret.product += [i for i in range(2,el,2)]+[el]
    ret.sort()
    set_squares_to_identity(ret)
    return ret


def convert_from_sigma(sigma):
    ret = MajoranaProduct(1, [])
    ret.scalar = sigma.scalar
    for el in sigma.product:
        if el % 2 == 0:
            ret.scalar *= -I
            ret.product += [el-1, el]
        else:
            ret.product += [i for i in range(1,el,1)]+[el]
            ret.scalar *= (-I)**((el-1)//2)
    ret.sort()
    set_squares_to_identity(ret)
    return ret


# Generally, to sum NCProducts just put them in a list.

# This class attempts to implement a translationally invariant sum.
# That is, everything inside is assummed repeated indefinitely.
# PLEASE IGNORE FTM.
class TransInvSum:
    """NOT TESTED Can be either MajoranaProduct or SigmaProduct, despite code names"""
    def __getitem__(self, ind):
        return self.ncprods[ind]

    def __init__(self, ncprods_list, check_base = False):
        """If width is given, you are vouching the sum is correct,
        check_base will be effectively False.
        """
        if isinstance(ncprods_list, TransInvSum):
            #Copy constructor
            self.ncprods = ncprods_list.ncprods[:]
        else:
            if isinstance(ncprods_list, MajoranaProduct) or isinstance(ncprods_list, SigmaProduct):
                ncprods_list = [ncprods_list]
            self.ncprods = ncprods_list
            if check_base:
                self.rebase()

    def get_width(self):
        return max((ncprod.product[-1] for ncprod in self.ncprods))

    def rebase(self):
        """Assumes products are ordered."""
        for ncprod in self.ncprods:
            if ncprod[0] != 1 and ncprod[0] != 2:
                ind_diff = ((ncprod[0]+1)//2) - 1
                ncprod[:] = [el - 2*ind_diff for el in ncprod]

    def __add__(self, other):
        """Adds translational inv. intelligently, otherwise assumes
        you just want to add to sum.
        """
        if isinstance(other, TransInvSum):
            return TransInvSum(self.ncprods + other.ncprods)
        else:
            return self+TransInvSum(other, check_base = True)

    def __radd__(self, other):
        """Adds translational inv. intelligently, otherwise assumes
        you just want to add to sum.
        """
        if isinstance(other, TransInvSum):
            return TransInvSum(other.ncprods+self.ncprods)
        else:
             return TransInvSum(other, check_base=True)+self

    def __sub__(self, other):
        """Adds translational inv. intelligently, otherwise assumes
        you just want to add to sum.
        """
        if isinstance(other, TransInvSum):
            return self + (-1)*other
        else:
            return self + (-1)*TransInvSum(other, check_base=True)

    def __rsub__(self, other):
        """Adds translational inv. intelligently, otherwise assumes
        you just want to add to sum.
        """
        if isinstance(other, TransInvSum):
            return other + (-1)*self
        else:
            return TransInvSum(other, check_base=True)+(-1)*self

    def __rmul__(self, other):
        return TransInvSum([other*ncprod for ncprod in self.ncprods])

    def __mul__(self, other):
        return TransInvSum([ncprod*other for ncprod in self.ncprods])

    def __repr__(self):
        return '\n'.join(repr(ncprod) for ncprod in self.ncprods)

    def conjugate(self):
        return TransInvSum([ncprod.conjugate() for ncprod in self.nprods])

    def __len__(self):
        return len(self.ncprods)

    def __str__(self):
        return ('\n'.join([str(ncprod) for ncprod in self.ncprods]))

    def texify(self, orders=None):
        return('\\begin{align*}\n \\sum_j '
               +'\\\\\n'.join(' &' + a.texify('j').replace('\\\\','\\')
                              for a in self.ncprods)+'\n\\end{align*}')
    def convert(self):
        newncprods = []
        if isinstance(self.ncprods[0], MajoranaProduct):
            for ncprod in self.ncprods:
                if len(ncprod.product) % 2 != 0:
                    raise ValueError("Non-even number of Majorana fermions"
                                     "will lead to trailing strings in spin basis,"
                                     "which is not supported by convert.")
                newncprods.append(convert_to_sigma(ncprod))
            return TransInvSum(newncprods, check_base = True)
        elif isinstance(self.ncprods[0], SigmaProduct):
            for ncprod in self.ncprods:
                if len(sum(el%2 for el in ncprod)) % 2 != 0:
                    raise ValueError("Non-even number of \\sigma^z operators"
                                     "will lead to trailing strings in fermion basis,"
                                     "which is not supported by convert.")
                newncprods.append(convert_from_sigma(ncprod))
            return TransInvSum(newncprods, check_base=True)
        else:
            raise ValueError('Unrecognised conversion asked for!')

    def simplify(self):
        """Returns self for legacy reasons with compatability with simplify_group"""
        remove_zeros(self.ncprods)
        for ncprod in self.ncprods:
            ncprod.sort()
            set_squares_to_identity(ncprod)
        self.rebase()
        self.ncprods = collect_terms(self.ncprods)
        remove_zeros(self.ncprods)
        return self


def commute_group_semi_inv(group_a, sumB, sign = 1):
    result = []
    prodtype = type(group_a[0])
    for a in group_a:
        for b in sumB:
            if not (a.is_identity() or b.is_identity()):
                maxB = (b[-1]+1)//2
                start_a = (a[0]+1)//2
                max_a = (a[-1]+1)//2
                for j in range(start_a-maxB, max_a):
                    result.append(a.commute(prodtype(sign*b.scalar, [el+j*2 for el in b])))
    return result

def commute_group_inv(sumA, sumB):
    result = []
    prodtype = type(sumA[0])
    for a in sumA:
        for b in sumB:
            if not (a.is_identity() or b.is_identity()):
                maxA = (a[-1]+1)//2
                maxB = (b[-1]+1)//2
                if maxA <= maxB:
                    fixed = prodtype(b.scalar, [el+2*(maxA-1) for el in b])
                    sign = 1
                    var = a
                if maxA > maxB:
                    fixed = prodtype(a.scalar, [el+2*(maxB-1) for el in a])
                    sign = -1
                    var = b
                for j in range(maxA+maxB-1):
                    result.append(prodtype(sign*var.scalar, [el+j*2 for el in var]).commute(fixed))
    return TransInvSum(result)



# Functions for dealing with real-space distances and supports.

def get_lattice_index(a):
    return (a+1)//2

def get_lattice_distance(a1,a2, L_periodic = None):
    if a1 ==a2:
        return 0
    n1 = get_lattice_index(a1)
    n2 = get_lattice_index(a2)
    n1, n2 = (n1, n2) if n1 < n2 else (n2, n1)
    return n2-n1 if not L_periodic else min([n2-n1, n1+L_periodic-n2])

def get_lattice_distance_only_wrapped(a1,a2, L_periodic):
    if a1 == a2:
        return 0
    n1 = get_lattice_index(a1)
    n2 = get_lattice_index(a2)
    n1, n2 = (n1, n2) if n1 < n2 else (n2, n1)
    return n1+L_periodic-n2

def get_support(ncproduct, L_periodic = None):
    product = ncproduct.product
    product = sorted(product)
    if len(product) < 2:
        return len(product)
    if not L_periodic:
        return get_lattice_distance(product[0], product[-1])+1
    site_indices = list(dict.fromkeys([(a+1)//2 for a in product]))
    return 1 + min([L_periodic+site_indices[i]-site_indices[i+1]
                   for i in range(len(site_indices)-1)]+
                   [site_indices[-1]-site_indices[0]])

class ind_converter:
    def __init__(self, row, col):
        self.row = row
        self.col = col
    def get_ind(self, x, y):
        return x*self.col + y
    def get_xy(self, ind):
        y = ind % self.col
        x = (ind - y)/self.col
        return x, y




# Class NCSum
# representing **sums** of NCProducts. An instance of an NCSum is called
# 'group' in the code.

# After any operation which might zero elements of a sum,
# or add duplicates which need to be merged simplify_group should be called.
# All of the group utility functions (operators +,-, * or functions
# multiply_groups, add_groups, etc.) do this automatically.
# If you wish to add groups without this .add_no_simplify can be used.

# full_simplify_group calls sympy.simplfy on each element. This is **slow**
# and should only be done when necesssary. mathematica_simplify_group does the
# same for mathematica.


class NCSum:
    def __init__(self, group):
        if isinstance(group, list):
            self.group = group
        elif isinstance(group, NCSum):
            self = group
        elif isinstance(group, NCProduct):
            self.group = [group]
        elif group == 0:
            self.group = []
        else:
            TypeError(str(type(group)) + " not convertable to NCSum")

    def __getitem__(self, ind):
        return self.group[ind]

    def __setitem__(self, ind, value):
        self.group[ind] = value

    def __len__(self):
        return len(self.group)

    def __iter__(self):
        return iter(self.group)

    def extend(self, group):
        if isinstance(group, NCSum):
            self.group += group.group
        else:
            self.group += group

    def add_no_simplify(self, group):
        self.group += group.group

    def is_zero(self):
        return len(self.group) == 0

    def __add__(self, other):
        if isinstance(other, NCSum):
            return simplify_group(self.group + other.group)
        elif isinstance(other, NCProduct):
            return simplify_group(self.group + [other])
        elif other == 0:
            return self
        else:
            raise TypeError('Not an NCSum or NCProduct')

    def __radd__(self, other):
        if isinstance(other, NCSum):
            return simplify_group(other.group + self.group)
        elif isinstance(other, NCProduct):
            return simplify_group([other] + self.group)
        elif other == 0:
            return self
        else:
            raise TypeError('Not an NCSum or NCProduct added to NCSum')

    def __sub__(self, other):
        if isinstance(other, NCSum):
            return simplify_group(self.group + [-a for a in other])
        elif isinstance(other, NCProduct):
            return simplify_group(self.group + [-other])
        elif other == 0:
            return self
        else:
            raise TypeError('Not an NCSum or NCProduct added to NCSum')

    def __rsub__(self, other):
        if isinstance(other, NCSum):
            return simplify_group(other.group + [-a for a in self])
        elif isinstance(other, NCProduct):
            return simplify_group([other] + [-a for a in self])
        elif other == 0:
            return -self
        else:
            raise TypeError('Not an NCSum or NCProduct added to NCSum')

    def __neg__(self):
        return (-1)*self

    def __mul__(self, other):
        if isinstance(other, NCSum):
            return multiply_groups(self, other)
        elif isinstance(other, NCProduct):
            return multiply_groups(self, NCSum(other))
        else:
            return NCSum([a*other for a in self])

    def __rmul__(self, other):
        if isinstance(other, NCSum):
            return multiply_groups(other, self)
        elif isinstance(other, NCProduct):
            return multiply_groups(NCSum(other), self)
        else:
            return NCSum([a*other for a in self])

    def __repr__(self):
        return repr(self.group)

    def __str__(self):
        if self._is_zero_():
            return '0'
        return (' + '.join([str(ncprod) for ncprod in self]))

    def translate(self, nSites):
        return NCSum([a.translate(nSites) for a in self])

    def texify(self, orders=None):
        return('\\begin{align*}\n \\sum_j '
               +'\\\\\n'.join(' &' + a.texify('j').replace('\\\\','\\')
                              for a in self)+'\n\\end{align*}')

    def __eq__(self, other):
        if isinstance(other, NCSum):
            return (self-other).is_zero()
        elif isinstance(other, NCProduct):
            return (self-NCSum(other)).is_zero()
        if other == 0:
            return self.is_zero()
        return False


# Simplification of groups

def collect_terms(group):
    from collections import defaultdict
    D = defaultdict(list)
    for i,ncprod in enumerate(group):
        D[tuple(ncprod.product)].append(i)
    return NCSum([type(group[0])(
        sympy.expand(sum([group[i].scalar for i in D[key]])), list(key))
                       for key in D])

def full_collect_terms(group):
    from collections import defaultdict
    D = defaultdict(list)
    for i,ncprod in enumerate(group):
        D[tuple(ncprod.product)].append(i)
    return NCSum([type(group[0])(
        sympy.simplify(sum([group[i].scalar for i in D[key]])), list(key))
                                for key in D])

def math_collect_terms(group):
    from collections import defaultdict
    D = defaultdict(list)
    for i,ncprod in enumerate(group):
        D[tuple(ncprod.product)].append(i)
    return NCSum([type(group[0])(
        mathematica_simplify(sum([group[i].scalar for i in D[key]])), list(key))
                       for key in D])

def remove_zeros(group):
    group[:] = (a for a in group if a.scalar != 0)

def simplify_remove_zeros(group):
    for ncprod in group:
        ncprod.scalar = sympy.simplify(ncprod.scalar)
    group[:] = (a for a in group if a.scalar != 0)

def simplify_group(group):
    remove_zeros(group)
    for ncprod in group:
        ncprod.sort()
        set_squares_to_identity(ncprod)
    group = collect_terms(group)
    remove_zeros(group)
    return group

def full_simplify_group(group):
    remove_zeros(group)
    for ncprod in group:
        ncprod.sort()
        set_squares_to_identity(ncprod)
    group = full_collect_terms(group)
    remove_zeros(group)
    return group

def mathematica_simplify_group(group, use_tempfile = False):
    remove_zeros(group)
    for ncprod in group:
        ncprod.sort()
        set_squares_to_identity(ncprod)
    group = math_collect_terms(group, use_tempfile)
    remove_zeros(group)
    return group



# Basic Group Mathematical Functions--------------------

def conjugate_group(group):
    return NCSum([a.conjugate() for a in group])

def hermitian_conjugate_group(group):
    return NCSum([a.hermitian_conjugate() for a in group])

#Legacy functions, use syntatic sugar instead (i.e just +,*,-)-------

def postmultiply(group,a):
     return group*a

def premultiply(a, group):
    return a*group

def multiply_groups(group_a, group_b):
    return simplify_group([a*b for a in group_a for b in group_b])

def add_groups(group_a, group_b):
    return simplify_group(group_a+group_b)

def subtract_groups(group_a, group_b):
    return simplify_group(group_a+premultiply(-1, group_b))

#End legacy functions---------------------------------

def commute_group(group_a, group_b):
    """For known groups. Prefer using calculate_commutator instead if
    could be individual NCProducts as well.
    """
    result = []
    for a in group_a:
        for b in group_b:
            if not (a.is_identity() or b.is_identity()):
                result.append(a.commute(b))
    return NCSum(result)

def calculate_commutator(group_a, group_b):
    """Works with single NCProducts and TransInvSum as well"""
    if len(group_a) == 0 or len(group_b) == 0:
        return type(group_a)([])
    if isinstance(group_a, TransInvSum):
        if isinstance(group_b, TransInvSum):
            group = commute_group_inv(group_a, group_b)
        else:
            if not isinstance(group_b, list):
                group_b = [group_b]
            group = commute_group_semi_inv(group_a, group_b, sign =-1)
    elif isinstance(group_b, TransInvSum):
        if not isinstance(group_a, list):
                group_a = [group_a]
        group = commute_group_semi_inv(group_a, group_b, sign =1)
    else:
        if isinstance(group_a, NCProduct):
                group_a = [group_a]
        if isinstance(group_b, NCProduct):
                group_b = [group_b]
        group = commute_group(group_a, group_b)
    return simplify_group(group)


def square_to_find_identity(group):
    return simplify_group([a*a for a in collect_terms(group)])

def square_to_find_identity_scalar(group):
    result = 0
    for ncprod in collect_terms(group):
        result += (sympy.expand(ncprod.scalar*ncprod.scalar))*ncprod.reverse_sign()
    return result

def trace_inner_product(group_a, group_b):
    """Tr(A B) not Tr(A^dagger B)"""
    from collections import defaultdict
    DA = defaultdict(list)
    DB = defaultdict(list)
    result = 0
    for i,ncprod in enumerate(group_a):
        DA[tuple(ncprod.product)].append(i)
    for i,ncprod in enumerate(group_b):
        DB[tuple(ncprod.product)].append(i)
    ind_pairs = [ [i, j] for key in set(DA.keys()).intersection(DB.keys())
                  for i in DA[key] for j in DB[key]]
    for i, j in ind_pairs:
        result += (group_a[i].scalar
                   *group_a[i].reverse_sign()*group_b[j].scalar)
    return result

def hermitian_trace_inner_product(group_a, group_b):
    """Tr(A^dagger B)"""
    from collections import defaultdict
    DA = defaultdict(list)
    DB = defaultdict(list)
    result = 0
    for i,ncprod in enumerate(group_a):
        DA[tuple(ncprod.product)].append(i)
    for i,ncprod in enumerate(group_b):
        DB[tuple(ncprod.product)].append(i)
    ind_pairs = [ [i, j] for key in set(DA.keys()).intersection(DB.keys())
                  for i in DA[key] for j in DB[key]]
    for i, j in ind_pairs:
        result += (sympy.conjugate(group_a[i].scalar)
                   *group_b[j].scalar)
    return result

def trace_hermitian_square(group):
    result = 0
    for ncprod in collect_terms(group):
        result += (sympy.expand(sympy.conjugate(ncprod.scalar)*ncprod.scalar))
    return result




# Order functions for perturbation theory.

# orders is a dict: keys = Sympy symbols; values = integer order in pertubation theory

# split_orders is a list of intervals associated with a group
# such that elements in the interval [split_orders[n], split_orders[n+1])
# have order n in pertubation theory.

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
    """Sort group by orders in perturbation theory"""
    return NCSum(sorted(group, key = lambda a: find_order(a.scalar,orders)))

def check_group_at_least_order(group, order, orders):
    """Check that group only has elements of order order or higher.
    Use to check that perturbation theory is being satisfied.
    """
    for ncprod in group:
        if find_order(ncprod.scalar, orders) < order:
            test = sympy.simplify(ncprod.scalar)
            if test != 0 and find_order(test, orders) < order:
                print('Error: ' + str(test) + ': ' + str(ncprod.product))
                return False
    return True

def commute_up_to_order(group_a, group_b, order, split_orders_a, split_orders_b):
    result = []
    a_order_max = len(split_orders_a)-1
    b_order_max = len(split_orders_b)-1
    if a_order_max + b_order_max <= order:
        return simplify_group(commute_group(group_a, group_b))
    for aorder in range(a_order_max):
        for border in range(min([order-aorder, b_order_max])):
            result += commute_group(group_a[split_orders_a[aorder]:split_orders_a[aorder+1]],
                                    group_b[split_orders_b[border]:split_orders_b[border+1]])
    return simplify_group(result)

def square_group_to_order(group, order, split_orders):
    return simplify_group([(a*b)
                           for i in range(order+1)
                           for a in group[split_orders[i]:split_orders[i+1]]
                           for b in group[split_orders[(order-i)]:split_orders[(order-i+1)]]
                           ])

def square_to_find_identity_scalar_up_to_order(group, order, split_orders):
    D = defaultdict(dict)
    typ = type(group[0])
    for i,ncprod in enumerate(group):
        D[tuple(ncprod.product)].update({bisect_right(split_orders,i)-1: i})
    result = 0
    for torder in range(order+1):
        for product, positions in D.items():
            for iorder, index in positions.items():
                jorder = torder - iorder
                if jorder in positions:
                    result += (group[positions[iorder]].scalar
                               *typ.reverse_sign_product(product)*group[positions[jorder]].scalar)
    return sympy.expand(result)


def unitary_transform(to_trans, Gs, max_order, inverse = False):
    """Unitary transform of form U = e^i(Gs[0]+Gs[1]+Gs[2]+...),
    where Gs[n] has order n+1. max_order is max_order of output,
    inverse = true calculates -i rather than +i.
    """
    result = NCSum(to_trans)
    for torder in range(1,max_order+1):
        for orderset in accel_asc(torder):
            numcomms = len(orderset)
            if inverse:
                taylorscalar = ((-I)**numcomms)/factorial(numcomms)
            else:
                taylorscalar = (I**numcomms)/factorial(numcomms)
            while True:
                cumul = NCSum(to_trans)
                for order in orderset:
                    cumul = calculate_commutator(cumul, Gs[order-1])
                result.extend(taylorscalar*cumul)
                if not next_permutationS(orderset):
                    break
    return simplify_group(result)

def unitary_transform_to_order(to_trans, Gs, torder,
                               not_single_comm = False,
                               inverse = False):
    """Unitary transform of form U = e^i(Gs[0]+Gs[1]+Gs[2]+...)
    at order torder, where Gs[n] has order n+1. Assumes to_trans is zeroth order.
    inverse = true calculates with a - sign. not_single_comm = true avoids using
    Gs[torder-1] for when this is not known yet.
    """
    if torder == 0:
        return NCSum(to_trans)
    result = NCSum([])
    for orderset in accel_asc(torder):
        numcomms = len(orderset)
        if numcomms != 1 or not not_single_comm:
            if inverse:
                taylorscalar = ((-I)**numcomms)/factorial(numcomms)
            else:
                taylorscalar = (I**numcomms)/factorial(numcomms)
            while True:
                cumul = NCSum(to_trans)
                for order in orderset:
                    cumul = calculate_commutator(cumul, Gs[order-1])
                result.extend(taylorscalar*cumul)
                if not next_permutationS(orderset):
                    break
    return simplify_group(result)

def exponentiate_to_order(Gs, torder, inverse = False):
    """Calculate U = e^i(Gs[0]+Gs[1]+Gs[2]+...), see above."""
    if torder == 0:
        return NCSum([type(Gs[0][0])(1,[])])
    result = NCSum([])
    for orderset in accel_asc(torder):
        numcomms = len(orderset)
        if inverse:
            taylorscalar = ((-I)**numcomms)/factorial(numcomms)
        else:
            taylorscalar = (I**numcomms)/factorial(numcomms)
        while True:
            cumul = [1]
            for order in orderset:
                cumul = multiply_groups(cumul, Gs[order-1])
            result.extend(taylorscalar*cumul)
            if not next_permutationS(orderset):
                break
    return simplify_group(result)



# Functions for dealing with free variables in groups.

def substitute_group(group, subs_rules, split_orders = None):
    """ Replace Sympy symbols in scalars in groups with the replacements subs_rules,
    keeping split_orders updated if need be.
    """
    temp = [type(group[0])(sympify(ncprod.scalar).xreplace(subs_rules),
                     ncprod.product) for ncprod in group]
    remove_zeros(temp)
    if split_orders is not None:
        split_orders[-1] = len(temp)
    return type(group)(temp)

def zero_free_vars(group, fvars):
    return substitute_group(group, dict(zip(fvars, [0]*len(fvars))))

def split_free_vars(group, fvars):
    """Split up group into a list of groups, each with only one fvar not zero"""
    zero_list = [0]*(len(fvars)-1)
    return [substitute_group(group, dict(zip([fvar2 for fvar2 in fvars if fvar2 != fvar],zero_list))) for fvar in fvars]





# Group Utility functions: conversion Sigma <--> Majorana, printing, saving, loading...


def convert_group(group):
    if isinstance(group, TransInvSum):
        return group.convert()
    if not isinstance(group, NCSum):
        group = NCSum(group)
    if group.is_zero():
        return []
    if isinstance(group[0], SigmaProduct):
        return NCSum([convert_from_sigma(el) for el in group])
    elif isinstance(group[0], MajoranaProduct):
        return NCSum([convert_to_sigma(el) for el in group])
    else:
        raise ValueError('Unrecognised conversion asked for!')



def print_group(group, breaks = True, neglect_order = None):
    if not hasattr(print_group, 'orders'):
        print_group.orders = {}
    if isinstance(group, NCSum):
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
    if isinstance(group, NCSum):
        if len(group) > 1:
            if not newlines:
                group = order_group(group, print_group.orders)
                return('$$'+' + '.join(a.texify() for a in group).replace('+ -', '-').replace('\\\\','\\')+'$$')
            else:
                group = order_group(group, print_group.orders)
                return('\\begin{align*}\n'+'\\\\\n'.join('&' + a.texify().replace('\\\\','\\') for a in group)+'\n\\end{align*}')
        else:
            return('$$'+group[0].texify().replace('\\\\','\\')+'$$')
    elif isinstance(group, TransInvSum):
        return group.texify(print_group.orders)
    else:
        return('$$'+group.texify().replace('\\\\','\\')+'$$')

def save_group(group, filename,
               iofvars = None,
               split_orders = None,
               normdict = None,
               zeroth_order = None):
    ext = '.p'
    if filename[-2:] != ext:
            filename = filename + ext
    if iofvars is None:
        iofvars = []
    if split_orders is None:
        split_orders = []
    if normdict is None:
        normdict = []
    if zeroth_order is None:
        zeroth_order = []
    data = OrderedDict([('group', group),
                        ('iofvars', iofvars),
                        ('split_orders', split_orders),
                        ('normdict', normdict),
                        ('zeroth_order', zeroth_order)])
    with open(filename, "wb" ) as f:
        pickle.dump(data, f, protocol=pickle.HIGHEST_PROTOCOL)

def load_group(filename,
               iofvars = None,
               split_orders = None,
               normdict = None,
               zeroth_order = None):
    ext = '.p'
    if filename[-2:] != ext:
            filename = filename + ext
    with open(filename, "rb") as f:
        parsed = pickle.load(f)
    if iofvars is not None:
        iofvars[:] = parsed['iofvars']
    if split_orders is not None:
        split_orders[:] = parsed['split_orders']
    if normdict is not None:
        normdict.update(parsed['normdict'])
    if zeroth_order is not None:
        zeroth_order[:] = parsed['zeroth_order']
    return parsed['group']


# Subspace functions: find the subspace of operator-string "vector space"
# to invert over to find solutions X to the operator equation:

# [Jpart, X] + to_cancel = 0

# Also find the projections of the operators as vectors into this space, and the
# projections of superoperators, in particular the commutator, as matrices.


def sparse_fill_subspace_rows(to_cancel, matrixrows, subspace, Jpart, ind_col):
    """Helper function for sparse_find_subspace"""
    ind_cols = [ind_col]
    to_cancels = [to_cancel]
    count = 0
    while count < len(ind_cols):
        comm = calculate_commutator(Jpart, to_cancels[count].get_unit())
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

def sparse_find_subspace(to_cancel, Jpart, verbose = False):
    """ Find the appropriate vector subspace to solve [Jpart, X] + to_cancel = 0
    by recursively commutating Jpart with to_cancel. In the process, also
    build the rows of the matrix [J, . ] into matrixrows, in triplet form:
    the indices of matrixrows are the row numbers, and the values are
    tuple (column number, scalar value).
    """
    subspace = OrderedDict()
    matrixrows = {}
    if verbose:
        print("Finding subspace to invert over...")
    length = len(to_cancel)
    for i, ncprod in enumerate(to_cancel):
        if not tuple(ncprod.product) in subspace:
            ind = len(subspace)
            subspace[tuple(ncprod.product)] = ind
            matrixrows[ind] = []
            sparse_fill_subspace_rows(ncprod, matrixrows, subspace, Jpart, ind)
        if verbose:
            print_progress(i, length)
    return subspace, matrixrows

def build_vector(to_cancel, subspace):
    """ Project an operator to_cancel into the vector subspace."""
    cvector = [0]*len(subspace)
    for ncprod in to_cancel:
        cvector[subspace[tuple(ncprod.product)]] = ncprod.scalar
    return cvector

def print_subspace(subspace, typ = MajoranaProduct):
    for key, item in subspace.items():
        print(str(item)+ ': ' + ' '.join([typ.stringify(a) for a in key]))



##Experimental TransInvSum stuff PLEASE IGNORE
def sparse_fill_subspace_rows_double_comm(to_cancel, matrixrows, subspace, Jpart, ind_col):
        ind_cols = [ind_col]
        to_cancels = [to_cancel]
        count = 0
        while count < len(ind_cols):
            comm = calculate_commutator(Jpart, calculate_commutator(I*Jpart, TransInvSum(to_cancels[count].get_unit())))
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

##Experimental TransInvSum stuff PLEASE IGNORE
def sparse_find_subspace_double_comm(to_cancel, Jpart):
    subspace = OrderedDict()
    matrixrows = {}
    length = len(to_cancel)
    for i, ncprod in enumerate(to_cancel):
        if not tuple(ncprod.product) in subspace:
            ind = len(subspace)
            subspace[tuple(ncprod.product)] = ind
            matrixrows[ind] = []
            sparse_fill_subspace_rows_double_comm(ncprod, matrixrows, subspace, Jpart, ind)
        print_progress(i, length)
    return subspace, matrixrows



# Functions to call symbolic (or numeric) linear algebra solvers to actually find
# the solution to the equation we have set up. We split the subspace we have projected
# to into further non-interacting blocks called sub_subspaces, then solve their
# resulting matrix equation on each sub_subspace by calling an external solver.

def solve_for_sub_subspace(matrixrows, sub_sub_space, coeffs, cvector,
                           iofvars,fvargen, newfvars, vardict,
                           numeric_dict = None, homogeneous = False):
    """Helper function for sparse_linear_solve which calls the external solver,
    for a given sub_subspace. Do NOT call directly, use sparse_linear_solve.
    """
    sspacedict = dict(zip(sub_sub_space, range(len(sub_sub_space))))
    length = len(sub_sub_space)
    sparse_mat_rep = OrderedDict()
    sub_cvector = []
    rownumstore = [] #for debugging
    row_count = 0
    for rownum, row in matrixrows.items():
        if row and row[0][0] in sub_sub_space:
            rownumstore.append(rownum)
            sub_cvector.append(cvector[rownum])
            for el in row:
                sparse_mat_rep[(row_count, sspacedict[el[0]])] = el[1]
            row_count += 1
    fvars = [coeffs[ind] for ind in sub_sub_space]
    if homogeneous:
        solver = n_h_linear_solve if numeric_dict else h_linear_solve
    else:
        solver = n_linear_solve if numeric_dict else linear_solve
    argdict = numeric_dict if numeric_dict else vardict
    sols = solver(sparse_mat_rep, sub_cvector, length, fvars,
                        iofvars, fvargen, newfvars, argdict)
    if not sols:
        print(sparse_mat_rep)
        print(sub_cvector)
        print(fvars)
        print(rownumstore)
        print(iofvars)
        raise ValueError("Failure. No solutions.")
    return sols



def find_sub_subspaces(matrixrows):
    return [list(space) for space in merge([[el[0] for el in row] for rownum, row in matrixrows.items()])]


def sparse_linear_solve(cvector, matrixrows, subspace, vardict,
                        typ = MajoranaProduct,
                        fvarname = 'A',
                        iofvars = None,
                        numeric_dict = None,
                        homogeneous = False,
                        verbose = False):
    """ Finds the solution to matrixrows * X + cvector = 0. over a subspace,
    by solving over non-interacting sub_subspaces.

    fvarname is the prefix for any free variables generated if the solution is underdetermined. These
    variables will overwrite the list iofvars if it is not none. (Currently there
    is no support for input iofvars.) Symbolic unless numeric_dict is not None,
    which should be a list of replacement rules for symbols with floats.
    Set homogenous = True if cvector = 0 for numeric matrices.
    """

    fvar_gen = sympy.numbered_symbols('fvar')
    fvars = [next(fvar_gen) for i in range(len(subspace))]
    sub_sub_spaces = find_sub_subspaces(matrixrows)
    if verbose:
        print("Solving operator matrix equation...")
        print("Subsubspaces: ")
        print(sub_sub_spaces) #Quite verbose (should make it a level)
        print("Number: " + str(len(sub_sub_spaces)))
        print("Dimensions: " +str([len(ss) for ss in sub_sub_spaces]))
    solutions = {}
    length_ss = len(sub_sub_spaces)
    fvar_gen = sympy.numbered_symbols(fvarname)
    newfvars = []
    for i, ss_space in enumerate(sub_sub_spaces):
        solutions.update(solve_for_sub_subspace(matrixrows, ss_space,
                                                fvars, cvector, iofvars,
                                                fvar_gen, newfvars, vardict,
                                                numeric_dict = numeric_dict,
                                                homogeneous = homogeneous))
        print_progress(i, length_ss)
    solvector = []
    subs_rules = {}
    for fvar in fvars:
        try:
            solvector.append(solutions[fvar])
        except KeyError:
            newfvars.append(next(fvar_gen))
            solvector.append(newfvars[-1])
            subs_rules[fvar] = newfvars[-1]
    if newfvars:
        if iofvars is not None:
            iofvars[:] = newfvars
        solvector[:] = (vec.xreplace(subs_rules) for vec in solvector)
    return simplify_group([typ(-solvector[i], list(key))
                                      for i,key in enumerate(subspace)])



def solve_commutator_equation(Jpart,
                              to_cancel,
                              vardict,
                              fvarname = 'F',
                              iofvars = None,
                              numeric_dict = None,
                              homogeneous = False,
                              verbose = False,
                              delete_to_cancel = True):
    """ Solve the equation [Jpart, X] + to_cancel = 0 for X. See above for
    description of arguments.
    """
    typ = type(Jpart[0])
    subspace, matrixrows = sparse_find_subspace(to_cancel, Jpart, verbose = verbose)
    cvector = build_vector(to_cancel, subspace)
    if delete_to_cancel:
        del to_cancel
    if verbose:
        print_subspace(subspace, typ = typ)
    if iofvars is None:
        iofvars = []
    return sparse_linear_solve(cvector,
                               matrixrows,
                               subspace,
                               vardict,
                               typ = typ,
                               iofvars = iofvars,
                               fvarname = fvarname,
                               numeric_dict = numeric_dict,
                               homogeneous = homogeneous,
                               verbose = verbose)



def _clean_to_cancel_of_fvars(check_unsolvable, to_cancel, fvars, vardict, verbose):
    """ Removes free variables from expression to cancel. Takes a
    pragmatic approach: sets the necessary free variables so that the
    operator equation *can* be solved at all, then sets all others to 0.
    """
    cvector = []
    solutions = {}
    matrixrows = {}
    ind = 0
    if verbose:
        #print(to_cancel) # Too verbose even for verbose.
        print("Setting free variables in solution...")
    for ncprod in to_cancel:
        if check_unsolvable(ncprod):
            term = sympy.expand(ncprod.scalar)
            matrixrows[ind] = []
            row = matrixrows[ind]
            for ind_col, coeff in enumerate(fvars):
                product = term.coeff(coeff)
                if product != 0:
                    row.append((ind_col, product))
            const_term = term.as_coeff_add(*fvars)[0]
            cvector.append(const_term)
            ind+=1
    sub_sub_spaces = find_sub_subspaces(matrixrows)
    length_ss = len(sub_sub_spaces)
    for i, ss_space in enumerate(sub_sub_spaces):
        solutions.update(solve_for_sub_subspace(matrixrows, ss_space,
                                                fvars, cvector, None,
                                                None, None, vardict))
        print_progress(i, length_ss)
    for fvar in fvars:
        if fvar not in solutions:
            if verbose:
                print("Zeroing fvar: " +str(fvar))
            solutions[fvar] = 0
    return substitute_group(to_cancel, solutions)


def _solve_single_nofvar(single, to_cancel):
    result = type(to_cancel)([])
    for ncprod in to_cancel:
        comm = single[0].commute(ncprod)
        if comm.scalar == 0:
            raise ValueError("Not invertible: " + str(ncprod))
        result.extend([-comm*sympy.Rational(1, 4)])
    return result

def solve_single_term_commutator_equation(single, to_cancel, fvars, vardict, verbose=False):
    """ Solve the equation [single, X] + to_cancel = 0 for X, where single
    is a single operator string, so we can use the fact that [single, .] is
    idempotent if it does not vanish.
    """
    single = NCSum(single)
    if len(single) > 1:
        raise ValueError("Single term solver called but multiple terms in commutator.")
    simplify_remove_zeros(to_cancel)
    if fvars:
        singleunit = single[0].get_unit()
        def check_unsolvable(ncprod):
            return singleunit.commute(ncprod.get_unit()).scalar == 0
        to_cancel = _clean_to_cancel_of_fvars(check_unsolvable, to_cancel, fvars, vardict, verbose)
    return _solve_single_nofvar(single, to_cancel)


# Functions to invert [a1, Gs[-1]] = to_invert to find Gs[-1],
# where a1 is the first Majorana = sz1.
# A significantly easier problem than the general case above.

def _check_a1_unsolvable(ncprod):
    """See if invertible by checking if length even w/out 1 if present"""
    first = (ncprod.product[0] == 1)
    return (len(ncprod)-first*1) % 2 == 0

def _solve_a1_G_nofvar(to_cancel, Gs):
    for ncprod in to_cancel:
        ncprod.scalar = sympy.simplify(ncprod.scalar) #unfortunately this seems necessary atm.
        if ncprod.scalar != 0:
            if _check_a1_unsolvable(ncprod):
                raise ValueError("Not invertible: " + str(ncprod))
            else:
                if ncprod.product[0] == 1:
                    Gs[-1] += MajoranaProduct(-I/2*ncprod.scalar, ncprod.product[1:])
                else:
                    Gs[-1] += MajoranaProduct(-I/2*ncprod.scalar, [1]+ncprod.product[:])

def solve_a1_G(to_cancel, Gs, fvars, vardict, verbose = False):
    Gs.append([])
    if fvars:
        to_cancel = _clean_to_cancel_of_fvars(_check_a1_unsolvable, to_cancel, fvars, vardict, verbose)
    _solve_a1_G_nofvar(to_cancel, Gs)


# Functions to build the entire vector subspace or the entire subspace
# on some finite support/parity of length of operator strings.
# Useful for comparison to recursive method (which is much more efficient).

def build_entire_subspace(L, typ = MajoranaProduct):
    return [typ(1, sorted(list(el))) for el in powerset(range(1, L*2+1))]

def build_norm_subspace(L, typ = MajoranaProduct):
    return [typ(1, sorted(list(el))) for el in powerset(range(2, L*2+1))]

def build_odd_subspace(L, typ = MajoranaProduct):
    return [ncprod for ncprod in build_entire_subspace(L, typ = typ) if len(ncprod.product) % 2 == 1]

def build_odd_norm_subspace(L, typ = MajoranaProduct):
    return [ncprod for ncprod in build_norm_subspace(L, typ = typ) if len(ncprod.product) % 2 == 1]

def build_finite_support_subspace(L, L_supp, typ = MajoranaProduct, periodic = False):
    return [ncprod for ncprod in build_entire_subspace(L, typ = typ)
            if get_support(ncprod, L_periodic = L if periodic else None) <= L_supp]

#Find the maximum possible subspace for Jpart at order order using the recursive algorithm.
def fill_subspace(Jpart, order):
    L = 2*order+1
    typ = type(Jpart[0])
    subspace_ops = build_odd_norm_subspace(L)
    len_subspace = len(subspace_ops)
    matrixrows = defaultdict(list)
    subspace = OrderedDict(zip([tuple(el.product) for el in subspace_ops],
                                [i for i in range(len_subspace)]))
    for product, ind_col in subspace.items():
        comm = calculate_commutator(Jpart, typ(1, list(product)))
        for ncprod in comm:
           ind_row = subspace[tuple(ncprod.product)]
           matrixrows[ind_row].append((ind_col, ncprod.scalar))
    return subspace, matrixrows

def fill_subspace_norm(Jpart, to_cancel, order):
    L = 2*order+1
    subspace_ops = build_odd_norm_subspace(L)
    return sparse_find_subspace(to_cancel+subspace_ops, Jpart)


# Try to find all operators which commutes with H in the maximally brute force way
# that is, solve [H, X] = 0 over the entire space. Obviously extremely inefficient!

def solve_at_once(H, L, iofvars = None, entire = False, numeric_dict = None):
    typ = type(H[0])
    if entire:
        subspace_ops = build_entire_subspace(L, typ = typ)
    else:
        subspace_ops = build_odd_subspace(L, typ = typ)
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
    return sparse_linear_solve(cvector,
                               matrixrows,
                               subspace,
                               homogeneous = True,
                               fvarname = 'F',
                               iofvars=iofvars,
                               typ= typ,
                               numeric_dict = numeric_dict)






# Functions for saving the groups to C functions.
def get_yaml_poly(poly, gens):
    s = sympy.Poly(poly, gens)
    pows = [[float(e) for e in el] for el in s.monoms()]
    coeffs = [float(sympy.N(coeff)) for coeff in s.coeffs()]
    return {'pows':pows,'coeffs':coeffs}

def get_c_function(scalar, couplings):
    args_replace = dict(zip([str(el) for el in couplings],
                            ['args['+str(i)+']' for i in range(len(couplings))]))
    return multiple_replace(ccode(scalar), args_replace)


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


def save_to_c_general(group,
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
    with open(filename+'.h', 'w') as header:
        header.write('mpreal norm(std::vector<mpreal> args){\n'
                     'return '+get_c_function(norm, couplings)+';\n}\n'
                     'std::vector<std::function<mpreal(std::vector<mpreal>)> > psi_coeffs = {')
        data = OrderedDict([('Psi',[])])
        for i, ncprod in enumerate(reversed(group)):
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
                              })
            header.write('[](std::vector<mpreal> args){return '
                         +get_c_function(sympy.simplify(sigprod.scalar), couplings)
                         +(';},' if i < (len(group)-1) else ';}};'))
    write_yaml(data, filename)
