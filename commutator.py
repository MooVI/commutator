import tempfile
import os
from collections import OrderedDict, defaultdict
from itertools import chain, combinations
import re
from sympy import  I, sympify, S, factorial
from sympy.printing.ccode import ccode
from bisect import bisect_right
import yaml
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
    Returns true unless wrapped around so result is lexicographically smaller. '''
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
    def __init__(self, scalar, product):
        self.scalar = scalar
        if isinstance(product, list):
            self.product = product
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
        typ = type(self)
        if isinstance(other, typ):
            return typ(self.scalar*other.scalar, self.product+other.product)
        else:
            return typ(self.scalar*other, self.product)

    def __rmul__(self, other):
        typ = type(self)
        if isinstance(other, typ):
            return typ(self.scalar*other.scalar, other.product+self.product)
        else:
            return typ(self.scalar*other, self.product)

    def __repr__(self):
        return str(self.scalar) + ' : ' +str(self.product)

    def __eq__(self, other):
        """Note == ignores the scalar field!"""
        try:
            return self.product == other.product
        except:
            return False

    def __str__(self):
        strings, scalar = self.self_stringify()
        return '\u202f'.join([str(self.scalar*scalar).replace('**','^').replace('*','\u202f').replace('I','i')]
                             + strings)


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
        for a pre-SORTED, NON-DUPLICATE product."""
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

def count_ys_in_paul_list(a):
    l = len(a) - 1
    n = 0
    i = 0
    while i < l:
        n += (a[i+1] == a[i] + 1)
        i += 2
    return n

class SigmaProduct(NCProduct):

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
        for a pre-SORTED, NON-DUPLICATE product."""
        return 1-2*(count_ys_in_paul_list(product)%2)

    def reverse_sign(self):
        return 1-2*(count_ys_in_paul_list(self.product)%2)

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




# Functions for dealing with 'groups' = lists of NCProducts
# representing **sums**.

# After any operation which might zero elements of a sum,
# or add duplicates which need to be merged simplify_group should be called.

# full_simplify_group calls sympy.simplfy on each element. This is **slow**
# and should only be done when necesssary. mathematica_simplify_group does the
# same for mathematica.

# All of the group utility functions (multiply_groups, add_groups, etc.)
# do this automatically

# TODO: Make this a class with proper syntatic sugar.


# Simplification of groups

def collect_terms(group):
    from collections import defaultdict
    D = defaultdict(list)
    for i,ncprod in enumerate(group):
        D[tuple(ncprod.product)].append(i)
    return type(group)(type(group[0])(
        sympy.expand(sum([group[i].scalar for i in D[key]])), list(key))
                       for key in D)

def full_collect_terms(group):
    from collections import defaultdict
    D = defaultdict(list)
    for i,ncprod in enumerate(group):
        D[tuple(ncprod.product)].append(i)
    return type(group)(type(group[0])(
        sympy.simplify(sum([group[i].scalar for i in D[key]])), list(key))
                       for key in D)

def math_collect_terms(group):
    from collections import defaultdict
    D = defaultdict(list)
    for i,ncprod in enumerate(group):
        D[tuple(ncprod.product)].append(i)
    return type(group)(type(group[0])(
        mathematica_simplify(sum([group[i].scalar for i in D[key]])), list(key))
                       for key in D)

def remove_zeros(group):
    group[:] = (a for a in group if a.scalar != 0)

def simplify_group(group):
    if isinstance(group, TransInvSum):
        return group.simplify()
    remove_zeros(group)
    for ncprod in group:
        ncprod.sort()
        set_squares_to_identity(ncprod)
    group = collect_terms(group)
    remove_zeros(group)
    return group

def full_simplify_group(group):
    if isinstance(group, TransInvSum):
        print("Full simplify not supported for TransInv.")
        return group.simplify()
    remove_zeros(group)
    for ncprod in group:
        ncprod.sort()
        set_squares_to_identity(ncprod)
    group = full_collect_terms(group)
    remove_zeros(group)
    return group

def mathematica_simplify_group(group, use_tempfile = False):
    if isinstance(group, TransInvSum):
        print("Mathematica simplify not supported for TransInv.")
        return group.simplify()
    remove_zeros(group)
    for ncprod in group:
        ncprod.sort()
        set_squares_to_identity(ncprod)
    group = math_collect_terms(group, use_tempfile)
    remove_zeros(group)
    return group



# Basic Group Mathematical Functions--------------------

def conjugate_group(group):
    return [a.conjugate() for a in group]

def hermitian_conjugate_group(group):
    return [a.hermitian_conjugate() for a in group]

def postmultiply(group,a):
     return type(group)([b*a for b in group])

def premultiply(a, group):
    return type(group)([a*b for b in group])

def multiply_groups(group_a, group_b):
    return simplify_group([a*b for a in group_a for b in group_b])

def add_groups(group_a, group_b):
    return simplify_group(group_a+group_b)

def subtract_groups(group_a, group_b):
    return simplify_group(group_a+premultiply(-1, group_b))

def commute_group(group_a, group_b):
    """For known groups. Prefer using calculate_commutator instead if
    could be individual NCProducts as well."""
    result = []
    for a in group_a:
        for b in group_b:
            if not (a.is_identity() or b.is_identity()):
                result.append(a.commute(b))
    return result

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
        if not isinstance(group_a, list):
                group_a = [group_a]
        if not isinstance(group_b, list):
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
    return sorted(group, key = lambda a: find_order(a.scalar,orders))

def check_group_at_least_order(group, order, orders):
    """Check that group only has elements of order order or higher.
    Use to check that perturbation theory is being satisfied."""
    for ncprod in group:
        if find_order(ncprod.scalar, orders) <= order:
            test = sympy.simplify(ncprod.scalar)
            if test != 0 and find_order(test, orders) <= order:
                print('Error: ' + str(test) + ': ' + str(ncprod.product))
                return False
    return True

def commute_up_to_order(group_a, group_b, order, split_orders_a, split_orders_b):
    result = []
    a_order_max = len(split_orders_a)-1
    b_order_max = len(split_orders_b)-1
    if a_order_max + b_order_max <= order:
        return commute_group(group_a, group_b)
    for aorder in range(a_order_max):
        for border in range(min([order-aorder, b_order_max])):
            result += commute_group(group_a[split_orders_a[aorder]:split_orders_a[aorder+1]],
                                    group_b[split_orders_b[border]:split_orders_b[border+1]])
    return result

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
    group_type = type(to_trans)
    result = group_type(to_trans)
    for torder in range(1,max_order+1):
        for orderset in accel_asc(torder):
            numcomms = len(orderset)
            if inverse:
                taylorscalar = ((-I)**numcomms)/factorial(numcomms)
            else:
                taylorscalar = (I**numcomms)/factorial(numcomms)
            while True:
                cumul = group_type(to_trans)
                for order in orderset:
                    cumul = calculate_commutator(cumul, Gs[order-1])
                result += premultiply(taylorscalar,cumul)
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
    group_type = type(to_trans)
    if torder == 0:
        return group_type(to_trans)
    result =  group_type([])
    for orderset in accel_asc(torder):
        numcomms = len(orderset)
        if numcomms != 1 or not not_single_comm:
            if inverse:
                taylorscalar = ((-I)**numcomms)/factorial(numcomms)
            else:
                taylorscalar = (I**numcomms)/factorial(numcomms)
            while True:
                cumul = group_type(to_trans)
                for order in orderset:
                    cumul = calculate_commutator(cumul, Gs[order-1])
                result += premultiply(taylorscalar,cumul)
                if not next_permutationS(orderset):
                    break
    return simplify_group(result)

def exponentiate_to_order(Gs, torder, inverse = False):
    """Calculate U = e^i(Gs[0]+Gs[1]+Gs[2]+...), see above."""
    if torder == 0:
        return [type(Gs[0][0])(1,[])]
    result = []
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
            result += premultiply(taylorscalar,cumul)
            if not next_permutationS(orderset):
                break
    return simplify_group(result)



# Functions for dealing with free variables in groups.

def substitute_group(group, subs_rules, split_orders = None):
    """ Replace Sympy symbols in scalars in groups with the replacements subs_rules,
    keeping split_orders updated if need be. """
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
    if not isinstance(group, list):
        group = [group]
    if isinstance(group[0], SigmaProduct):
        return [convert_from_sigma(el) for el in group]
    elif isinstance(group[0], MajoranaProduct):
        return [convert_to_sigma(el) for el in group]
    else:
        raise ValueError('Unrecognised conversion asked for!')



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
    elif isinstance(group, TransInvSum):
        return group.texify(print_group.orders)
    else:
        return('$$'+group.texify().replace('\\\\','\\')+'$$')

def save_group(group, filename, iofvars=None, split_orders = None, normdict = None):
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

def sparse_find_subspace(to_cancel, Jpart):
    """ Find the appropriate vector subspace to solve [Jpart, X] + to_cancel = 0
    by recursively commutating Jpart with to_cancel. In the process, also
    build the rows of the matrix [J, . ] into matrixrows, in triplet form:
    the indices of matrixrows are the row numbers, and the values are
    tuple (column number, scalar value).
    """
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
                           iofvars,fvargen, newfvars,
                           numeric_dict = None, homogeneous = False):
    """Helper function for sparse_linear_solve which calls the external solver,
    for a give sub_subspace. Do NOT call directly, use sparse_linear_solve."""
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
    argdict = numeric_dict if numeric_dict else mathematica_parser.vardict
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


def sparse_linear_solve(cvector, matrixrows, subspace,
                        typ = MajoranaProduct,
                        group_type = list,
                        fvarname = 'A',
                        iofvars = None,
                        numeric_dict = None,
                        homogeneous = False):
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
    print("Subsubspaces: ")
    print(sub_sub_spaces)
    print("Number: " + str(len(sub_sub_spaces)))
    print("Lengths: " +str([len(ss) for ss in sub_sub_spaces]))
    solutions = {}
    length_ss = len(sub_sub_spaces)
    fvargen = sympy.numbered_symbols(fvarname)
    newfvars = []
    for i, ss_space in enumerate(sub_sub_spaces):
        solutions.update(solve_for_sub_subspace(matrixrows, ss_space,
                                                fvars, cvector, iofvars,
                                                fvargen, newfvars,
                                                numeric_dict = numeric_dict,
                                                homogeneous = homogeneous))
        print_progress(i, length_ss)
    solvector = []
    for fvar in fvars:
        try:
            solvector.append(solutions[fvar])
        except KeyError:
            newfvars.append(next(fvargen))
            solvector.append(newfvars[-1])
    if newfvars:
        if iofvars is not None:
            iofvars[:] = newfvars
    return simplify_group(group_type([typ(-solvector[i], list(key))
                                      for i,key in enumerate(subspace)]))



# Functions to invert [a1, Gs[-1]] = to_invert to find Gs[-1] --
# a significantly easier problem than the general case above.

def _check_sz1_uninvertible(ncprod):
    """See if invertible by checking if length even w/out 1 if present"""
    first = (ncprod.product[0] == 1)
    return (len(ncprod)-first*1) % 2 == 0

def _invert_sz1_G_nofvar(to_invert, Gs):
    for ncprod in to_invert:
        if sympy.simplify(ncprod.scalar) != 0:
            if _check_sz1_uninvertible(ncprod):
                raise ValueError("Not invertible: " + str(ncprod))
            else:
                if ncprod.product[0] == 1:
                    Gs[-1] += MajoranaProduct(-I/2*ncprod.scalar, ncprod.product[1:])
                else:
                    Gs[-1] += MajoranaProduct(-I/2*ncprod.scalar, [1]+ncprod.product[:])


def _clean_to_invert_of_fvars(to_invert, fvars):
    cvector = []
    solutions = {}
    matrixrows = {}
    ind = 0
    print(to_invert)
    for ncprod in to_invert:
        if _check_sz1_uninvertible(ncprod):
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
                                                None, None))
        print_progress(i, length_ss)
    for fvar in fvars:
        if fvar not in solutions:
            print("Zeroing fvar: " +str(fvar))
            solutions[fvar] = 0
    return substitute_group(to_invert, solutions)


def invert_sz1_G(to_invert, Gs, fvars):
    Gs.append([])
    if fvars:
        to_invert = _clean_to_invert_of_fvars(to_invert, fvars)
    _invert_sz1_G_nofvar(to_invert, Gs)


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
    subspace_ops = build_odd_norm_subspace(L)
    len_subspace = len(subspace_ops)
    matrixrows = defaultdict(list)
    subspace = OrderedDict(zip([tuple(el.product) for el in subspace_ops],
                                [i for i in range(len_subspace)]))
    for product, ind_col in subspace.items():
        comm = calculate_commutator(Jpart, MajoranaProduct(1, list(product)))
        for ncprod in comm:
           ind_row = subspace[tuple(ncprod.product)]
           matrixrows[ind_row].append((ind_col, ncprod.scalar))
    return subspace, matrixrows

def fill_subspace_norm(Jpart, to_cancel, order):
    L = 2*order+1
    subspace_ops = build_odd_norm_subspace(L)
    return sparse_find_subspace(to_cancel+subspace_ops, Jpart)


# Try to find all operators which commutes with H in the maximally brute force way
# that is, solve [H,X] = 0 over the entire space. Obviously extremely inefficient!

def solve_at_once(H, L, iofvars = None, entire = False, numeric_dict = None):
    group_type = type(H[0])
    if entire:
        subspace_ops = build_entire_subspace(L, typ = group_type)
    else:
        subspace_ops = build_odd_subspace(L, typ = group_type)
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
                               typ= group_type,
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
