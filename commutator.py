from sympy import symbols, I, pretty, sympify
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

    def destringify(self, string):
        result = []
        for op in string.split(' '):
            if op[0] == 'a':
                result.append(int(op[1:])*2-1)
            elif op[0] == 'b':
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

def simplify_group(group):
    from collections import defaultdict
    D = defaultdict(list)
    for i,ncprod in enumerate(group):
        D[tuple(ncprod.product)].append(i)
    return [Ncproduct(sum([group[i].scalar for i in D[key]]), list(key)) for key in D]


def calculate_commutator(group_a,group_b):
    if not isinstance(group_a, list):
        group_a = [group_a]
    if not isinstance(group_b, list):
        group_b = [group_b]
    group = commute_group(group_a, group_b)
    remove_zeros(group)
    for ncprod in group:
        sort_anticommuting_product(ncprod)
        set_squares_to_identity(ncprod)
    group = simplify_group(group)
    remove_zeros(group)
    return group


def find_order(expr,orders):
    """Where order is power of small quantity, and orders a dict of
    symbols with their order."""
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

        
V,J,f = symbols('V J f')
L = 20
orders = {V:1,f:1,J:-1}
print_group.orders = orders
fpart = [Ncproduct(I*f, [2*j+1,2*j+2]) for j in range(L)]
Vpart = [Ncproduct(V, [2*j+1,2*j+2,2*j+3,2*j+4]) for j in range(L-2)]
small = fpart+Vpart
Jpart = [Ncproduct(I*J, [2*j+2,2*j+3]) for j in range(L-1)]
H = fpart + Vpart +Jpart
