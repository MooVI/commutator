import commutator as comm
from sympy import I, symbols, sympify, sqrt, Add, factorial, S
from collections import defaultdict
from itertools import permutations

p = comm.print_group
c = comm.calculate_commutator
N = comm.Ncproduct


PSI_NAME = ""
NORM_NAME = ""
OUTPUT = ""
V, f = symbols('V f', real=True)
vardict = {'f': f, 'V':V}
orders = {V:1,f:1}

p.orders = orders
comm.mathematica_parser.vardict = vardict

START_ORDER = 1
END_ORDER = 5



def accel_asc(n):
    """From http://jeromekelleher.net/generating-integer-partitions.html
    Generates integer partions"""
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


def get_G_contributions(Gs, torder):
    result = 0
    for orderset in accel_asc(torder):
        numcomms = len(orderset)
        if numcomms > 1:
            taylorscalar = -(I**numcomms)/factorial(numcomms)
            #minus as we take from lhs to rhs
            for orderperm in permutations(orderset):
                cumul = N(1, [1])
                for order in orderperm:
                    cumul = c(cumul, Gs[order-1])
                result += comm.premultiply(taylorscalar,cumul)
    return comm.simplify_group(result)



with open(NORM_NAME, 'r') as norm_file:
    norm = sympify(norm_file.readline().strip(), vardict)
    taylorexpand = comm.mathematica_series(1/sqrt(norm), [V,f], END_ORDER)
    collectnorms = defaultdict(list)
    for expr in taylorexpand.args:
        collectnorms[comm.find_order(expr,orders)] = expr
    norms = {order: Add(*exprs) for order, exprs in collectnorms}

normdict = {}
iofvars = []
split_orders = []
psi = comm.load_group(PSI_NAME, normdict=normdict, iofvars=iofvars,
                 split_orders=split_orders)
psi = comm.substitute_group(psi, normdict)
Gs = [] #U = exp(i (G[0]+G[1]+G[2]+...)


for torder in range(START_ORDER, END_ORDER+1):
    to_cancel = []
    Gs[torder-1] = []
    for psiorder in range(torder+1):
        to_cancel += comm.premultiply(
            norms[torder-psiorder],
            psi[split_orders[psiorder]:split_orders[psiorder+1]])
    to_cancel = comm.simplify_group(to_cancel
                                    + get_G_contributions(Gs, torder))
    for ncprod in to_cancel:
        if comm.simplify(ncprod.scalar) != 0:
            first = (ncprod.product[0] == 1)
            #See if inveritble by checking if length even w/out 1 if present
            if (len(ncprod)-first*1) % 2 == 0:
                raise ValueError("Not invertible: " + str(ncprod))
            else:
                if first:
                    Gs[torder-1] += N(1/S(2)*ncprod.scalar, ncprod.product[1:])
                else:
                    Gs[torder-1] += N(1/S(2)*ncprod.scalar, [1]+ncprod.product[:])
    comm.save_group(Gs, OUTPUT+'_r'+str(torder))
