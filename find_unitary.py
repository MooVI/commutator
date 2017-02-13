import commutator as comm
from sympy import I, symbols, sympify, sqrt, Add, factorial, S
from collections import defaultdict
from itertools import permutations
import sympy

p = comm.print_group
c = comm.calculate_commutator
N = comm.Ncproduct


#PSI_NAME = "jto1psi_matlab_r8.yaml"
#NORM_NAME = "jto1psi_matlab_norm"
#OUTPUT = "jto1uni"
#V, f = symbols('V f', real=True)
#vardict = {'f': f, 'V':V}
#orders = {V:1,f:1}

PSI_NAME = "psi_j2_matlab_large_r10"
NORM_NAME = "psi_j2_matlab_large_norm"
OUTPUT = "j2_smalluni"
J2, f = symbols('J2 f', real=True)
vardict = {'f': f, 'J2':J2}
orders = {J2:1,f:1}


# PSI_NAME = "xyzfast_r7.yaml"
# NORM_NAME = "xyzfast_norm"
# OUTPUT = "xyzuni"
# X, Y = symbols('X Y', real=True)
# ivars = [X,Y]
# vardict = {'X': X, 'Y':Y}
# orders = {X:1,Y:1}


p.orders = orders
comm.mathematica_parser.vardict = vardict
ivars = list(vardict.values())

START_ORDER = 1
END_ORDER = 10

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
    result = []
    for orderset in accel_asc(torder):
        numcomms = len(orderset)
        if numcomms > 1:
            taylorscalar = -(I**numcomms)/factorial(numcomms)
            #minus as we take from lhs to rhs
            while True:
                cumul = N(1, [1])
                for order in orderset:
                    cumul = c(cumul, Gs[order-1])
                result += comm.premultiply(taylorscalar,cumul)
                if not next_permutationS(orderset):
                    break
    return comm.simplify_group(result)


with open(NORM_NAME, 'r') as norm_file:
    norm = sympify(norm_file.readline().strip(), vardict)
    taylorexpand = comm.mathematica_series(1/sqrt(norm), ivars, END_ORDER)
    collectnorms = defaultdict(list)
    for expr in taylorexpand.args:
        collectnorms[comm.find_order(expr,orders)].append(expr)
    norms = {order: Add(*exprs) for order, exprs in collectnorms.items()}

normdict = {}
iofvars = []
split_orders = []
psi = comm.load_group(PSI_NAME, normdict=normdict, iofvars=iofvars,
                 split_orders=split_orders)
psi = comm.substitute_group(psi, normdict)
Gs = [] #U = exp(i (G[0]+G[1]+G[2]+...)
#We are solving U^\dagger \sigma^z_1 U = \Psi

for torder in range(START_ORDER, END_ORDER+1):
    to_cancel = []
    Gs.append([])
    for norder, norm in norms.items():
        if norder <= torder:
            psiorder = torder-norder
            to_cancel += comm.premultiply(
                norm,
                psi[split_orders[psiorder]:split_orders[psiorder+1]])
    to_cancel = comm.simplify_group(to_cancel
                                    + get_G_contributions(Gs, torder))
    for ncprod in to_cancel:
        if sympy.simplify(ncprod.scalar) != 0:
            first = (ncprod.product[0] == 1)
            #See if invertible by checking if length even w/out 1 if present
            if (len(ncprod)-first*1) % 2 == 0:
                raise ValueError("Not invertible: " + str(ncprod)
                                 + " at order " + str(torder))
            else:
                if first:
                    Gs[-1] += N(-I/2*ncprod.scalar, ncprod.product[1:])
                else:
                    Gs[-1] += N(-I/2*ncprod.scalar, [1]+ncprod.product[:])
    #comm.print_group(Gs[-1])
    comm.save_group(Gs, OUTPUT+'_r'+str(torder))
