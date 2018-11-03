from sympy import I, symbols
from mathematica_printer import mstr
from commutator import load_group, substitute_group, calculate_commutator
from commutator import Ncproduct, square_to_find_identity

N = Ncproduct
NAME = "jto1psi_matlab_r4.yaml"
ORDER = 1
END_ORDER = 3

L=25
V,J,f = symbols('V J f')
orders = {V:1,f:1}

fpart = [N(I*f, [2*j+1,2*j+2]) for j in range(L)]
Vpart = [N(V, [2*j+1,2*j+2,2*j+3,2*j+4]) for j in range(L-1)]
small = fpart+Vpart

iofvars = []
split_orders = []
normdict = {}
psi = load_group(NAME, iofvars = iofvars, split_orders = split_orders,
                 normdict = normdict)
psi = substitute_group(psi, normdict)
old_error = [-1*a for a in calculate_commutator(small, psi[split_orders[ORDER-1]:split_orders[ORDER]])]
for order in range(ORDER, END_ORDER+1):
    error = calculate_commutator(small, psi[split_orders[order]:split_orders[order+1]])
    #print(psi[split_orders[order]:split_orders[order+1]])
    #print(error)
    error_norm = square_to_find_identity(error+old_error)[0].scalar
    old_error = [-1*a for a in error]
    print(mstr(error_norm))
    #print(order)
