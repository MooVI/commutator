
from sympy import I
from mathematica_printer import mstr
from commutator import load_group, substitute_group, calculate_commutator, square_to_find_identity
from commutator import Ncproduct

N = Ncproduct
NAME = "psi_j2_matlab_r2"
ORDER = 1
END_ORDER = 2

L=40
largeH =  [N(I, [2*j+2,2*j+3]) for j in range(L-1)]
iofvars = []
split_orders = []
normdict = {}
psi = load_group(NAME, iofvars = iofvars, split_orders = split_orders,
                 normdict = normdict)
psi = substitute_group(psi, normdict)
for order in range(ORDER, END_ORDER+1):
    error = calculate_commutator(largeH, psi[split_orders[order]:split_orders[order+1]])
    print(psi[split_orders[order]:split_orders[order+1]])
    print(error)
    error_norm = square_to_find_identity(error)[0].scalar
    print(mstr(error_norm))
    print(order)
