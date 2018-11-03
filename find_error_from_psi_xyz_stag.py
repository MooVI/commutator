from sympy import I, symbols
from mathematica_printer import mstr
from commutator import load_group, substitute_group, calculate_commutator, square_to_find_identity_scalar
from commutator import Ncproduct

N = Ncproduct
NAME = "xyz_stagDsym_psi_r8"
ORDER = 0
END_ORDER = 8

X, D = symbols('X D', real=True)

L = 40
orders = {X:1,D:1}

Ypart = [N(I*((X-D) if j%2==0 else (X+D)), [2*j+1,2*j+4]) for j in range(L-2)]
Xpart = [N(-((X-D) if j%2==0 else (X+D)), [2*j+1,2*j+2,2*j+3,2*j+4]) for j in range(L-2)]
small = Xpart+Ypart

iofvars = []
split_orders = []
normdict = {}
psi = load_group(NAME, iofvars = iofvars, split_orders = split_orders,
                 normdict = normdict)
psi = substitute_group(psi, normdict)
for order in range(ORDER, END_ORDER+1):
    error = calculate_commutator(small, psi[split_orders[order]:split_orders[order+1]])
    #print(psi[split_orders[order]:split_orders[order+1]])
    #print(error)
    error_norm = square_to_find_identity_scalar(error)
    print(mstr(error_norm))
    #print(order)
