from sympy import I, symbols
from mathematica_printer import mstr
from commutator import load_group, substitute_group, calculate_commutator, square_to_find_identity_scalar
from commutator import Ncproduct

N = Ncproduct
NAME = "j2_small_gs_r11.yaml"
ORDER = 0
END_ORDER = 11

L=40
V,J,f, J2 = symbols('V J f J2')
orders = {V:1,f:1, J2:1}

fpart = [N(I*f, [2*j+1,2*j+2]) for j in range(L)]
Vpart = [N(V, [2*j+1,2*j+2,2*j+3,2*j+4]) for j in range(L-2)]
V2part = [Ncproduct(J2, [2*j+2,2*j+3,2*j+4,2*j+5]) for j in range(L-2)]
small = fpart+V2part

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
