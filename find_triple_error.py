from sympy import I, symbols, simplify
from mathematica_printer import mstr
from commutator import load_group, substitute_group, calculate_commutator, square_to_find_identity
from commutator import Ncproduct, trace_inner_product

N = Ncproduct
NAME = "jto1gpsi_r6.yaml"
ORDER = 0
END_ORDER = 6

L=40
V,J,f = symbols('V J f')
orders = {V:1,f:1}

fpart = [N(I*f, [2*j+1,2*j+2]) for j in range(L)]
Vpart = [N(V, [2*j+1,2*j+2,2*j+3,2*j+4]) for j in range(L-2)]
Jpart = [N(I, [2*j+2,2*j+3]) for j in range(L-1)]
small = fpart+Vpart
H = small + Jpart

iofvars = []
split_orders = []
normdict = {}
psi = load_group(NAME, iofvars = iofvars, split_orders = split_orders,
                 normdict = normdict)
psi = substitute_group(psi, normdict)
for order in range(ORDER, END_ORDER+1):
    error = calculate_commutator(small, psi[split_orders[order]:split_orders[order+1]])
    error =  calculate_commutator(H, error)
    error =  calculate_commutator(H, error)
    #print(psi[split_orders[order]:split_orders[order+1]])
    #print(error)
    error_norm = trace_inner_product(error, psi[0:split_orders[order+1]])
    print(mstr(simplify(error_norm)))
    #print(order)
