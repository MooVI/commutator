from commutator import *
p = print_group
c = calculate_commutator
N = Ncproduct
fs = full_simplify_group
ms = mathematica_simplify_group
V,J,f = symbols('V J f')
L = 4
orders = {V:1,f:1}
print_group.orders = orders
fpart = [Ncproduct(I*f, [2*j+1,2*j+2]) for j in range(L-1)]
Vpart = [Ncproduct(V, [2*j+1,2*j+2,2*j+3,2*j+4]) for j in range(L-2)]
small = fpart+Vpart
Jpart = [Ncproduct(I*J, [2*j+2,2*j+3]) for j in range(L-1)]
H = fpart + Vpart +Jpart
iofvars = []
psi = solve_at_once(H,L, iofvars)
save_group(psi, 'solve_psi_' +str(L), iofvars=iofvars)
