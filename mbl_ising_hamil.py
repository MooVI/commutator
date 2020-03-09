from sympy import I, symbols, numbered_symbols
import commutator as comm

p = comm.print_group
N = comm.SigmaProduct


L = 20
V, J, J2, f = symbols('V J J2 f', real=True)
h_gen = numbered_symbols('h', real=True)
hs = [next(h_gen) for i in range(L+1)]
vardict = {'J2':J2, 'f': f, 'V':V, 'J':J}
vardict.update(zip([str(h) for h in hs],hs))
orders = {f:1}
Js = [J, J2]
p.orders = orders

fpart = [N(-f, ("x", [j])) for j in range(1, L+1)]
hpart = [N(-hs[j], ("z", [j])) for j in range(1, L+1)]
Jpart = [N(-J, ("zz", [j, j+1])) for j in range(1, L)]

small = fpart
large = Jpart + hpart
H = large + small

zeroth_order = [N(1, "z10")]

print("Hamiltonian:")
p(H, breaks = False)
