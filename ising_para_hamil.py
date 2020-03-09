from sympy import I, symbols
import commutator as comm

p = comm.print_group
N = comm.SigmaProduct

V, J, J2, f = symbols('V J J2 f', real=True)
L = 20
vardict = {'J2':J2, 'f': f, 'V':V, 'J':J}
orders = {J:1}
Js = [J, J2]
p.orders = orders

fpart = [N(-f, ("x", [j])) for j in range(1, L+1)]
Jpart = [N(-J, ("zz", [j, j+1])) for j in range(1, L)]

small = Jpart
large = fpart
H = large + small

zeroth_order = [N(1, "x10")]

print("Hamiltonian:")
p(H, breaks = False)
