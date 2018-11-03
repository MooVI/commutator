import commutator as comm
from sympy import I, symbols, latex

p = comm.print_group
c = comm.calculate_commutator
N = comm.Ncproduct
Ncproduct = N

V, J, J2, f = symbols('V J J2 f', real=True)
L = 20
comm.mathematica_parser.vardict = {'J2':J2, 'f': f, 'V':V, 'J':J}
orders = {V:1,f:1}
Js = [J, J2]
p.orders = orders
fpart = [Ncproduct(I*f, [2*j+1,2*j+2]) for j in range(L)]
Vpart = [Ncproduct(V, [2*j+1,2*j+2,2*j+3,2*j+4]) for j in range(L-2)]
small =  Vpart+fpart
J1part = [Ncproduct(I*Js[j%2], [2*j+2,2*j+5]) for j in range(L-1)]
J2part = [Ncproduct(I*J2, [2*j+2,2*j+4]) for j in range(L-2)]
Jpart =  J1part
H = Jpart+Vpart+fpart

START_ORDER = 1
END_ORDER = 8

#START_PSI = [N(1, 'a1')]

START_GS = []

START_IOFVARS = []
START_SPLIT_ORDERS = [0,1]
START_NORMDICT = {}

#START_PSI = comm.load_group('psi_j2_r8', START_IOFVARS, START_SPLIT_ORDERS, START_NORMDICT)

orders.update(zip(START_IOFVARS,[START_ORDER-1]*len(START_IOFVARS)))

FILEHEAD = 'Chris_f_gs'
NORM_AS_YOU_GO = True

# START_NORMDICT = comm.check_normalisable(START_PSI,
#                                          START_IOFVARS,
#                                          START_ORDER-1,
#                                          orders,
#                                          START_SPLIT_ORDERS,
#                                          update_splits = False)
