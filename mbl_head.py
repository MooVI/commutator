import commutator as comm
from sympy import I, symbols, latex, numbered_symbols

p = comm.print_group
c = comm.calculate_commutator
N = comm.Ncproduct
cv = comm.convert_group
Ncproduct = N

V, J, f = symbols('V J f', real=True)
h_gen = numbered_symbols('h', real=True)
L = 20
h = [next(h_gen) for i in range(L+1)]
comm.mathematica_parser.vardict = {'J':J, 'f': f}
comm.mathematica_parser.vardict.update(zip([str(hq) for hq in h],h))
orders = {V:1, f:1}
p.orders = orders
fpart = [Ncproduct(I*f, [2*j+1,2*j+2]) for j in range(L)]
V1part = [Ncproduct(V, [2*j+1,2*j+2,2*j+3,2*j+4]) for j in range(L-2)]
hpart = [h[j+1]*(cv(comm.SigmaProduct(1,2*j+1))[0]) for j in range(L-2)]
small = fpart
J1part = [Ncproduct(I*J, [2*j+2,2*j+3]) for j in range(L-1)]
Jpart =  J1part+hpart
H = small+ Jpart

START_ORDER = 1
END_ORDER = 4

START_PSI = cv(comm.SigmaProduct(1,'z5'))

START_IOFVARS = []
START_SPLIT_ORDERS = [0, 1]
START_NORMDICT = {}

#START_PSI = comm.load_group('psi_j2_matlab_large_r10', START_IOFVARS, START_SPLIT_ORDERS, START_NORMDICT)

orders.update(zip(START_IOFVARS,[START_ORDER-1]*len(START_IOFVARS)))

FILEHEAD = 'mbl_5'
NORM_AS_YOU_GO = True

# START_NORMDICT = comm.check_normalisable(START_PSI,
#                                          START_IOFVARS,
#                                          START_ORDER-1,
#                                          orders,
#                                          START_SPLIT_ORDERS,
#                                          update_splits = False)
