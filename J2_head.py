import commutator as comm
from sympy import I, symbols, latex

p = comm.print_group
c = comm.calculate_commutator
N = comm.Ncproduct
Ncproduct = N

J2, f = symbols('J2 f', real=True)
L = 40
orders = {f:1, J2:1}
comm.mathematica_parser.vardict = {'J2':J2, 'f': f}
p.orders = orders
fpart = [Ncproduct(I*f, [2*j+1,2*j+2]) for j in range(L-1)]
V2part = [Ncproduct(J2, [2*j+2,2*j+3,2*j+4,2*j+5]) for j in range(L-2)]
small = fpart+V2part
Jpart = [Ncproduct(I, [2*j+2,2*j+3]) for j in range(L-1)]
H = small +Jpart

START_ORDER = 1
END_ORDER = 14

START_PSI = [N(1, 'a1')]

START_IOFVARS = []
START_SPLIT_ORDERS = [0, 1]
START_NORMDICT = {}

#START_GS = comm.load_group('j2_smalluni_r8.yaml')

START_GS = []

#START_PSI = comm.load_group('psi_j2_matlab_large_r10', START_IOFVARS, START_SPLIT_ORDERS, START_NORMDICT)

orders.update(zip(START_IOFVARS,[START_ORDER-1]*len(START_IOFVARS)))

FILEHEAD = 'j2_small_gs'
NORM_AS_YOU_GO = True

# START_NORMDICT = comm.check_normalisable(START_PSI,
#                                          START_IOFVARS,
#                                          START_ORDER-1,
#                                          orders,
#                                          START_SPLIT_ORDERS,
#                                          update_splits = False)
