import commutator as comm
from sympy import I, symbols, latex

p = comm.print_group
c = comm.calculate_commutator
N = comm.Ncproduct


X, Y = symbols('X Y', real=True)
comm.mathematica_parser.vardict = {'X': X, 'Y':Y}
L = 12
orders = {X:1,Y:1}
p.orders = orders
Ypart = [N(I*Y, [2*j+1,2*j+4]) for j in range(L-2)]
Xpart = [N(-X, [2*j+1,2*j+2,2*j+3,2*j+4]) for j in range(L-2)]
small = Xpart+Ypart
Zpart = [N(-I, [2*j+2,2*j+3]) for j in range(L-1)]
Jpart = Zpart
H = small +Jpart

START_ORDER = 1
END_ORDER = 20

START_PSI = N(1, 'a1')
START_GS = []

START_IOFVARS = []
START_SPLIT_ORDERS = [0, 1]
START_NORMDICT = {}

#START_PSI = comm.load_group('xyzfast_r5', START_IOFVARS, START_SPLIT_ORDERS, START_NORMDICT)

orders.update(zip(START_IOFVARS,[START_ORDER-1]*len(START_IOFVARS)))

FILEHEAD = 'xyzfast_L_12_gs'
NORM_AS_YOU_GO = True

# START_NORMDICT = comm.check_normalisable(START_PSI,
#                                         START_IOFVARS,
#                                         START_ORDER-1,
#                                         orders,
#                                         START_SPLIT_ORDERS,
#                                         update_splits = False)
