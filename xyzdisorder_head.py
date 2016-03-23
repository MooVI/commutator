import commutator as comm
from sympy import I, symbols, latex, numbered_symbols

p = comm.print_group
c = comm.calculate_commutator
N = comm.Ncproduct


X_gen = numbered_symbols('X')
Y_gen = numbered_symbols('Y')
L = 20
X = [next(X_gen) for i in range(L+1)]
Y = [next(Y_gen) for i in range(L+1)]
orders = dict(zip(X+Y,[1]*2*(L+1)))
p.orders = orders
Ypart = [N(I*Y[j+1], [2*j+1,2*j+4]) for j in range(L-2)]
Xpart = [N(-X[j+1], [2*j+1,2*j+2,2*j+3,2*j+4]) for j in range(L-2)]
small = Xpart+Ypart
Zpart = [N(-I, [2*j+2,2*j+3]) for j in range(L-1)]
Jpart = Zpart
H = small +Jpart

START_ORDER = 1
END_ORDER = 8

START_PSI = N(1, 'a1')

START_IOFVARS = []
START_SPLIT_ORDERS = [0, 1]
START_NORMDICT = {}

#START_PSI = comm.load_group('jto1psi2_r8', START_IOFVARS, START_SPLIT_ORDERS, START_NORMDICT)

orders.update(zip(START_IOFVARS,[START_ORDER-1]*len(START_IOFVARS)))

FILEHEAD = 'xyzfastdisorder'
NORM_AS_YOU_GO = True

# START_NORMDICT = comm.check_normalisable(START_PSI,
#                                         START_IOFVARS,
#                                         START_ORDER-1,
#                                         orders,
#                                         START_SPLIT_ORDERS,
#                                         update_splits = False)
