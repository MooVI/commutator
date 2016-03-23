import commutator as comm
from sympy import I, symbols, latex

p = comm.print_group
c = comm.calculate_commutator
N = comm.Ncproduct


X, Y, f, g = symbols('X Y f g')
L = 20
orders = {X:1,Y:1,f:1,g:1}
p.orders = orders
Ypart = [N(I*Y, [2*j+1,2*j+4]) for j in range(L-2)]
Xpart = [N(-X, [2*j+1,2*j+2,2*j+3,2*j+4]) for j in range(L-2)]
linpart = [N(I*f, [2*j+1,2*j+2]) for j in range(L-1)]
lin2part = [N(-g, [2*j+1,2*j+2,2*j+3,2*j+4,2*j+5]) for j in range(L-2)]
small = Xpart+Ypart+linpart+lin2part
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

FILEHEAD = 'xyzpypx'
NORM_AS_YOU_GO = True

# START_NORMDICT = comm.check_normalisable(START_PSI,
#                                         START_IOFVARS,
#                                         START_ORDER-1,
#                                         orders,
#                                         START_SPLIT_ORDERS,
#                                         update_splits = False)
