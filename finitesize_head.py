import commutator as comm
from sympy import I, symbols, latex

p = comm.print_group
c = comm.calculate_commutator
N = comm.Ncproduct


V,J,f = symbols('V J f')
L = 4
orders = {V:1,f:1}
p.orders = orders
fpart = [N(I*f, [2*j+1,2*j+2]) for j in range(L-1)]
Vpart = [N(V, [2*j+1,2*j+2,2*j+3,2*j+4]) for j in range(L-2)]
small = fpart+Vpart
Jpart = [N(I*J, [2*j+2,2*j+3]) for j in range(L-1)]
H = fpart + Vpart +Jpart

START_ORDER = 2
END_ORDER = 35

START_PSI = (N(1, 'a1')
             + N(f/J, 'a2') + N(V/(I*J), 'b1 a2 a3'))

START_IOFVARS = []
START_SPLIT_ORDERS = [0, 1, 3]
START_NORMDICT = {}

#START_PSI = comm.load_group('finitesize_r25', START_IOFVARS, START_SPLIT_ORDERS, START_NORMDICT)

orders.update(zip(START_IOFVARS,[START_ORDER-1]*len(START_IOFVARS)))

FILEHEAD = 'finitesizetrunc'
NORM_AS_YOU_GO = False

# START_NORMDICT = comm.check_normalisable(START_PSI,
#                                         START_IOFVARS,
#                                         START_ORDER-1,
#                                         orders,
#                                         START_SPLIT_ORDERS,
#                                         update_splits = False)
