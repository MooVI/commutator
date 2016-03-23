import commutator as comm
from sympy import I, symbols, latex

p = comm.print_group
c = comm.calculate_commutator
N = comm.Ncproduct


V,J,f = symbols('V J f')
L = 20
orders = {V:1,f:1}
p.orders = orders
fpart = [N(I*f, [2*j+1,2*j+2]) for j in range(L)]
Vpart = [N(V, [2*j+1,2*j+2,2*j+3,2*j+4]) for j in range(L-2)]
small = fpart+Vpart
Jpart = [N(I*J, [2*j+2,2*j+3]) for j in range(L-1)]
H = fpart + Vpart +Jpart

START_ORDER = 2
END_ORDER = 5

START_PSI = (N(1, 'a1')
             + N(f/J, 'a2') + N(V/(I*J), 'b1 a2 a3'))

START_IOFVARS = []
START_SPLIT_ORDERS = [0, 1, 3]

#START_PSI = comm.load_group('speedpsi3_r6', START_IOFVARS, START_SPLIT_ORDERS)

orders.update(zip(START_IOFVARS,[START_ORDER-1]*len(START_IOFVARS)))

FILEHEAD = 'testpsi'
NORM_AS_YOU_GO = True

START_NORMDICT = comm.check_normalisable(START_PSI,
                                         START_IOFVARS,
                                         START_ORDER-1,
                                         orders,
                                         START_SPLIT_ORDERS,
                                         update_splits = False)
