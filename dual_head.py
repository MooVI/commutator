import commutator as comm
from sympy import I, symbols, latex

p = comm.print_group
c = comm.calculate_commutator
N = comm.Ncproduct
Ncproduct = N

V1, V2, J,f = symbols('V1 V2 J f')
L = 20
orders = {V1:1, V2:1, f:1}
p.orders = orders
fpart = [Ncproduct(I*f, [2*j+1,2*j+2]) for j in range(L-1)]
V1part = [Ncproduct(V1, [2*j+1,2*j+2,2*j+3,2*j+4]) for j in range(L-2)]
V2part = [Ncproduct(V2, [2*j+2,2*j+3,2*j+4,2*j+5]) for j in range(L-2)]
small = fpart+V1part+V2part
Jpart = [Ncproduct(I*J, [2*j+2,2*j+3]) for j in range(L-1)]
H = fpart + V1part + V2part +Jpart

START_ORDER = 7
END_ORDER = 8

START_PSI = (N(1, 'a1')
             + N(f/J, 'a2') + N(V1/(I*J), 'b1 a2 a3'))

START_IOFVARS = []
START_SPLIT_ORDERS = [0, 1, 3]
START_NORMDICT = {}

START_PSI = comm.load_group('psidual_r6', START_IOFVARS, START_SPLIT_ORDERS)#, START_NORMDICT)

orders.update(zip(START_IOFVARS,[START_ORDER-1]*len(START_IOFVARS)))

FILEHEAD = 'psidual'
NORM_AS_YOU_GO = False

START_NORMDICT = comm.check_normalisable(START_PSI,
                                         START_IOFVARS,
                                         START_ORDER-1,
                                         orders,
                                         START_SPLIT_ORDERS,
                                         update_splits = False)
