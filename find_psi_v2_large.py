import commutator as comm
from sympy import I, symbols, latex

p = comm.print_group
c = comm.calculate_commutator
N = comm.Ncproduct
Ncproduct = N

V1, J2, J,f = symbols('V1 J2 J f')
L = 20
orders = {V1:1, f:1}
p.orders = orders
fpart = [Ncproduct(I*f, [2*j+1,2*j+2]) for j in range(L)]
V1part = [Ncproduct(V1, [2*j+1,2*j+2,2*j+3,2*j+4]) for j in range(L-2)]
J2part = [Ncproduct(J2, [2*j+2,2*j+3,2*j+4,2*j+5]) for j in range(L-2)]
small = fpart+V1part
J1part = [Ncproduct(I*J, [2*j+2,2*j+3]) for j in range(L-1)]
Jpart =  J1part + J2part
H = fpart +Jpart+V1part

START_ORDER = 1
END_ORDER = 2

START_PSI = N(1, 'a1')

START_IOFVARS = []
START_SPLIT_ORDERS = [0, 1]

#START_PSI = comm.load_group('speedpsi3_r6', START_IOFVARS, START_SPLIT_ORDERS)

orders.update(zip(START_IOFVARS,[START_ORDER-1]*len(START_IOFVARS)))

FILEHEAD = 'psi_v2_large'
NORM_AS_YOU_GO = True


psi = START_PSI
iofvars = START_IOFVARS
split_orders = START_SPLIT_ORDERS

for test_order in range(START_ORDER, END_ORDER+1):
    psi_sub = None
    psi_test_sub = None

    Hcomm = c(H, psi)
    if not comm.check_group_at_least_order(Hcomm, test_order-1, orders):
        raise ValueError('Psi does not to order '+str(test_order-1)+'!')
    subspace, matrixrows = comm.sparse_find_subspace(Hcomm, Jpart)
    comm.print_subspace(subspace)
    cvector = comm.build_vector_to_cancel(Hcomm, subspace)
    subs_rules = {}
    psi_test = comm.sparse_solve_for_commuting_term(cvector,
                                                    psi,
                                                    test_order,
                                                    orders,
                                                    matrixrows,
                                                    subspace,
                                                    norm = False,
                                                    subs_rules = subs_rules,
                                                    iofvars = iofvars,
                                                    split_orders = split_orders,
                                                    fvarname = 'F' + str(test_order) + '_')

    for x in sorted(subs_rules.keys(), key = lambda x: int(str(x)[2+len(str(test_order)):])):
        print(str(x)+': ' +str(subs_rules[x]))
    print('\n')
    psi = comm.substitute_group(psi, subs_rules, split_orders)

    orders.update(zip(iofvars,[test_order]*len(iofvars)))
    normdict = comm.check_normalisable(psi+psi_test, iofvars, test_order, orders, split_orders)
    for x in sorted(normdict.keys(), key = lambda x: int(str(x)[2+len(str(test_order)):])):
        print(str(x)+': ' +str(normdict[x]))
    psi_test_sub = comm.substitute_group(psi_test, normdict)
    psi_sub = psi + psi_test_sub
    comm.save_group(psi_sub, FILEHEAD)

    psi += psi_test
    comm.save_group(psi,
                    FILEHEAD + '_r' + str(test_order), iofvars=iofvars, split_orders=split_orders)

    if NORM_AS_YOU_GO:
        prop = comm.square_to_find_identity(psi_sub)[0].scalar
        with open(FILEHEAD+'_norm', mode = 'w') as f:
            f.write(str(prop)+'\n\n')
            f.write(latex(prop).replace('\\\\', '\\'))



if not NORM_AS_YOU_GO:
    prop = comm.square_to_find_identity(psi_sub)[0].scalar
    with open(FILEHEAD+'_norm', mode = 'w') as f:
        f.write(str(prop)+'\n\n')
        f.write(latex(prop).replace('\\\\', '\\'))

print('Done!')