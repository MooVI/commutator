psi = START_PSI
iofvars = START_IOFVARS
split_orders = START_SPLIT_ORDERS
normdict = START_NORMDICT

norm_fail = False

try:
    FAIL_ON_NO_NORM
except:
    FAIL_ON_NO_NORM = True

for test_order in range(START_ORDER, END_ORDER+1):
    Hcomm = c(H, psi)
    if not comm.check_group_at_least_order(Hcomm, test_order-1, orders):
        raise ValueError('Psi does not to order '+ str(test_order-1) + '!')
    subspace, matrixrows = comm.sparse_find_subspace(Hcomm, Jpart)


    comm.print_subspace(subspace)
    cvector = comm.build_vector_to_cancel(comm.substitute_group(Hcomm, normdict), subspace)
    del Hcomm
    psi = comm.substitute_group(psi, normdict, split_orders = split_orders)
    del normdict
    iofvars = []
    psi_test = comm.sparse_solve_for_commuting_term(cvector,
                                                    psi,
                                                    test_order,
                                                    orders,
                                                    matrixrows,
                                                    subspace,
                                                    norm = False,
                                                    iofvars = iofvars,
                                                    split_orders = split_orders,
                                                    fvarname = 'F' + str(test_order) + '_')

    print('\n')

    orders.update(zip(iofvars,[test_order]*len(iofvars)))
    psi += psi_test
    try:
        normdict = comm.check_normalisable(psi, iofvars, test_order, orders, split_orders)
        for x in sorted(normdict.keys(), key = lambda x: int(str(x)[2+len(str(test_order)):])):
            print(str(x)+': ' +str(normdict[x]))
    except Exception as e:
        print(str(e))
        norm_fail = True
        normdict = {}


    comm.save_group(psi,
                    FILEHEAD + '_r' + str(test_order),
                    iofvars=iofvars,
                    split_orders=split_orders,
                    normdict=normdict)

    if(norm_fail and FAIL_ON_NO_NORM):
        raise ValueError("Not normalisable")

    iofvars = []

    if NORM_AS_YOU_GO:
        prop = comm.square_to_find_identity_scalar_up_to_order(psi, test_order, split_orders)
        with open(FILEHEAD+'_norm', mode = 'w') as f:
            f.write(str(prop)+'\n\n')
            f.write(latex(prop).replace('\\\\', '\\'))


if not NORM_AS_YOU_GO:
    prop = comm.square_to_find_identity_scalar_up_to_order(psi, test_order, split_orders)
    with open(FILEHEAD+'_norm', mode = 'w') as f:
        f.write(str(prop)+'\n\n')
        f.write(latex(prop).replace('\\\\', '\\'))

print('Done!')
