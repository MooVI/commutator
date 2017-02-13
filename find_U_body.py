Gs = START_GS
sz1 = [comm.Ncproduct(1,[1])]
if START_GS:
    psi_order = comm.unitary_transform_to_order(sz1, Gs, START_ORDER-1)
else:
    psi_order = sz1[:]

for test_order in range(START_ORDER, END_ORDER+1):
    Hcomm = c(small, psi_order)
    psi_order = comm.unitary_transform_to_order(sz1, Gs, test_order, not_single_comm = True)
    Hcomm = comm.simplify_group(c(Jpart, psi_order) + Hcomm)
    subspace, matrixrows = comm.sparse_find_subspace(Hcomm, Jpart)
    comm.print_subspace(subspace)
    cvector = comm.build_vector_to_cancel(Hcomm, subspace)
    del Hcomm
    subs_rules = {}
    iofvars = []
    g_comm = comm.sparse_solve_for_commuting_term(cvector,
                                                    None,
                                                    test_order,
                                                    orders,
                                                    matrixrows,
                                                    subspace,
                                                    subs_rules = subs_rules,
                                                    iofvars = iofvars,
                                                    fvarname = 'F' + str(test_order) + '_')

    for x in sorted(subs_rules.keys(), key = lambda x: int(str(x)[2+len(str(test_order)):])):
        print(str(x)+': ' +str(subs_rules[x]))
    print('\n')
    print(g_comm)
    comm.invert_sz1_G(g_comm, Gs, iofvars)

    #orders.update(zip(iofvars,[test_order]*len(iofvars)))
    psi_order += comm.premultiply(I, comm.commute_group(sz1, Gs[-1]))

    comm.save_group(Gs, FILEHEAD + '_r' + str(test_order))


print('Done!')
