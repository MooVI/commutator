Gs = START_GS
sz1 = [comm.Ncproduct(1,[1])]
if START_GS:
    psi_order = comm.unitary_transform_to_order(sz1, Gs, START_ORDER-1)
else:
    psi_order = sz1[:]

for test_order in range(START_ORDER, END_ORDER+1):
    Hcomm = c(small, psi_order)
    psi_order = comm.unitary_transform_to_order(sz1, Gs, test_order, not_single_comm = True)
    Hcomm = comm.add_groups(c(Jpart, psi_order), Hcomm)
    subspace, matrixrows = comm.sparse_find_subspace(Hcomm, Jpart)
    comm.print_subspace(subspace)
    cvector = comm.build_vector(Hcomm, subspace)
    del Hcomm
    iofvars = []
    g_comm = comm.sparse_linear_solve(cvector,
                                matrixrows,
                                subspace,
                                iofvars = iofvars,
                                fvarname = 'F' + str(test_order) + '_')

    print('\n')
    #print(g_comm)
    comm.invert_sz1_G(g_comm, Gs, iofvars)

    #orders.update(zip(iofvars,[test_order]*len(iofvars)))
    psi_order += comm.premultiply(I, comm.commute_group(sz1, Gs[-1]))

    comm.save_group(Gs, FILEHEAD + '_r' + str(test_order))
    comm.save_group(psi_order, FILEHEAD + '_psi_r' + str(test_order))


print('Done!')
