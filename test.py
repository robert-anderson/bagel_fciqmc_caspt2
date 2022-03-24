from main import *

molecule = {
    'title': 'molecule',
    'geometry': [
        {'xyz': [0.0, 0.0, 0.0], 'atom': 'N'},
        {'xyz': [1.4, 0.0, 0.0], 'atom': 'N'}
    ],
    'basis': '/users/k1507071/code/bagel/src/basis/svp.json',
    'df_basis': '/users/k1507071/code/bagel/src/basis/svp-jkfit.json',
    'angstrom': True
}

trel = False
ncas, ncore = 6, 4

bagel_numthreads, bagel_mpiranks = 4, 1

bagel(exact_caspt2_config(molecule, trel, ncas, ncore), 'exact', numthreads=bagel_numthreads, mpiranks=bagel_mpiranks)

neci_mpiranks = 12
neci_data = deepcopy(default_neci_data)
neci_data['nw'] = 1e5
neci_data['rdm']['mrpt']=None
neci_data['rdm']['samplingiters'] = 200
neci_data['rdm']['startiter'] = 500
neci_data['twf'] = None
neci_data['ss']['popscore'] = 100
neci_data['ss']['startiter'] = 1

niter_casscf = 15
conv_grad = 1e-6

auto_init(molecule, neci_data, neci_mpiranks, bagel_numthreads, bagel_mpiranks, trel, ncas, ncore)
for i in range(niter_casscf):
    grad = auto_iter(neci_mpiranks, bagel_numthreads, bagel_mpiranks)
    if abs(grad)<conv_grad:
        print('CASSCF CONVERGED WITHIN {}'.format(conv_grad))
        break
neci_data['nw'] = 333333
auto_iter(neci_mpiranks, bagel_numthreads, bagel_mpiranks, neci_data=neci_data)

auto_fockmat(neci_mpiranks, bagel_numthreads, bagel_mpiranks, i+1)
auto_caspt2(neci_mpiranks, bagel_numthreads, bagel_mpiranks)
