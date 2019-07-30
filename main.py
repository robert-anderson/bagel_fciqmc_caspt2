import pickle, os
from copy import deepcopy
import tempfile
import shutil
import json
import time
import re
from subprocess import Popen, PIPE

def get_line_metadata(line, token, type_=int):
    if type_==int:
        regex = re.compile(r'{}\s?=\s?[0-9]+'.format(token))
    elif type_==bool:
        regex = re.compile(r'{}\s?=\s?\..+\.'.format(token))
    else:
        assert 0, 'Unsupported type'
    regex = regex.search(line)
    if regex is None: return None
    if type_ is bool: return regex.group().split('=')[1].strip()=='.TRUE.'
    return type_(regex.group().split('=')[1].strip())

def get_fcidump_metadata(token, type_=int, fname='FCIDUMP'):
    with open(fname, 'r') as f:
        for line in f.readlines():
            tmp = get_line_metadata(line, token, type_)
            if tmp is not None: return tmp

def extract(lines, anchor_pattern, split_offset, line_offset=0, type_=float):
    if isinstance(lines, str) and '\n' in lines: lines = lines.split('\n')
    for i, line in enumerate(lines):
        if re.search(anchor_pattern, line):
            return type_(lines[i+line_offset].split()[split_offset])

def file_extract(fname, anchor_pattern, split_offset, line_offset=0, type_=float):
    with open(fname, 'r') as f: 
        return extract(f.readlines(), anchor_pattern, split_offset, line_offset, type_)

def storeobj(obj, fname):
    f = tempfile.NamedTemporaryFile(dir='.', delete=False)
    pickle.dump(obj, f)
    shutil.move(f.name, fname)

def loadobj(fname):
    with open(fname, 'rb') as f: return pickle.load(f)

def link(src, dst):
    assert os.path.exists(src)
    shell('rm {1}; ln -s {0} {1}'.format(os.path.abspath(src), os.path.abspath(dst)))

def copy(src, dst):
    assert os.path.exists(src)
    shell('cp {0} {1}'.format(os.path.abspath(src), os.path.abspath(dst)))

def link_file(fname, srcdir, dstdir):
    link('{}/{}'.format(srcdir, fname), '{}/{}'.format(dstdir, fname))

def copy_file(fname, srcdir, dstdir):
    copy('{}/{}'.format(srcdir, fname), '{}/{}'.format(dstdir, fname))

def iterable_or_numeric_to_string(obj, delimiter=' '):
    try: return str(obj+0)
    except ValueError: return ' '.join(tuple(map(str, self['rdm']['mrpt']['promotionfractions'])))
ionts = iterable_or_numeric_to_string

class StructuredFile(list):
    def __init__(self, data):
        self.data = data
        list.__init__(self)
        self.body()
    def __getitem__(self, key):
        if isinstance(key, int): return list.__getitem__(self, key)
        else: return self.data[key]
        # it's convenient to refer to elements of self.data directly from self
    def __lshift__(self, string):
        self.append(string)
        return self
    def body(self):
        pass
    def render(self, fname):
        with open(fname, 'w') as f: f.write('\n'.join(self))

class NeciInputFile(StructuredFile):
    def body(self):
        self.extend([
        'title',
        'system read noorder',
        'electrons {}'.format(self['nelec']),
        'symignoreenergies',
        'freeformat',
        'sym 0 0 0 0',
        'nobrillouintheorem',
        'nonuniformrandexcits {}'.format('4ind-weighted' if not self['trel'] else 'pick-uniform-mag'),
        'endsys',
        'calc',
        'methods',
        'method vertex fcimc',
        'endmethods',
        'totalwalkers {}'.format(self['nw']),
        'memoryfacpart {}'.format(self['memoryfacpart']),
        'memoryfacspawn {}'.format(self['memoryfacspawn']),
        'seed {}'.format(self['seed']),
        'startsinglepart'.format(self['startsinglepart']),
        'shiftdamp {}'.format(self['shiftdamp']),
        'truncinitiator',
        'addtoinitiator {}'.format(self['addtoinitiator']),
        'allrealcoeff',
        'realspawncutoff 0.4',
        'jump-shift',
        'stepsshift {}'.format(self['stepsshift']),
        'maxwalkerbloom 3',
        'load-balance-blocks off',
        'nmcyc {}'.format(self['nmcyc'])
        ])
        if self['rdm']:
            self << 'rdmsamplingiters {}'.format(self['rdm']['samplingiters'])
        if self['ss']:
            self << 'semi-stochastic {}'.format(self['ss']['startiter'])
            if self['ss']['readcore']:
                # initialise from previous semistochastic calculation
                self << 'readcore'
            elif self['ss']['fcicore']:
                # initialise entire from previous semistochastic calculation
                self << 'fci-core'
            else:
                # otherwise we need the number of core determinant to be specified
                self << 'pops-core {}'.format(self['ss']['popscore'])
        if self['twf']:
            self << 'trial-wavefunction {}'.format(self['twf']['startiter'])
            self << 'pops-trial {}'.format(self['twf']['popstrial'])
        self.extend([
        'endcalc',
        'integral',
        'freeze 0 0',
        'endint',
        'logging'
        ])
        if self['ss']:
           if self['ss']['writecore']:
               self << 'write-core'
        if self['rdm']:
            self << 'calcrdmonfly 3 {} {}'.format(self['rdm']['startiter'], self['rdm']['energyiter'])
            if self['rdm']['mrpt']:
                self << self['rdm']['mrpt']['type']
                self << 'hbrdm-ghost-granularity {}'.format(self['rdm']['mrpt']['granularity'])
                self << '4rdm-promotion-fractions {}'.format(ionts(self['rdm']['mrpt']['promotionfractions']))
            if not self['trel']: self << 'write-spin-free-rdm'
            self << 'printonerdm'
            self << 'rdm-main-size-fac {}'.format(ionts(self['rdm']['mainsizefac']))
            self << 'rdm-spawn-size-fac {}'.format(ionts(self['rdm']['spawnsizefac']))
            self << 'rdm-recv-size-fac {}'.format(ionts(self['rdm']['recvsizefac']))
        self << 'binarypops'
        self << 'endlog'
        self << 'end'

default_neci_data = {
    'memoryfacpart': 10.0,
    'memoryfacspawn': 20.0,
    'seed': 14,
    'startsinglepart': 100,
    'nw': 10000,
    'shiftdamp': 0.05,
    'addtoinitiator': 3.0,
    'stepsshift': 5,
    'nmcyc': -1,
    'rdm':{
        'samplingiters': 100,
        'startiter': 0,
        'energyiter': 50,
        'mrpt': {
            'type': 'nevpt2',
            'granularity': 1,
            'promotionfractions': 1,
        },
        'mainsizefac': 1,
        'spawnsizefac': 1,
        'recvsizefac': 1,
    },
    'ss':{
        'startiter': 0,
        'readcore': None,
        'writecore': None,
        'fcicore': None,
        'popscore': 10,
        'energyiter': 50,
    },
    'twf':{
        'startiter': 0,
        'popstrial': 10
    },
    'trel': False
}



class BagelInputFile(dict):
    def __init__(self, molecule):
        dict.__init__(self)
        self['bagel'] = []
        if 'title' not in molecule.keys(): molecule['title']='molecule'
        self['bagel'].append(molecule)

class MeanFieldBagelInputFile(BagelInputFile):
    def __init__(self, molecule, trel=False, gaunt=None, breit=None, thresh=1e-10, maxiter=100, robust=True):
        self.molecule, self.trel, self.gaunt, self.breit = molecule, trel, gaunt, breit
        if not trel: assert not gaunt and not breit
        BagelInputFile.__init__(self, molecule)
        tmp = {'thresh':thresh, 'maxiter':maxiter, 'robust':robust}
        if not trel:
            tmp['title'] = 'hf'
        else:
            tmp['title'] = 'dhf'
            tmp['gaunt'] = False if gaunt is None else gaunt
            tmp['breit'] = False if breit is None else breit
        self['bagel'].append(tmp)

class CasscfBagelInputFile(MeanFieldBagelInputFile):
    def __init__(self, molecule, norb, ncore, trel=False, gaunt=None, breit=None, thresh=1e-10, maxiter=100, robust=True, canonical=False):
        MeanFieldBagelInputFile.__init__(self, molecule, trel, gaunt, breit, thresh, maxiter, robust)
        self['bagel'].append({
            'thresh': thresh, 
            'maxiter': maxiter,
            'nact': norb,
            'nclosed': ncore,
            'canonical': canonical,
            'robust': robust
        })
        if trel:
            self['bagel'][-1]['title'] = 'zcasscf'
            self['bagel'][-1]['state'] = [1]
        else:
            self['bagel'][-1]['title'] = 'casscf'
            self['bagel'][-1]['nstate'] = [1]

class InitCasscfBagelInputFile(CasscfBagelInputFile):
    def __init__(self, molecule, norb, ncore, trel=False, gaunt=None, breit=None, thresh=1e-10, maxiter=100, robust=True, canonical=False):
        CasscfBagelInputFile.__init__(self, molecule, norb, ncore, trel, gaunt, breit, thresh, maxiter, robust, canonical)
        self['bagel'][0]['external_rdm'] = 'noref'
        self['bagel'].append({
            "only_ints" : True,
            "ncore" : ncore,
            "norb" :  norb,
            "frozen" : False,
            "state" : [1]
        })
        if trel:
            self['bagel'][-1]['title'] = 'zfci'
            self['bagel'][-1]['state'] = [1]
        else:
            self['bagel'][-1]['title'] = 'fci'
            self['bagel'][-1]['nstate'] = 1

class Caspt2BagelInputFile(MeanFieldBagelInputFile):
    def __init__(self, molecule, norb, ncore, trel=False, gaunt=None, breit=None, thresh=1e-10, maxiter=100, robust=True, canonical=False):
        CasscfBagelInputFile.__init__(self, molecule, norb, ncore, trel, gaunt, breit, thresh, maxiter, robust, canonical)
        self['bagel'].append({
            'method': 'caspt2',
            'frozen': False
        })
        self['bagel'][-1]['title'] = 'relsmith' if trel else 'smith'

def shell(*lines):
    if len(lines): lines = ';'.join(lines)
    start = time.time()
    o, e = Popen(lines, shell=1, stdout=PIPE, stderr=PIPE).communicate()
    end = time.time()
    return o, e, end-start

class NeciInstance:
    infile = 'neci.inp'
    outfile = 'neci.out'
    pklfile = 'neci.pkl'
    def __init__(self, input_data, bagel_dirname='.', dirname=None):
        if dirname is None: dirname = bagel_dirname
        assert not os.path.exists(self.pklfile), \
                'NECI working directory must not contain a previous calculation'
        print(input_data)
        self.input_data = deepcopy(input_data)
        self.norb  = get_fcidump_metadata('NORB')
        self.nelec = get_fcidump_metadata('NELEC')
        self.trel  = get_fcidump_metadata('TREL', bool)
        self.input_data['nelec'] = self.nelec
        self.input_data['trel'] = self.trel
        self.input_file = NeciInputFile(self.input_data)

        self.dirname = dirname
        self.infile = os.path.abspath('{}/{}'.format(self.dirname, self.infile))
        self.outfile = os.path.abspath('{}/{}'.format(self.dirname, self.outfile))
        self.bagel_dirname = bagel_dirname
        prev_bagel_obj = loadobj('{}/bagel.pkl'.format(bagel_dirname))
    def execute(self, neci_bin_dirname='/users/k1507071/code/neci_merge/build_ib/bin', mpiranks=1):
        self.input_file.render(self.infile)
        if not os.path.exists('{}/FCIDUMP'.format(self.dirname)):
            link_file('FCIDUMP', self.bagel_dirname, self.dirname)
        print('Invoking NECI...')
        stdout, stderr, time = shell(
            'cd {}'.format(self.dirname),
            'module purge',
            'module load /users/k1507071/opt/modules/gcc/7.4.0',
            'module load /users/k1507071/opt/modules/openmpi/4.0.1/gnu_7.4.0',
            'mpirun -np {} {}/{} {} > {}'.format(
                mpiranks, neci_bin_dirname, 'kdneci' if self.trel else 'dneci', self.infile, self.outfile)
        )
        print(stderr)
        print('NECI instance completed ({:.2f} seconds)'.format(time))
        self.duration = time
        storeobj(self, self.pklfile)

class BagelInstance:
    infile = 'bagel.json'
    outfile = 'bagel.out'
    pklfile = 'bagel.pkl'
    def __init__(self, input_file, dirname='.'):
        self.input_file = input_file
        self.dirname = dirname
        self.infile = os.path.abspath('{}/{}'.format(self.dirname, self.infile))
        self.outfile = os.path.abspath('{}/{}'.format(self.dirname, self.outfile))
        self.pklfile = os.path.abspath('{}/{}'.format(self.dirname, self.pklfile))
        assert not os.path.exists(self.infile), \
                'BAGEL working directory must not contain a previous calculation'
    def execute(self, bagel_exe='/users/k1507071/code/rja_bagel/obj/src/BAGEL', numthreads=1, mpiranks=1):
        with open(self.infile, 'w') as f: json.dump(self.input_file, f)
        print('Invoking BAGEL...')
        stdout, stderr, time = shell(
            'cd {}'.format(self.dirname),
            'export OMP_NUM_THREADS={}'.format(numthreads),
            'module purge',
            'module load /users/k1507071/opt/modules/gcc/7.4.0',
            'module load /users/k1507071/opt/modules/openmpi/4.0.1/gnu_7.4.0',
            'module load /users/k1507071/opt/modules/boost/1.70/gcc7.4.0/python2.7',
            '. $HOME/intel/mkl/bin/mklvars.sh intel64',
            'mpirun -np {} {} {} > {}'.format(mpiranks, bagel_exe, self.infile, self.outfile)
        )
        print('BAGEL instance completed ({:.2f} seconds)'.format(time))
        self.duration = time
        storeobj(self, self.pklfile)

class MeanFieldBagelInstance(BagelInstance):
    def __init__(self, molecule, trel=False, gaunt=None, breit=None, thresh=1e-10, maxiter=100, robust=True, dirname='.'):
        input_file = MeanFieldBagelInputFile(molecule, trel, gaunt, breit, thresh, maxiter, robust)
        BagelInstance.__init__(self, input_file, dirname)
    def get_energy(self, bagel_exe='/users/k1507071/code/rja_bagel/obj/src/BAGEL', numthreads=1, mpiranks=1):
        if not os.path.exists(self.outfile):
            self.execute(bagel_exe, numthreads, mpiranks)
        return file_extract(self.outfile, 'SCF iteration converged', 1, -2)

class CasscfBagelInstance(MeanFieldBagelInstance):
    def __init__(self, molecule, norb, ncore, trel=False, gaunt=None, breit=None, thresh=1e-10, maxiter=100, robust=True, canonical=False, dirname='.'):
        input_file = CasscfBagelInputFile(molecule, norb, ncore, trel, gaunt, breit, thresh, maxiter, robust, canonical)
        BagelInstance.__init__(self, input_file, dirname)
    def get_energy(self, bagel_exe='/users/k1507071/code/rja_bagel/obj/src/BAGEL', numthreads=1, mpiranks=1):
        if not os.path.exists(self.outfile):
            self.execute(bagel_exe, numthreads, mpiranks)
        return file_extract(self.outfile, 'optimization converged', 2, -2)
    def get_gradient(self, bagel_exe='/users/k1507071/code/rja_bagel/obj/src/BAGEL', numthreads=1, mpiranks=1):
        if not os.path.exists(self.outfile):
            self.execute(bagel_exe, numthreads, mpiranks)
        return file_extract(self.outfile, 'optimization converged', 3, -2)

class InitCasscfBagelInstance(CasscfBagelInstance):
    def __init__(self, molecule, norb, ncore, trel=False, gaunt=None, breit=None, thresh=1e-10, maxiter=100, robust=True, canonical=False, dirname='.'):
        input_file = InitCasscfBagelInputFile(molecule, norb, ncore, trel, gaunt, breit, thresh, maxiter, robust, canonical)
        BagelInstance.__init__(self, input_file, dirname)

class CasscfNeciInstance(NeciInstance):
    def __init__(self, input_data, dirname='.'):
        assert not os.path.exists('{}/neci.pkl'.format(dirname)), \
                'BAGEL working directory must not contain a previous calculation'

class IterCasscfBagelInstance(CasscfBagelInstance):
    def __init__(self, prev_bagel_dirname, prev_neci_dirname=None):
        assert os.path.exists(prev_bagel_dirname), \
                'BAGEL CASSCF iteration must be based on a previous instance'
        assert os.path.exists('{}/bagel.pkl'.format(prev_bagel_dirname)),\
                'BAGEL CASSCF iteration must be based on a previously executed instance'
        prev_bagel_obj = loadobj('{}/bagel.pkl'.format(prev_bagel_dirname))
        if self.prev_bagel_obj.trel:
            assert os.path.exists('{}/relref.archive'.format(prev_bagel_dirname)),\
                'Previous relativistic BAGEL instance lacks a serialised state'
            link_file('relref.archive', prev_bagel_dirname, dirname)
        else:
            assert os.path.exists('{}/ref.archive'.format(prev_bagel_dirname)),\
                'Previous nonrelativistic BAGEL instance lacks a serialised state'
            link_file('ref.archive', prev_bagel_dirname, dirname)
        if prev_neci_dirname is None: prev_neci_dirname=prev_bagel_dirname
        assert os.path.exists(prev_neci_dirname), \
                'NECI CASSCF iteration must be based on a previous instance'
        assert os.path.exists('{}/neci.pkl'.format(prev_neci_dirname)),\
                'NECI CASSCF iteration must be based on a previously executed instance'
        if self.prev_bagel_obj.trel:
            assert os.path.exists('{}/1RDM.1'.format(prev_neci_dirname)),\
                'Previous NECI instance lacks a relativistic 1RDM file'
            link('{}/1RDM.1'.format(prev_neci_dirname), '{}/fciqmc_0_0.rdm1'.format(dirname))
            assert os.path.exists('{}/2RDM.1'.format(prev_neci_dirname)),\
                'Previous NECI instance lacks a relativistic 2RDM file'
            link('{}/2RDM.1'.format(prev_neci_dirname), '{}/fciqmc_0_0.rdm2'.format(dirname))
        else:
            assert os.path.exists('{}/spinfree_OneRDM.1'.format(prev_neci_dirname)),\
                'Previous NECI instance lacks a spinfree 1RDM file'
            link('{}/spinfree_OneRDM.1'.format(prev_neci_dirname), '{}/fciqmc_0_0.rdm1'.format(dirname))
            assert os.path.exists('{}/spinfree_TwoRDM.1'.format(prev_neci_dirname)),\
                'Previous NECI instance lacks a spinfree 2RDM file'
            link('{}/spinfree_TwoRDM.1'.format(prev_neci_dirname), '{}/fciqmc_0_0.rdm2'.format(dirname))

        input_file = IterCasscfBagelInputFile(molecule, norb, ncore, trel, gaunt, breit, thresh, maxiter, robust, canonical)
        BagelInstance.__init__(self, input_file, dirname)

class Caspt2BagelInstance(CasscfBagelInstance):
    def __init__(self, molecule, norb, ncore, trel=False, gaunt=None, breit=None, thresh=1e-10, maxiter=100, robust=True, canonical=False, dirname='.'):
        input_file = Caspt2BagelInputFile(molecule, norb, ncore, trel, gaunt, breit, thresh, maxiter, robust, canonical)
        BagelInstance.__init__(self, input_file, dirname)
    def get_energy(self, bagel_exe='/users/k1507071/code/rja_bagel/obj/src/BAGEL', numthreads=1, mpiranks=1):
        if not os.path.exists(self.outfile):
            self.execute(bagel_exe, numthreads, mpiranks)
        return file_extract(self.outfile, 'CASPT2 energy evaluation', 6, 2)

molecule = {
    'geometry': [
        {'xyz': [0.0,    0.0, 0.0], 'atom': 'Be'},
        {'xyz': [1.6788, 0.0, 0.0], 'atom': 'Be'}
    ],
    'basis': '/users/k1507071/code/bagel/src/basis/tzvpp.json',
    'df_basis': '/users/k1507071/code/bagel/src/basis/tzvpp-jkfit.json',
    'angstrom': True
}



tmp = InitCasscfBagelInstance(molecule, 4, 2, canonical=1)
print(tmp.get_energy())
tmp = NeciInstance(default_neci_data, bagel_dirname='.')
tmp.execute()

assert 0

class SolverInstance:
    pass

class EnvConfig:
    pass

class Calc:
    def __init__(self, low_wn, high_wn, chain_pops=False):
        # walker number for first iterations
        self.low_wn = low_wn
        # walker number for tight-convergence
        self.high_wn = high_wn
        # should the WF POPSFILES be chained together in the tight-convergence regime?
        self.chain_pops = chain_pops
    def casscf(self, wn, ediff_thresh=1e-6, grad_thresh=1e-5, max_iters=30, chain_pops=False):
        assert self.mean_field_data is not None



calc = Calc()

calc.setup()
calc.mean_field()
calc.casscf(1e5, chain_pops=False)
calc.casscf(wn=1e6, chain_pops=True)
