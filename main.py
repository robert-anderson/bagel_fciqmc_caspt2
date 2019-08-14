import pickle, os
from copy import deepcopy
import tempfile
import shutil
import json
import time
import re
from subprocess import Popen, PIPE

def shell(*lines):
    if len(lines): lines = ';'.join(lines)
    start = time.time()
    o, e = Popen(lines, shell=1, stdout=PIPE, stderr=PIPE).communicate()
    end = time.time()
    return o, e, end-start

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

def file_extract(fname, anchor_pattern, split_offset, line_offset=0, type_=float):
    if line_offset<0:
        stdout, stderr, _ = shell('grep "{}" {} -B{} | head -n1'.format(anchor_pattern, fname, -line_offset))
    elif line_offset>0:
        stdout, stderr, _ = shell('grep "{}" {} -A{} | tail -n1'.format(anchor_pattern, fname, line_offset))
    else:
        stdout, stderr, _ = shell('grep "{}" {}'.format(anchor_pattern, fname))
    return type_(str(stdout, 'utf-8').split()[split_offset])

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

def reformat_fockmat(fname):
    EPS = 1e-8
    lines = None
    with open(fname, 'r') as f:
        lines = f.readlines()

    with open(fname, 'w') as f:
        i = 1
        for line in lines:
            j = 1
            for val in line.split():
                try:
                    num_val = float(val)
                except ValueError:
                    num_val = complex(*map(float, val[1:-1].split(',')))
                if j!=i and abs(num_val)>EPS:
                    print("WARNING: large off-diagonal Fock matrix elements")
                if j>=i:
                    f.write('{}    {}    {}\n'.format(val, i, j))
                if j==len(lines):
                    break
                j+=1
            i+=1


def check_hardlink(src, dst):
    if os.path.exists(src):
        os.link(src, dst)

def check_copyfile(src, dst):
    if os.path.exists(src):
        shutil.copyfile(src, dst)
        
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
        if self['ss'] is not None:
            self << 'semi-stochastic {}'.format(self['ss']['startiter'] if self['ss']['startiter'] else '')
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

def hf_block(trel, thresh=1e-10, maxiter=100, robust=False, gaunt=None, breit=None):
    block = {
        'thresh':thresh,
        'maxiter':maxiter,
        'robust':robust
    }
    if not trel:
        block['title'] = 'hf'
    else:
        block['title'] = 'dhf'
        block['gaunt'] = False if gaunt is None else gaunt
        block['breit'] = False if breit is None else breit
    return block

def casscf_block(trel, ncas, ncore, topt=True, thresh=1e-10, maxiter=100,
        robust=False, canonical=False, external_rdm=None):
    block = {
        'thresh': thresh, 
        'maxiter': maxiter,
        'nact': ncas,
        'nclosed': ncore,
        'canonical': canonical,
        'robust': robust
    }
    if trel:
        block['title'] = 'zcasscf'
        block['state'] = [1]
    else:
        block['title'] = 'casscf'
        block['nstate'] = [1]
    if not topt:
        block['algorithm'] = 'noopt'
    if external_rdm:
        assert isinstance(external_rdm, str)
        block['external_rdm'] = external_rdm
    return block

def caspt2_block(trel, ovlp=1e-9, external_rdm=None):
    block = {
        'method': 'caspt2',
        'thresh_overlap': ovlp,
        'frozen': False
    }
    block['title'] = 'relsmith' if trel else 'smith'
    if external_rdm:
        assert isinstance(external_rdm, str)
        block['external_rdm'] = external_rdm
    return block

def fci_block(trel, ncas, ncore, only_ints=False):
    block = {
        "only_ints": only_ints,
        "ncore": ncore,
        "norb":  ncas,
        "frozen": False,
        "state": [1]
    }
    if trel:
        block['title'] = 'zfci'
        block['state'] = [1]
    else:
        block['title'] = 'fci'
        block['nstate'] = 1
    return block

rdm_map = {
    '1RDM.1': 'fciqmc_0_0.rdm1',
    '2RDM.1': 'fciqmc_0_0.rdm2',
    '3RDM.1': 'fciqmc_0_0.rdm3',
    'CASPT2_AUX.1': 'fciqmc_0_0.rdm4f',
    'spinfree_OneRDM.1': 'fciqmc_0_0.rdm1',
    'spinfree_TwoRDM.1': 'fciqmc_0_0.rdm2',
    'spinfree_ThreeRDM.1': 'fciqmc_0_0.rdm3',
    'spinfree_CASPT2_AUX.1': 'fciqmc_0_0.rdm4f'
}
def bagel(config, dirname, bagel_dirname=None, neci_dirname=None, exe='/users/k1507071/code/rja_bagel/obj/src/BAGEL', 
        numthreads=1, mpiranks=1, destructive=False):
    if not os.path.exists(dirname): os.makedirs(dirname)
    relocator = check_hardlink if destructive else check_copyfile
    infile = os.path.abspath('{}/bagel.json'.format(dirname))
    outfile = os.path.abspath('{}/bagel.out'.format(dirname))
    assert not len(os.listdir(dirname)), 'Must have a clean directory to execute BAGEL'
    with open(infile, 'w') as f: json.dump(config, f, indent=4)
    if bagel_dirname is not None:
        # restarting from earlier invokation
        for f in ('ref.archive', 'relref.archive'):
            relocator('{}/{}'.format(bagel_dirname, f), '{}/{}'.format(dirname, f))
    if neci_dirname is not None:
        for fsrc, fdst in rdm_map.items():
            relocator('{}/{}'.format(neci_dirname, fsrc), '{}/{}'.format(dirname, fdst))

    print('Invoking BAGEL...')
    stdout, stderr, time = shell(
        'cd {}'.format(dirname),
        'export OMP_NUM_THREADS={}'.format(numthreads),
        'module purge',
        'module load /users/k1507071/opt/modules/gcc/7.4.0',
        'module load /users/k1507071/opt/modules/openmpi/4.0.1/gnu_7.4.0',
        'module load /users/k1507071/opt/modules/boost/1.70/gcc7.4.0/python2.7',
        '. $HOME/intel/mkl/bin/mklvars.sh intel64',
        'mpirun -np {} {} {} > {}'.format(mpiranks, exe, infile, outfile)
    )
    # if this BAGEL instance did not emit an FCIDUMP, link to previous  
    if not os.path.exists('{}/FCIDUMP'.format(dirname)):
        check_hardlink('{}/FCIDUMP'.format(bagel_dirname), '{}/FCIDUMP'.format(dirname))
    print('BAGEL instance completed ({:.2f} seconds)'.format(time))
    return time

#'/users/k1507071/code/neci_merge/build_ib/bin'
def neci(config, dirname, bagel_dirname=None, neci_dirname=None, neci_bin_dirname='/users/k1507071/lustre/code/neci_link/bin',
        mpiranks=1, destructive=False, quiet=False):
    if not os.path.exists(dirname): os.makedirs(dirname)
    relocator = check_hardlink if destructive else check_copyfile
    infile = os.path.abspath('{}/neci.inp'.format(dirname))
    outfile = os.path.abspath('{}/neci.out'.format(dirname))
    errfile = os.path.abspath('{}/neci.err'.format(dirname))
    jsonfile = os.path.abspath('{}/neci.json'.format(dirname))
    fcidumpfile = os.path.abspath('{}/FCIDUMP'.format(dirname))
    assert not len(os.listdir(dirname)), 'Must have a clean directory to execute NECI'
    with open(jsonfile, 'w') as f: json.dump(config, f, indent=4)
    for fname in ('FCIDUMP', 'FOCKMAT'):
        relocator('{}/{}'.format(bagel_dirname, fname), '{}/{}'.format(dirname, fname))
    if neci_dirname is not None:
        for fname in ['POPSFILEHEAD', 'POPSFILEBIN', '2RDM_POPSFILE', '3RDM_POPSFILE', 'CASPT2_AUX_POPSFILE']:
            relocator('{}/{}'.format(neci_dirname, fname), '{}/{}'.format(dirname, fname))
    norb  = get_fcidump_metadata('NORB', fname=fcidumpfile)
    nelec = get_fcidump_metadata('NELEC', fname=fcidumpfile)
    trel  = get_fcidump_metadata('TREL', bool, fname=fcidumpfile)
    config['nelec'] = nelec
    NeciInputFile(config).render(infile)
    print('Invoking NECI...')
    stdout, stderr, time = shell(
        'cd {}'.format(dirname),
        'module purge',
        'module load /users/k1507071/opt/modules/gcc/7.4.0',
        'module load /users/k1507071/opt/modules/openmpi/4.0.1/gnu_7.4.0',
        'mpirun -np {} {}/{} {} > {} 2> {}'.format(
            mpiranks, neci_bin_dirname, 'kdneci' if trel else 'dneci', infile, outfile, errfile)
    )
    energy = file_extract(outfile, 'REDUCED DENSITY MATRICES', -1)
    if not quiet:
        print('NECI instance completed ({:.2f} seconds)'.format(time))
        print('2RDM energy {:.8f}'.format(energy))
    return energy, time

def exact_caspt2_config(molecule, trel, ncas, ncore, gaunt=None, breit=None):
    return {
        'bagel':[
            molecule,
            hf_block(trel, thresh=1e-10, maxiter=100, robust=False, gaunt=gaunt, breit=breit),
            casscf_block(trel, ncas, ncore),
            caspt2_block(trel)
        ]
    }

def init_casscf_config(molecule, trel, ncas, ncore, gaunt=None, breit=None):
    return {
        'bagel':[
            molecule,
            hf_block(trel, thresh=1e-10, maxiter=100, robust=False, gaunt=gaunt, breit=breit),
            casscf_block(trel, ncas, ncore, topt=1, maxiter=0, external_rdm='noref'),
            fci_block(trel, ncas, ncore, only_ints=True)
        ]
    }

def init_to_iter_casscf_config(config):
    config = deepcopy(config)
    # remove HF block
    del config['bagel'][1]
    config['bagel'][1]['external_rdm'] = 'fciqmc'
    config['bagel'][1]['maxiter'] = 1
    return config

def init_to_iter_casscf_config(config):
    config = deepcopy(config)
    # remove HF block
    del config['bagel'][1]
    config['bagel'][1]['external_rdm'] = 'fciqmc'
    config['bagel'][1]['maxiter'] = 1
    return config

def iter_to_pseudocanonical_config(config):
    config = deepcopy(config)
    config['bagel'][1]['algorithm'] = 'noopt'
    config['bagel'][1]['canonical'] = True
    return config

def pseudocanonical_to_fockdump_config(config):
    config = deepcopy(config)
    config['bagel'][1]['canonical'] = False
    # remove integral output block
    del config['bagel'][-1]
    trel = config['bagel'][1]['title']=='zcasscf'
    config['bagel'].append(caspt2_block(trel, external_rdm='noref'))
    return config

def fockdump_to_caspt2_config(config):
    config = deepcopy(config)
    config['bagel'][-1]['external_rdm'] = 'fciqmc'
    return config

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
bagel(exact_caspt2_config(molecule, trel, ncas, ncore), 'exact')

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

def last_successful_casscf_iter():
    return None

def auto_init(molecule, trel, ncas, ncore):
    print('================== CASSCF iteration {:<3}=================='.format(0))
    init_bagel_config = init_casscf_config(molecule, trel, ncas, ncore)
    bagel(init_bagel_config, 'bagel_0')
    neci(neci_data, 'neci_0', 'bagel_0', mpiranks=neci_mpiranks)
    iter_bagel_config = init_to_iter_casscf_config(init_bagel_config)

def auto_iter(prev_iter=None):
    if prev_iter is None: prev_iter = last_successful_casscf_iter()
    print('================== CASSCF iteration {:<3}=================='.format(i))
    bagel(iter_bagel_config, 'bagel_'+str(prev_iter+1), 'bagel_'+str(prev_iter), 'neci_'+str(prev_iter))
    neci(neci_data, 'neci_'+str(prev_iter+1), 'bagel_'+str(prev_iter+1), mpiranks=neci_mpiranks)

def auto_fockmat(prev_iter=None):
    if prev_iter is None: prev_iter = last_successful_casscf_iter()
    print('=============== PSEUDOCANONICAL iteration ===============')
    pseudocanonical_bagel_config = iter_to_pseudocanonical_config(iter_bagel_config)
    bagel(pseudocanonical_bagel_config, 'bagel_pseudocanonical', 'bagel_'+str(prev_iter), 'neci_'+str(prev_iter))
    neci(neci_data, 'neci_pseudocanonical', 'bagel_pseudocanonical', mpiranks=neci_mpiranks)

    print('=============== FOCK MATRIX iteration ===============')
    fockdump_bagel_config = pseudocanonical_to_fockdump_config(pseudocanonical_bagel_config)
    bagel(fockdump_bagel_config, 'bagel_fockdump', 'bagel_pseudocanonical', 'neci_pseudocanonical')
    reformat_fockmat('bagel_fockdump/FOCKMAT')

def auto_caspt2(ovlps=(4, 5, 6, 7, 8, 9)):
    print('================= CASPT2 iteration =================')
    neci_data['rdm']['samplingiters'] = 500
    neci_data['rdm']['mrpt'] = {
        'type': 'caspt2',
        'granularity': 1,
        'promotionfractions': 1
    }
    neci(neci_data, 'neci_caspt2', 'bagel_fockdump', mpiranks=neci_mpiranks)
    caspt2_bagel_config = fockdump_to_caspt2_config(fockdump_bagel_config)
    for ovlp in ovlps:
        caspt2_bagel_config['bagel'][-1]['thresh_overlap'] = 10**-ovlp
        bagel(caspt2_bagel_config, 'bagel_caspt2/thresh_{}'.format(ovlp), 'bagel_pseudocanonical', 'neci_caspt2')

auto_init(molecule, trel, ncas, ncore)
for i in range(niter_casscf):
    auto_iter(i)
auto_fockmat(i)
auto_caspt2()

