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
    try:
        return type_(str(stdout, 'utf-8').split()[split_offset])
    except IndexError:
        return float('nan')

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

class BatchFile(StructuredFile):
    def body(self):
        self.extend([
        '#!/bin/bash -l',
        '#SBATCH --nodes={}'.format(self['nodes']),
        '#SBATCH --partition="morty,nms_research"',
        '#SBATCH --ntasks-per-node={}'.format(self['taskspernode']),
        '#SBATCH --time={}'.format(self['walltime']),
        '#SBATCH --job-name={}'.format(self['jobname'])
        ])
        self.extend(self['commandlines'])

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
        'nonuniformrandexcits {}'.format('4ind-weighted' if not self['trel'] else 'pick-virt-uniform-mag')
        ])
        try:
            if self['spinrestrict']:
                self << 'spin-restrict {}'.format(self['spinrestrict'])
        except KeyError:
            self << 'spin-restrict 0'
        self.extend([
        'endsys',
        'calc',
        'methods',
        'method vertex fcimc',
        'endmethods',
        'totalwalkers {}'.format(self['nw']),
        'memoryfacpart {}'.format(self['memoryfacpart']),
        'memoryfacspawn {}'.format(self['memoryfacspawn']),
        'seed {}'.format(self['seed']),
        'startsinglepart {}'.format(self['startsinglepart']),
        'shiftdamp {}'.format(self['shiftdamp']),
        'diagshift {}'.format(self['diagshift']),
        'truncinitiator',
        'addtoinitiator {}'.format(self['addtoinitiator']),
        'allrealcoeff',
        'realspawncutoff 0.4',
        'stepsshift {}'.format(self['stepsshift']),
        'maxwalkerbloom 3',
        'load-balance-blocks off',
        'nmcyc {}'.format(self['nmcyc'])
        ])
        try:
            if self['walkcontgrow']: self << 'walkcontgrow'
        except KeyError: pass
        try:
            if self['jumpshift']: self << 'jump-shift'
        except KeyError: pass
        if self['tau']:
            self << 'tau {}'.format(self['tau'])
        if self['readpops']:
            self << 'readpops'
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
            #try:
            #    self << 'max-rank-hbrdm-semi-stoch-fill {}'.format(self['rdm']['maxrankhbrdmsemistochfill'])
            #except KeyError: pass

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
    'diagshift': 0.1,
    'orthogonalisereplicas': None,
    'addtoinitiator': 3.0,
    'stepsshift': 5,
    'nmcyc': -1,
    'readpops': False,
    'jumpshift': False,
    'tau': None,
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
        'maxrankhbrdmsemistochfill':4,
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

def hf_block(trel, thresh=1e-8, maxiter=100, charge=0, robust=False, gaunt=None, breit=None):
    block = {
        'thresh':thresh,
        'maxiter':maxiter,
        'charge':charge,
        'robust':robust,
        'pop':True
    }
    if not trel:
        block['title'] = 'hf'
    else:
        block['title'] = 'dhf'
        block['gaunt'] = False if gaunt is None else gaunt
        block['breit'] = False if breit is None else breit
    return block

def casscf_block(trel, ncas, ncore, topt=True, thresh=1e-10, maxiter=100,
        robust=False, canonical=False, external_rdm=None, charge=None):
    block = {
        'thresh': thresh, 
        'maxiter': maxiter,
        'nact': ncas,
        'fci_algorithm': 'hz',
        'nclosed': ncore,
        'canonical': canonical,
        'robust': robust
    }
    if charge is not None:
        block['charge'] = charge
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
def bagel(config, dirname, bagel_dirname=None, neci_dirname=None, exe='/users/k1507071/lustre/code/rja_bagel/obj/src/BAGEL', 
        numthreads=1, mpiranks=1, destructive=False, quiet=False):
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
    gradient = file_extract(outfile, 'Using the second-order algorithm', 3, 2)
    if not quiet:
        print('BAGEL instance completed ({:.2f} seconds)'.format(time))
        print('CASSCF Gradient {}'.format(gradient))
    return gradient, time

#'/users/k1507071/code/neci_merge/build_ib/bin'
nbd='/users/k1507071/lustre/code/neci_merge/build_ib/bin'
#nbd='/users/k1507071/lustre/code/neci_master/build/bin'
def neci(config, dirname, bagel_dirname=None, neci_dirname=None, neci_bin_dirname=nbd,
        mpiranks=1, destructive=False, quiet=False, batch_opts=None, openmpiv4=True):
    if not os.path.exists(dirname): os.makedirs(dirname)
    relocator = check_hardlink if destructive else check_copyfile
    submitfile = os.path.abspath('{}/submit.sh'.format(dirname))
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
    config['trel'] = trel
    NeciInputFile(config).render(infile)
    if config['rdm'] is None and config['orthogonalisereplicas'] is None:
        bin_name = 'kneci' if trel else 'neci'
    else:
        bin_name = 'kdneci' if trel else 'dneci'
    if openmpiv4: 
        mpimod = 'libs/openmpi/4.0.0-gcc7.3.0'
    else:
        mpimod = 'libs/openmpi/3.1.4/gcc/7.3.0'
    print('Invoking NECI...')
    command_lines = [
        'cd {}'.format(dirname),
        'module purge',
        'module load {} libs/lapack/3.8.0/gcc-7.3.0'.format(mpimod),
        'module li',
        'which mpirun',
        '. $HOME/intel/mkl/bin/mklvars.sh intel64',
        'mpirun -np {} {}/{} {} > {} 2> {}'.format(
            mpiranks, neci_bin_dirname, bin_name, infile, outfile, errfile)
    ]

    if batch_opts is None:
        stdout, stderr, time = shell(*command_lines)
        energy = file_extract(outfile, 'REDUCED DENSITY MATRICES', -1)
        if not quiet:
            print('NECI instance completed ({:.2f} seconds)'.format(time))
            print('2RDM energy {:.8f}'.format(energy))
    else:
        energy = None
        tmp = {'commandlines':command_lines}
        tmp.update(batch_opts)
        BatchFile(tmp).render(submitfile)
        stdout, stderr, time = shell('cd {}'.format(dirname), 'sbatch {}'.format(submitfile))
        jobid = int(str(stdout, 'utf8').strip().split()[-1])
        print('NECI batch job id: {}'.format(jobid))
    return energy, time

def exact_caspt2_config(molecule, trel, ncas, ncore, charge=0, gaunt=None, breit=None):
    return {
        'bagel':[
            molecule,
            hf_block(trel, maxiter=100, charge=charge, robust=False, gaunt=gaunt, breit=breit),
            casscf_block(trel, ncas, ncore, charge=charge),
            caspt2_block(trel)
        ]
    }

def init_casscf_config(molecule, trel, ncas, ncore, charge=0, gaunt=None, breit=None):
    return {
        'bagel':[
            molecule,
            hf_block(trel, maxiter=100, robust=False, charge=charge, gaunt=gaunt, breit=breit),
            casscf_block(trel, ncas, ncore, topt=1, maxiter=0, external_rdm='noref', charge=charge),
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


def last_successful_casscf_iter():
    i = 0
    while 1:
        if not os.path.exists('neci_{}'.format(i)):
            return i-1
        i+=1

def get_bagel_config_from_disk(it):
    with open('bagel_{}/bagel.json'.format(it), 'r') as f: return json.load(f)

def get_neci_data_from_disk(it):
    with open('neci_{}/neci.json'.format(it), 'r') as f: return json.load(f)

def get_configs_from_disk(it):
    return get_bagel_config_from_disk(it), get_neci_data_from_disk(it)

def auto_init(molecule, neci_data, neci_mpiranks, bagel_numthreads, bagel_mpiranks, trel, ncas, ncore, charge):
    print('================== CASSCF iteration {:<3}=================='.format(0))
    init_bagel_config = init_casscf_config(molecule, trel, ncas, ncore, charge=charge)
    bagel(init_bagel_config, 'bagel_0', numthreads=bagel_numthreads, mpiranks=bagel_mpiranks)
    neci(neci_data, 'neci_0', 'bagel_0', mpiranks=neci_mpiranks)

def auto_iter(neci_mpiranks, bagel_numthreads, bagel_mpiranks, prev_iter=None, neci_data=None, pops_restart=False):
    if prev_iter is None: prev_iter = last_successful_casscf_iter()
    if neci_data is None:
        bagel_config, neci_data = get_configs_from_disk(prev_iter)
    else:
        pops_restart = neci_data['readpops']
        bagel_config = get_bagel_config_from_disk(prev_iter)
    if prev_iter==0: bagel_config = init_to_iter_casscf_config(bagel_config)
    print('================== CASSCF iteration {:<3}=================='.format(prev_iter+1))
    grad, _ = bagel(bagel_config, 'bagel_'+str(prev_iter+1), 'bagel_'+str(prev_iter), 'neci_'+str(prev_iter), 
            numthreads=bagel_numthreads, mpiranks=bagel_mpiranks)
    if not pops_restart:
        neci(neci_data, 'neci_'+str(prev_iter+1), 'bagel_'+str(prev_iter+1), mpiranks=neci_mpiranks)
    else:
        neci(neci_data, 'neci_'+str(prev_iter+1), 'bagel_'+str(prev_iter+1), mpiranks=neci_mpiranks, 
                neci_dirname='neci_'+str(prev_iter))
    return grad

def auto_casscf(neci_mpiranks, bagel_numthreads, bagel_mpiranks, conv_grad, niter_casscf=1, neci_data=None):
    for i in range(niter_casscf):
        grad = auto_iter(neci_mpiranks, bagel_numthreads, bagel_mpiranks, neci_data=neci_data)
        if abs(grad)<conv_grad:
            print('CASSCF CONVERGED WITHIN {}'.format(conv_grad))
            return

def auto_fockmat(neci_mpiranks, bagel_numthreads, bagel_mpiranks, prev_iter=None, neci_data=None, pops_restart=False):
    if prev_iter is None: prev_iter = last_successful_casscf_iter()
    if neci_data is None:
        bagel_config, neci_data = get_configs_from_disk(prev_iter)
    else:
        pops_restart = neci_data['readpops']
        bagel_config = get_bagel_config_from_disk(prev_iter)
    print('=============== PSEUDOCANONICAL iteration ===============')
    bagel_config = iter_to_pseudocanonical_config(bagel_config)
    bagel(bagel_config, 'bagel_pseudocanonical', 'bagel_'+str(prev_iter), 'neci_'+str(prev_iter),
            numthreads=bagel_numthreads, mpiranks=bagel_mpiranks)
    if not pops_restart:
        neci(neci_data, 'neci_pseudocanonical', 'bagel_pseudocanonical', mpiranks=neci_mpiranks)
    else:
        neci(neci_data, 'neci_pseudocanonical', 'bagel_pseudocanonical', mpiranks=neci_mpiranks,
                neci_dirname='neci_'+str(prev_iter))

    print('=============== FOCK MATRIX iteration ===============')
    bagel_config = pseudocanonical_to_fockdump_config(bagel_config)
    bagel(bagel_config, 'bagel_fockdump', 'bagel_pseudocanonical', 'neci_pseudocanonical',
            numthreads=bagel_numthreads, mpiranks=bagel_mpiranks)
    reformat_fockmat('bagel_fockdump/FOCKMAT')

def auto_caspt2(neci_mpiranks, bagel_numthreads, bagel_mpiranks, ovlps=(4, 5, 6, 7, 8, 9), neci_data=None, pops_restart=False):
    if neci_data is None: neci_data= get_neci_data_from_disk('pseudocanonical')
    else: pops_restart = neci_data['readpops']
    print('================= CASPT2 iteration =================')
    neci_data['rdm']['samplingiters'] = 500
    neci_data['rdm']['mrpt'] = {
        'type': 'caspt2',
        'granularity': 1,
        'promotionfractions': 1
    }
    if not pops_restart:
        neci(neci_data, 'neci_caspt2', 'bagel_fockdump', mpiranks=neci_mpiranks)
    else:
        neci(neci_data, 'neci_caspt2', 'bagel_fockdump', neci_dirname='neci_fockdump')
    bagel_config = get_bagel_config_from_disk('fockdump')
    bagel_config = fockdump_to_caspt2_config(bagel_config)
    for ovlp in ovlps:
        bagel_config['bagel'][-1]['thresh_overlap'] = 10**-ovlp
        bagel(bagel_config, 'bagel_caspt2/thresh_{}'.format(ovlp), 'bagel_pseudocanonical', 'neci_caspt2',
            numthreads=bagel_numthreads, mpiranks=bagel_mpiranks)

