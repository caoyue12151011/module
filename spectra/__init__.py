"""
Model initialization. Create some variables for latter use. 

Returns
-------
Molatalog: dict of parameters of molecules (linear or symmetric top only)
    All basic parameters, no derived parameters.

    name: e.g. 'NH3', 'H13CO+'

        'sw_ls': True for linear, False for symmetric top
        'shape': [symmetric top only] 'oblate' or 'prolate'
        'A0': [symmetric top only] Rigid rotor constants
        'B0': Rigid rotor constants
        'C0': [symmetric top only] Rigid rotor constants
        'mu': Electric dipole moment
        'sigma': [symmetric top only] Symmetry number

Splatalog: dict of parameters of transitions

    name: e.g. 'NH3_1-1', 'H13CO+_1-0'

        # Basic parameters
        'J': nominal J quantum number 
        'dJ': change of J 
        'K': nominal K quantum number 
        'E_u': Upper energy level
        'S': Line strength
        'nu0': Rest frequency
        'R': 1darray, Relative strength of hyperfine lines, sum(R)=1
        'v': 1darray, Velocity shifts of hyperfine lines w.r.t. 'nu0', 
            radio velocity convention used. I.e. v = (nu0-nu)/nu0 * c,
            while optical v = (lambda-lambda0)/lambda0 * c 

        # derived parameters
        'J_u': upper J quantum number
        'g_u': Degeneracy of the upper energy level
        'C_I', 'C_N', 'C_Q1', 'C_T1', 'C_T2':
                dimensionless constants used for calculation
        'C_Q2', 'C_Q3', 'C_Q4': [symmetric top only]
"""
import dill
import socket
import numpy as np
import astropy.units as u

# path of this package
hostname = socket.gethostname()
path = None
if hostname == 'Yues-MBP':
    path = '/Users/yuecao/Documents/coding/module/spectra'
elif hostname == 'yue-caos-ubuntu':
    path = '/home/dev/Documents/coding/module/spectra'

# Molatalog -------------------------------------------------------------------

# load data 
f = open(f'{path}/data/Molatalog.tab', 'r')
head = f.readline().split()
units = f.readline().split()
f.readline()  # skip the separator
data = f.readlines()
f.close()

Molatalog = {}
for line in data:
    line = line.split()
    name = line[0]
    Molatalog[name] = {}

    for i in range(1, len(head)):  # skip the first item in head
        key = head[i]

        unit = units[i]
        if '[' in unit:
            unit = u.Unit(unit[1:-1])  # remove '[]'
        else: 
            unit = 1
        
        value = line[i]
        if value == '-':
            continue  # skip empty values

        try:
            value = float(value) * unit  # handle floats and strings
        except:
            if value == 'T':   # handle bools
                value = True
            elif value == 'F':
                value = False

        Molatalog[name][key] = value

# Splatalog -------------------------------------------------------------------

# load data 
f = open(f'{path}/data/Splatalog.tab', 'r')
head = f.readline().split()
units = f.readline().split()
f.readline()  # skip the separator
data = f.readlines()
f.close()

Splatalog = {}
for line in data:
    line = line.split()
    name = line[0]
    Splatalog[name] = {}

    for i in range(1, len(head)):  # skip the first item in head
        key = head[i]

        unit = units[i]
        if '[' in unit:
            unit = u.Unit(unit[1:-1])  # remove '[]'
        else: 
            unit = 1
        
        value = line[i]
        if value == '-':
            continue  # skip empty values

        try:
            value = float(value) * unit  # handle floats and strings
        except:
            if value == 'T':   # handle bools
                value = True
            elif value == 'F':
                value = False

        Splatalog[name][key] = value


# load supplement data of hyperfine-line info
f = open(f'{path}/data/Splatalog_hyperfine.tab', 'r')
head = f.readline().split()
units = f.readline().split()
f.readline()  # skip the separator
data = f.readlines()
f.close()

for line in data:
    line = line.split()
    name = line[0]

    for i in range(1, len(head)):  # skip the first item in head
        key = head[i]

        unit = units[i]
        if '[' in unit:
            unit = u.Unit(unit[1:-1])  # remove '[]'
        else: 
            unit = 1
        
        value = line[i]
        if value == '-':
            continue  # skip empty values
        value = np.array([float(item) for item in value.split(',')]) * unit

        Splatalog[name][key] = value

# derived keys ................................................................
# 'J_u', 'g_u'
# 'C_I', 'C_N', 'C_Q1', 'C_Q2', 'C_Q3', 'C_Q4', 'C_T1', 'C_T2'

for line in Splatalog:
    mole = line.split('_')[0]
    sw_ls = Molatalog[mole]['sw_ls']

    tmp = Splatalog[line]
    J = tmp['J']
    dJ = tmp['dJ']
    K = tmp['K']
    E_u = tmp['E_u']
    S = tmp['S']
    nu0 = tmp['nu0']
    R = tmp['R']
    v = tmp['v']

    # J_u
    J_u = max(J, J+dJ)

    # g_u
    g_Ju = 2*J_u + 1
    g_Ku = 2 if (not sw_ls) and K != 0 else 1
    g_Iu = 1
    g_u = g_Ju * g_Ku * g_Iu



# saves
dill.dump(Molatalog, open(f'{path}/data/Molatalog.p','wb'))
dill.dump(Splatalog, open(f'{path}/data/Splatalog.p','wb'))


# Splatalog0 = dill.load(open('data/Splatalog0.p','rb'))
# for line in sorted(Splatalog):
#     R = Splatalog[line]['R']
#     v = Splatalog[line]['v']
#     R = ','.join([str(i) for i in R])
#     v = ','.join([str(i) for i in v])
#     print(f'{line}\t{R}\t{v}')