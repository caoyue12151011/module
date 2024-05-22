"""
Model initialization. Create some variables for latter use. 

Returns
-------
Molatalog: dict of parameters of molecules (linear or symmetric top only)
    All basic parameters, no derived parameters.

    name: e.g. 'NH3', 'H13CO+'

        'shape': 'linear', 'oblate' or 'prolate', latter two are for 
            symmetric top molecules
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
        'nu0': Rest frequency
        'R': 1darray, Relative strength of hyperfine lines, sum(R)=1
        'v': 1darray, Velocity shifts of hyperfine lines w.r.t. 'nu0', 
            radio velocity convention used. I.e. v = (nu0-nu)/nu0 * c,
            while optical v = (lambda-lambda0)/lambda0 * c 

        # derived parameters
        'J_u': upper J quantum number
        'g_u': Degeneracy of the upper energy level
        'S': Line strength
        'C_I', 'C_N' 'C_T1', 'C_T2': dimensionless constants used for 
                                     calculation
        'C_Q1': [linear only]
        'C_Q2', 'C_Q3', 'C_Q4': [symmetric top only]
"""
import os
import dill
import numpy as np
import astropy.units as u
import astropy.constants as ct

# path of this package
path = os.path.dirname(os.path.abspath(__file__))

# load original splatalog
Splatalog0 = dill.load(open(f'{path}/data/Splatalog0.p','rb'))

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

    # molecular constants
    shape = Molatalog[mole]['shape']
    mu = Molatalog[mole]['mu']
    B0 = Molatalog[mole]['B0']
    if shape != 'linear':
        A0 = Molatalog[mole]['A0']
        C0 = Molatalog[mole]['C0']
        sigma = Molatalog[mole]['sigma']

    # basic transitional constants
    tmp = Splatalog[line]
    J = tmp['J']
    dJ = tmp['dJ']
    K = tmp['K']
    E_u = tmp['E_u']
    nu0 = tmp['nu0']
    R = tmp['R']
    v = tmp['v']

    # J_u
    J_u = max(J, J+dJ)

    # g_u
    g_Ju = 2*J_u + 1
    g_Ku = 2 if (shape != 'linear') and K != 0 else 1
    g_Iu = 1
    g_u = g_Ju * g_Ku * g_Iu

    # S
    if dJ == -1:
        S = (J**2-K**2)/J/(2*J+1)
    elif dJ == 0:
        S = K**2/J/(J+1)
    elif dJ == 1:
        S = ((J+1)**2-K**2)/(J+1)/(2*J+1)

    # C_I
    C_I = (2*ct.h* nu0**3 / ct.c**2).to_value(u.Jy)

    # C_N
    C_N = 4*2**.5 * np.pi**2.5 * S * mu**2 * R * g_u / (3*ct.h)
    C_N = C_N.to_value(u.km*u.cm**2/u.s)

    # C_T
    C_T1 = (E_u/ct.k_B).to_value(u.K)
    C_T2 = (ct.h*nu0/ct.k_B).to_value(u.K)

    # C_Q1-4
    if shape == 'linear':
        C_Q1 = (ct.h*B0/ct.k_B).to_value(u.K)
    else:
        m = B0/A0 if shape=='prolate' else B0/C0
        C_Q2 = ((m*np.pi)**.5/sigma*(ct.k_B/ct.h/B0)**1.5).to_value(u.K**-1.5)
        C_Q3 = (ct.h*B0*(4-m)/(12*ct.k_B)).to_value(u.K)
        C_Q4 = ((ct.h*B0*(1-m)/ct.k_B)**2/90).to_value(u.K**2)

    # save
    Splatalog[line]['J_u'] = J_u
    Splatalog[line]['g_u'] = g_u
    Splatalog[line]['S'] = S
    Splatalog[line]['C_I'] = C_I
    Splatalog[line]['C_N'] = C_N
    Splatalog[line]['C_T1'] = C_T1
    Splatalog[line]['C_T2'] = C_T2
    if shape == 'linear':
        Splatalog[line]['C_Q1'] = C_Q1
    else:
        Splatalog[line]['C_Q2'] = C_Q2
        Splatalog[line]['C_Q3'] = C_Q3
        Splatalog[line]['C_Q4'] = C_Q4

# self-check ..................................................................
"""
Except for 'C18O_1-0', 'CO_2-1', which have changed B0 values, all constants 
have differences < 1e-6 but > 1e-7.
"""
# parameters 
tol = 1e-6
keys = ['J_u','g_u','S','C_I','C_N','C_T','C_Q1','C_Q2','C_Q3','C_Q4']

for line in sorted(Splatalog):
    tmp = Splatalog[line]
    tmp0 = Splatalog0[line]

    for key in keys:
        if key in tmp:
            value = np.array(tmp[key])
            value0 = np.array(tmp0[key])
            diff = np.max(np.abs(value-value0)/value0)
            if diff > tol:
                print(f'{line} {key}: diff > {tol:.1e}.')

# saves -----------------------------------------------------------------------

# dills
dill.dump(Molatalog, open(f'{path}/data/Molatalog.p','wb'))
dill.dump(Splatalog, open(f'{path}/data/Splatalog.p','wb'))