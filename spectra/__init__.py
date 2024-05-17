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
        'E_u': Upper energy level
        'J': nominal J quantum number 
        'J_u': upper J quantum number
        'K': nominal K quantum number 
        'R': 1darray, Relative strength of hyperfine lines, sum(R)=1
        'S': Line strength
        'dJ': change of J 
        'g_u': Degeneracy of the upper energy level
        'nu0': Rest frequency
        'v': 1darray, Velocity shifts of hyperfine lines w.r.t. 'nu0', 
            radio velocity convention used. I.e. 

"""







keys = ['A0', 'B0', 'C0', 'C_I', 'C_N', 'C_Q1', 'C_Q2', 'C_Q3', 'C_Q4', 
        'C_T1', 'C_T2', 'E_u', 'J', 'J_u', 'K', 'R', 'S', 'dJ', 'g_u', 
        'mu', 'nu0', 'shape', 'sigma', 'sw_ls', 'v']

for line in sorted(Splatalog):
    if 'sw_ls' in Splatalog[line]:
        print(line, Splatalog[line]['sw_ls'])


Basic 

Mol: 'A0', 'B0', 'C0',  
        'mu', 'shape', 'sigma', 'sw_ls', 

Sp: 'E_u', 'J', 'J_u', 'K', 'R', 'S', 'dJ', 'g_u', 
        'nu0', 'v'


Derived 

Sp: 'C_I', 'C_N', 'C_Q1', 'C_Q2', 'C_Q3', 'C_Q4', 
        'C_T1', 'C_T2'