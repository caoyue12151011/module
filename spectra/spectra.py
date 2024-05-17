''' 
To analyze the spectral data with the rotational transition models. Can handle
transitions of linear or symmetric top molecules only.

Functions
---------
gaussian: calculate the Gaussian function

find_peaks: find +/- peaks in spectra and estimate the parameters.

calc_Q: calculate the rotational partition function

specmod: model spectrum with {Tex, N, v, sigma_v} as parameters

fitter: fit the spectral line with the given model

calc_tau: calculate mean tau given lgN, tex, sigma_v

calc_lgN: calculate lgN given tau, tex, sigma_v

calc_lgN_ch: calculate lgN in a velocity range (without the assumption of
    Gaussian velocity profile) given tau, tex, dv. Useful for e.g. deriving 
    CO column density given tau (from 12CO/13CO data) and 12CO (2-1) intensity.

calc_tau_iso_ratio: To calculate the opacity given the line ratio & the column 
    density ratio of two isotopic species

test: to do some tests 

How to use
----------
1. Determine the number of apparent line number and the index of the reference
    apparent line (in ascending-v order). See image/spectra/model to see the 
    modeled lines.
2. Run find_peaks. Adjust the input parameters for good results. Use its 
    outputs together with the parameters in step 1 to determine the initial 
    guess for fitting. The initial guess is referred to the reference
    apparent line. 
3. Run fitter to fit the spectral lines.
'''
import time
import dill
import socket
import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u 
import astropy.constants as ct 
from scipy.optimize import curve_fit, fsolve
from scipy.ndimage import gaussian_filter1d

# path of this package
hostname = socket.gethostname()
path = None
if hostname == 'Yues-MBP':
    path = '/Users/yuecao/Documents/coding/module/spectra'
elif hostname == 'yue-caos-ubuntu':
    path = '/home/dev/Documents/coding/module/spectra'

# load spectra data
Splatalog = dill.load(open(f'{path}/data/Splatalog.p', 'rb'))


def gaussian(x, *para): 
    '''
    To calculate the 1D multi-Gaussian function y = a*exp(-(x-b)^2/(2*c^2)).
    
    Inputs
    ------
    x: 1D array, independent variable
    para: 1D array, parameters [a1, b1, c1, a2, b2, c2, ...]
    
    Returns
    -------
    y: 1D array
    '''
    a, b, c = np.reshape(para,(-1,3)).transpose()[...,np.newaxis]  # (para, x)
    return np.sum(a*np.exp(-(x-b)**2/(2*c**2)), axis=0)


def find_peaks(x, y, sign, lim): 
    '''
    To find peaks as 1D Gaussian in the spectrum.

    Parameters
    ----------
    x: 1D array of equally spaced floats, the independent variable
    y: 1D array, the spectrum
    sign: {-1|0|1}, 1 for positive peaks, -1 for negative peaks, 0 for both
    lim: list or 1D array, boundaries of the parameters of Gaussian peaks,
        in the order of [amp_min, amp_max, pos_min, pos_max, sigma_min, 
        sigma_max]

    Returns
    --------
    para0: list, estimated Gaussian parameters, [amp, pos, sigma, amp, ...]

    Algorithm
    ----------
    1. Calculate dy/dx and smooth it according to sigma_min. Investigation 
       shows that peaks with sigmas < length/5 will be smooth out efficiently.
    2. Find peaks by identifying zero-across in dy/dx.
    3. Check the amplitudes of the peaks.
    '''
    para0 = []

    # constants
    th_snr_dy = 1 # threshold for determining the signs of dy
    smth_ratio = 1 # smoothing kernel sigma = smth_ratio*sigma_min

    # length of one pixel
    dx = abs(x[1]-x[0])

    # smooth the spectrum
    y = gaussian_filter1d(y,smth_ratio*lim[4]/dx)

    # derivative
    dy = np.append(np.diff(y)/np.diff(x),0)
    if 'c_min' in lim:
        dy = gaussian_filter1d(dy,smth_ratio*lim[4]/dx)
    std_dy = np.std(dy)

    mono = np.zeros_like(dy,dtype=np.int8)
    mono[dy > th_snr_dy*std_dy] = 1
    mono[dy < -th_snr_dy*std_dy] = -1
    
    # find peaks & check
    if sign==1 or sign==0:
        i = 0
        while i < len(x)-1:
            find = False
            if mono[i]==1:
                j = i 
                while j < len(x)-1:
                    if mono[j+1] > mono[j]:
                        break
                    j += 1
                if mono[j]==-1:
                    amp = y[int((i+j)/2)]
                    pos = x[int((i+j)/2)]
                    wid = (x[j]-x[i])/2.*.5
                    # print(amp,pos,wid,x[i],x[j])
                    if (lim[0] < amp < lim[1] and lim[2] < pos < lim[3] 
                        and lim[4] < wid < lim[5]):
                        para0.extend([amp,pos,wid])
                        find = True
                i = j       
            if not find:
                i += 1

    if sign==-1 or sign==0:
        i = 0
        while i < len(x)-1:
            find = False
            if mono[i]==-1:
                j = i
                while j < len(x)-1:
                    if mono[j] < mono[j-1]:
                        break
                    j += 1
                if mono[j-1]==1:
                    amp = y[int((i+j)/2)]
                    pos = x[int((i+j)/2)]
                    wid = (x[j]-x[i])/2.*.5
                    if (lim[0] < amp < lim[1] and lim[2] < pos < lim[3] 
                        and lim[4] < wid < lim[5]):
                        para0.extend([amp,pos,wid])
                        find = True
                i = j
            if not find:
                i += 1

    ''' display
    plt.figure()
    plt.plot(x,y)
    plt.figure()
    plt.plot(x,dy,x,dy.mean()+mono*(dy.max()-dy.mean()))
    plt.show()
    print(para0)
    #'''
    return para0


def calc_Q(model, tex, sw_ls):
    '''
    To calculate the rotational partition function Q_rot(Tex).

    Inputs
    ------
    model: string, name of the transition
    tex: [K], scalar or ndarray, excitation temperature
    sw_ls: T for linear, F for symmetric top molecules

    Outputs
    -------
    Q_rot: the rotational partition function
    '''
    tmp = Splatalog[model]
    Q_rot = None
    if sw_ls:
        C_Q1 = tmp['C_Q1']
        Q_rot = 1/C_Q1 * tex * np.exp(C_Q1/3/tex)
    else:
        C_Q2 = tmp['C_Q2']
        C_Q3 = tmp['C_Q3']
        C_Q4 = tmp['C_Q4']
        Q_rot = C_Q2 * tex**1.5 * np.exp(C_Q3/tex) * (1+C_Q4/tex**2)
    return Q_rot


def specmod(x, model, T_bg, Switch, para_type, *para): 
    '''
    To generate the modeled spectrum.

    Parameters
    ----------
    x: [km/s], 1D array, line-of-sight velocity in the LSR frame. The 
        independent variable for calculating the spectrum

    model: string, name of the transition, see spectra/Splatalog.p for 
        available models

    T_bg: [K], scalar, background brightness temperature

    Switch: dict of switches controlling the returned spectrum. Bool values.

        'thin': whether to use the optically thin limit
        'intensity': whether to calculate intensity or brightness temperature
        'hyperfine': whether to output the spectra of hyperfine lines, only 
            valid when thin=True
        'collapsed_hyperfine': whether to set the hyperfine velocity offsets to 
            0, used for demonstrating only one line, shouldn't use for 
            quantitative analyses

    para_type: str, type of the spectrum function. Can be 'T-N', 'A-T-N', 
        'T-tau'. See the description of *para for how it influences 'para'.
        'T-N' uses temperature & colomn density. 'A-T-N' adds an additional 
        amplitude parameter in case that 'T-N' becomes wrongly optically 
        thick and the spectrum's amplitude is limited. 'T-tau' replaces N with
        opacity.

    *para: 1D array, parameters of the model, Format varies with 'para_type'
        'T-N':   [tex, lgN1, v0_1, sigma_v1, lgN2, v0_2, sigma_v2, ...]
        'A-T-N': [amp1, tex, lgN1, v0_1, sigma_v1, lgN2, v0_2, sigma_v2, ...]
        'T-tau': [tex, tau1, v0_1, sigma_v1, tau2, v0_2, sigma_v2, ...]
        where tex is excitation temperature in K, amp=B_nu(tex) is in Jy/sr,
        lgNi, v0_i, sigma_vi are log10(column density) [cm^-2] of the TRACER, 
        bulk velocity [km/s] and 1-sigma velocity dispersion [km/s] of the ith 
        velocity component.

    Returns
    -------
    y: the modeled spectrum.
        If Switch['thin'] & Switch['hyperfine'], y is 2D array of shape 
        (n_hyper+1, len(x)), where y[0] is the sum of all the hyperfine 
        spectra. Otherwise, y is 1D array with the same length as x's.
        If Switch['intensity'] y is the intensity in Jy/sr, else y is the 
        brightness temperature in K.

    Notices
    -------
    See spectra/model.md and other docs for details of the models.
    '''
    # handle parameters -------------------------------------------------------
    para = np.array(para)

    # spectral parameters
    tmp = Splatalog[model]
    sw_ls = tmp['sw_ls']
    C_N = tmp['C_N']
    C_T1 = tmp['C_T1']
    C_T2 = tmp['C_T2']
    C_I = tmp['C_I']
    v = tmp['v']
    R = tmp['R']

    # whether collapse hyperfine lines
    if sw_chyper:
        v = np.zeros_like(v)

    # reshape (n_hyper, n_component, x)
    C_N = C_N[..., np.newaxis, np.newaxis]
    v = v[..., np.newaxis, np.newaxis] 
    R = R[..., np.newaxis, np.newaxis]

    # switches
    sw_thin = Switch['thin']
    sw_I = Switch['intensity']
    sw_hyper = Switch['hyperfine']
    sw_chyper = Switch['collapsed_hyperfine']

    # macroscopic parameters, shape = (n_component, x)
    if para_type in ['T-N', 'T-tau']:
        tex = para[0]
        lgN, v0, sigma_v = np.reshape(para[1:],
                                      (-1,3)).transpose()[...,np.newaxis]
    elif para_type in ['A-T-N']:
        amp, tex = para[:2] 
        lgN, v0, sigma_v = np.reshape(para[2:],
                                      (-1,3)).transpose()[...,np.newaxis]

    # calculate the spectrum --------------------------------------------------

    # opacity
    if para_type in ['T-N', 'A-T-N']:
        # Q_rot
        Q_rot = calc_Q(model, tex, sw_ls)
        # (n_hyper, n_component, x)
        TAU = (C_N * 10**lgN / sigma_v / Q_rot *
            np.exp(-(x-v0-v)**2/(2*sigma_v**2)-C_T1/tex) *
            (np.exp(C_T2/tex)-1)) 
        
    elif para_type in ['T-tau']:
        TAU = 2*(np.log(2)/np.pi)**.5*R*tau*np.exp(-(x-v0-v)**2/(2*sigma_v**2))

    # collapse opacity
    if sw_thin and sw_hyper:
        tau = np.sum(TAU, axis=1)  # (n_hyper, len(x))
        tau = np.concatenate((np.sum(tau,axis=0)[np.newaxis], tau))
            # (n_hyper+1, len(x))
    else:
        tau = np.sum(TAU, axis=(0,1))  # shape = x.shape

    # terms
    term_tau = tau if sw_thin else 1-np.exp(-tau)
    if para_type in ['T-N', 'T-tau']:
        term_T = 1/(np.exp(C_T2/tex)-1) - 1/(np.exp(C_T2/T_bg)-1)
    elif para_type in ['A-T-N']:
        term_T = 1/(np.exp(C_T2/T_bg)-1)
    
    # modeled spectra
    if para_type in ['T-N', 'T-tau']:
        if sw_I:
            y = C_I * term_T * term_tau  # intensity [Jy/sr]
        else:
            y = C_T2 / np.log(1 + 1/term_tau/term_T)  # T_b [K]   
    
    elif para_type in ['A-T-N']:
        if sw_I:
            y = (amp - C_I*term_T)*term_tau 
        else:
            y = C_T2 / np.log(1 + 1/term_tau/(amp/C_I-term_T)) 

    return y


def fitter(f,x,y,model,para0,bounds,T_bg,Switch): # ---------------
    '''
    To fit the spectrum with the model in specmod.

    Parameters
    ----------
    f: the model function
    x: [km/s], 1D array, los velocity of the spectrum
    y: 1D array, the spectrum
    model: string, name of the model
    para0: 1D array, initial guesses
    bounds: ([p1_min, p2_min, ...], [p1_max, p2_max, ...])
    T_bg, Switch: cf. spectra_rot.specmod

    Returns
    --------
    success: {True|False}, whether fitting is successful
    para: 1D array, optimized parameters
    perr: 1D array, uncertainties of parameters
    '''

    try:
        para, pcov = curve_fit(lambda x,*para: f(x,model,T_bg,Switch,*para),
            x,y,para0,bounds=bounds) 
        perr = np.sqrt(np.diag(pcov))
        return True, para, perr

    except:
        return False, np.full(len(para0),np.nan), np.full(len(para0),np.nan)

    

def calc_tau(tex,lgN,sigma_v,model): 
    '''
    To calculate the mean opacity of a Gaussian velocity component, see 
    spectra_rot.md for the detailed algorithm.
    
    Inputs
    ------
    tex: [K], scalar or nD array, excitation temperature
    lgN: scalar or nD array with the same shape as tex, log10(N_tracer [cm^-2])
    sigma_v: [km/s], scalar or nD array with the same shape as tex, 
        1-sigma velocity dispersion
    model: string, name of the transition, see Splatalog for possible values
    
    Returns
    -------
    tau: mean opacity
    '''

    tmp = Splatalog[model]

    sw_ls = tmp['sw_ls']
    C_N = tmp['C_N']
    C_T1 = tmp['C_T1']
    C_T2 = tmp['C_T2']
    
    # Q_rot
    Q_rot = calc_Q(model,tex,sw_ls)

    # tau
    tau = (.5*(np.pi/np.log(2))**.5* np.sum(C_N)* 10**lgN/sigma_v/Q_rot*
        np.exp(-C_T1/tex)* (np.exp(C_T2/tex)-1))
    
    return tau

    

def calc_lgN(tex,tau,sigma_v,model): 
    '''
    To calculate lgN, see spectra_rot.md for the detailed algorithm.
    
    Inputs
    ------
    tex: [K], scalar or nD array, excitation temperature
    tau: scalar or nD array with the same shape as tex, mean opacity
    sigma_v: [km/s], scalar or nD array with the same shape as tex, 
        1-sigma velocity dispersion
    model: string, name of the transition, see Splatalog for possible values
    
    Returns
    -------
    lgN: log10(N/cm^-2)
    '''

    tmp = Splatalog[model]

    sw_ls = tmp['sw_ls']
    C_N = tmp['C_N']
    C_T1 = tmp['C_T1']
    C_T2 = tmp['C_T2']
    
    # Q_rot
    Q_rot = calc_Q(model,tex,sw_ls)

    # lgN
    lgN = np.log10(2*(np.log(2)/np.pi)**.5*tau/np.sum(C_N)*sigma_v*Q_rot*
        np.exp(C_T1/tex)/(np.exp(C_T2/tex)-1))

    return lgN

    

def calc_lgN_ch(tex,tau,dv,model): 
    '''
    To calculate lgN in a velocity range (without the assumption of
    Gaussian velocity profile) given tau, tex, dv. Useful for e.g. deriving 
    CO column density given tau (from 12CO/13CO data) and 12CO (2-1) intensity. 
    See spectra_rot.md for the detailed algorithm.
    
    Inputs
    ------
    tex: [K], scalar or nD array, excitation temperature
    tau: scalar or nD array with the same shape as tex, mean opacity
    dv: [km/s], scalar or nD array with the same shape as tex, velocity width
    model: string, name of the transition, see Splatalog for possible values
    
    Returns
    -------
    lgN: log10(N/cm^-2)

    Notes
    -----
    Actually, dv is equivalent to 2(2ln2)^.5*sigma_v, i.e.
    calc_lgN_ch(tex,tau,dv,model) = calc_lgN(tex,tau,dv/2(2ln2)^.5,model)
    '''

    return calc_lgN(tex,tau,dv/(2*(2*np.log(2))**.5),model)



def calc_tau_iso_ratio(r_I,r_N):
    '''
    To calculate the opacity given the line ratio & the column density ratio of 
    two isotopic species, i.e. to solve the equation 
    I1/I2 = (1-e^-tau1)/(1-e^-tau2). The two transitions must be identical so 
    that tau is proportional only to the column density.

    Assumptions made in this calculation:
    1. LTE for the two species along the LoS, i.e. the radiative transfer 
    equation is analytically integrable with a uniform Tex along LoS.
    2. tau is proportional only to the column density. Actually the parameters 
    of the transitions of two isotopic species have some slight differences. 
    
    Algorithm:
    Let r_I=I1/I2, r_N=tau1/tau2, x=e^(-tau1/r_N). Since I1, I2, tau1, tau2 >0, 
    we have r_I, r_N >0, 0<x<1. We want to solve the equation 
    f(x) = x^r_N - r_I*x + r_I - 1 = 0.
    The equation has solutions in these domains only
    if 1<r_I<r_N or r_N<r_I<1. Since x=1 is always a solution to the 
    equation, the initial guess 
    x0 = epsilon*x|[df/dx=0] = epsilon*(r_I/r_N)^(1/(r_N-1))
    to avoid this solution, where 0<epsilon<1.

    Inputs
    ------
    r_I: I1/I2
    r_N: tau1/tau2

    Outputs
    -------
    tau1: the opacity of the 1st species
    '''

    tau1 = np.nan 
    if 1<r_I<r_N or r_N<r_I<1:
        epsilon = .5
        x0 = epsilon*(r_I/r_N)**(1/(r_N-1)) # initial guess
        x = fsolve(lambda x, r_I, r_N: x**r_N-r_I*x+r_I-1,x0,args=(r_I,r_N))
        tau1 = -np.log(x)*r_N

    return tau1



def test():

    ''' calc_tau_iso_ratio with different initial guesses
    # Conclusion: solution is good as long as epsilon < 1.

    # parameters
    r_I = 1.1
    r_N = 8
    X0 = np.linspace(0,1,100) # initial guess

    xp = (r_I/r_N)**(1/(r_N-1))
    Tau1 = np.full_like(X0,np.nan)
    for i in range(len(X0)):
        x = fsolve(lambda x, r_I, r_N: x**r_N-r_I*x+r_I-1,X0[i],args=(r_I,r_N))
        Tau1[i] = -np.log(x)*r_N

    plt.figure()
    plt.scatter(X0,Tau1)
    plt.axvline(xp)
    plt.xlabel('epsilon')
    plt.ylabel('tau1')
    plt.grid()
    plt.show()
    #'''









