''' 
Purpose
------- 
To analyze the spectral data with the rotational transition models. 

 
Functions
---------
gaussian: calculate Gaussian function
find_peaks: to find +/- peaks in spectra and estimate the parameters.
calc_Q: calculate the rotational partition function
func: model spectrum with Tex, N, v, sigma_v as parameters
func_amp: Same as func but with amp, Tex, N, v, sigma_v as parameters
func_tau: Same as func but with Tex, Tau, v, sigma_v as parameters
fitter: spectral line fitter
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
import pickle
import numpy as np
import matplotlib.pyplot as plt
import astropy.constants as ct 
import astropy.units as u 
from scipy.optimize import curve_fit
from scipy.ndimage import gaussian_filter1d
from scipy.optimize import fsolve
from matplotlib.widgets import Slider


# load data (change the directory to your absolute directory of the file Splatalog.p)
Splatalog = pickle.load(open('/Users/caoyue/Documents/astro/modules/spectra/'
    'Splatalog.p','rb'))


def gaussian(x,*pa): 
    '''
    To calculate the 1D multi-Gaussian function.
    
    Inputs
    ------
    x: 1D array, independent variable
    pa: 1D array, parameters [a1,b1,c1,a2,...]
    
    Returns
    -------
    y: 1D array
    '''

    a, b, c = np.reshape(pa,(-1,3)).transpose()[...,np.newaxis]
    y = np.sum(a*np.exp(-(x-b)**2/(2*c**2)),axis=0)
    
    return y



def find_peaks(x,y,sign,lim): # ------------------------------------------------
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



def calc_Q(model,tex,sw_ls):
    '''
    To calculate the rotational partition function Q_rot(Tex).

    Inputs
    ------
    model: string, name of the transition
    tex: [K], scalar or ndarray, excitation temperature
    sw_ls: bool, whether the molecule is linear (T) or symmetric top (F)

    Outputs
    -------
    Q_rot: the rotational partition function
    '''

    tmp = Splatalog[model]

    Q_rot = {}
    if sw_ls:
        C_Q1 = tmp['C_Q1']
        Q_rot = 1/C_Q1*tex*np.exp(C_Q1/3/tex)

    else:
        C_Q2 = tmp['C_Q2']
        C_Q3 = tmp['C_Q3']
        C_Q4 = tmp['C_Q4']
        Q_rot = C_Q2*tex**1.5*np.exp(C_Q3/tex)*(1+C_Q4/tex**2)

    return Q_rot



def func(x,model,T_bg,Switch_dict,*para): # --------------------------
    '''
    To generate the spectrum given parameters Tex, N, v, sigma_v.


    Parameters
    ----------
    x: [km/s], 1D array, line-of-sight velocity in the LSR frame. The 
        independent variable for calculating the spectrum

    T_bg: [K], scalar, background brightness temperature

    model: string, name of the transition, see variable Splatalog for available
        models

    Switch_dict: dict of switches for controlling the output. Values are bools.

        'thin': whether to use the optically thin limit
        'intensity': whether to calculate intensity (T) or brightness T (F)
        'hyperfine': whether to output the spectra of hyperfine lines, only 
            valid when thin=True
        'collapse_hyperfine': whether to set the hyperfine velocity offsets to 
            0, used for demonstrating only one line, shouldn't use for 
            quantitative analyses

    *para: 1D array, parameters of the model, extend para if you have 
        multiple velocity components. Format: [tex, lgN1, v0_1, sigma_v1, 
        lgN2, v0_2, sigma_v2, ...], where the universal tex is in K, 
        lgNi, v0_i, sigma_vi are log10(column density) [cm^-2] of the TRACER, 
        bulk velocity [km/s] and 1-sigma velocity dispersion [km/s] of the ith 
        velocity component.


    Returns
    -------
    y: the modeled spectrum.
        If sw_thin & sw_hyper, y is 2D array with shape (n_hyper+1, len(x)),
        where y[0] is the sum of all the hyperfine spectra. Otherwise, 
        y is 1D array with the same length as x's.
        If sw_I=True y is the intensity in Jy/sr, else y is the brightness 
        temperature in K.


    Notices
    -------
    See documentation spectra_rot.md for details of the models.
    '''


    para = np.array(para)
    tmp = Splatalog[model]

    # switches
    sw_thin = Switch_dict['thin']
    sw_I = Switch_dict['intensity']
    sw_hyper = Switch_dict['hyperfine']
    sw_collapse_hyper = Switch_dict['collapse_hyperfine']


    # component parameters, shape = (n_component,1)
    tex = para[0] # [K]
    lgN, v0, sigma_v = np.reshape(para[1:],(-1,3)).transpose()[...,np.newaxis]

    # structure parameters
    sw_ls = tmp['sw_ls']
    C_N = tmp['C_N']
    C_T1 = tmp['C_T1']
    C_T2 = tmp['C_T2']
    C_I = tmp['C_I']
    v = tmp['v']

    # whether collapse hyperfine lines
    if sw_collapse_hyper:
        v = np.zeros_like(v)

    # reshape
    C_N = C_N[...,np.newaxis,np.newaxis] # shape = (n_hyper, 1, 1)
    v = v[...,np.newaxis,np.newaxis] # shape = (n_hyper, 1, 1)


    # Q_rot
    Q_rot = calc_Q(model,tex,sw_ls)

    # opacity
    TAU = (C_N*10**lgN/sigma_v/Q_rot*
        np.exp(-(x-v0-v)**2/(2*sigma_v**2)-C_T1/tex)*(np.exp(C_T2/tex)-1)) 
        # shape = (n_hyper, n_component, len(x))

    if sw_thin and sw_hyper:
        tau = np.sum(TAU,axis=1) # shape = (n_hyper, len(x))
        tau = np.concatenate((np.sum(tau,axis=0)[np.newaxis],tau))
             # shape = (n_hyper+1, len(x))
    else:
        tau = np.sum(TAU,axis=(0,1)) # shape = x.shape

    # terms
    term_tau = {}
    if sw_thin:
        term_tau = tau 
    else:
        term_tau = (1-np.exp(-tau))

    term_T = 1/(np.exp(C_T2/tex)-1) - 1/(np.exp(C_T2/T_bg)-1)

    # line intensity
    y = {}
    if sw_I:
        y = C_I*term_T*term_tau
            # intensity in Jy/sr
    else:
        y = C_T2/np.log(1+1/term_tau/term_T)
            # Brightness temperature in K

    return y



def func_amp(x,model,T_bg,sw_thin,sw_I,sw_hyper,*para): # ----------------------
    '''
    To generate the spectrum given parameters amp, Tex, N, v, sigma_v. 


    Parameters
    ----------
    x: [km/s], 1D array, line-of-sight velocity in the LSR frame. The 
        independent variable for calculating the spectrum
    T_bg: [K], scalar, background brightness temperature
    model: string, name of the transition, see variable Splatalog for available
        models
    sw_thin: bool, whether to use the optically thin limit
    sw_I: bool, whether to calculate intensity (T) or brightness T (F)
    sw_hyper: bool, whether to output the spectra of hyperfine lines. Only 
        valid when sw_thin=True
    *para: 1D array, parameters of the model, extend para if you have 
        multiple velocity components. Format: [amp1, tex, lgN1, v0_1, sigma_v1, 
        lgN2, v0_2, sigma_v2, ...], where the universal tex is in K, 
        lgNi, v0_i, sigma_vi are log10(column density) [cm^-2] of the TRACER, 
        bulk velocity [km/s] and 1-sigma velocity dispersion [km/s] of the ith 
        velocity component, the universal amp=B_nu(tex') is in Jy/sr.


    Returns
    -------
    y: the modeled spectrum.
        If sw_thin & sw_hyper, y is 2D array with shape (n_hyper+1, len(x)),
        where y[0] is the sum of all the hyperfine spectra. Otherwise, 
        y is 1D array with the same length as x's.
        If sw_I=True y is the intensity in Jy/sr, else y is the brightness 
        temperature in K.


    Notices
    -------
    See documentation spectra_rot.md for details of the models.
    '''


    para = np.array(para)
    tmp = Splatalog[model]


    amp, tex = para[:2] 

    # component parameters, shape = (n_component,1)
    lgN, v0, sigma_v = np.reshape(para[2:],(-1,3)).transpose()[...,np.newaxis]

    # structure parameters
    sw_ls = tmp['sw_ls']
    C_N = tmp['C_N'][...,np.newaxis,np.newaxis] # shape = (n_hyper, 1, 1)
    C_T1 = tmp['C_T1']
    C_T2 = tmp['C_T2']
    C_I = tmp['C_I']
    v = tmp['v'][...,np.newaxis,np.newaxis] # shape = (n_hyper, 1, 1)

    # Q_rot
    Q_rot = calc_Q(model,tex,sw_ls)

    # opacity
    TAU = (C_N*10**lgN/sigma_v/Q_rot*
        np.exp(-(x-v0-v)**2/(2*sigma_v**2)-C_T1/tex)*(np.exp(C_T2/tex)-1)) 
        # shape = (n_hyper, n_component, len(x))

    if sw_thin and sw_hyper:
        tau = np.sum(TAU,axis=1) # shape = (n_hyper, len(x))
        tau = np.concatenate((np.sum(tau,axis=0)[np.newaxis],tau))
             # shape = (n_hyper+1, len(x))
    else:
        tau = np.sum(TAU,axis=(0,1)) # shape = x.shape

        
    # terms
    term_tau = {}
    if sw_thin:
        term_tau = tau 
    else:
        term_tau = (1-np.exp(-tau))

    term_T = 1/(np.exp(C_T2/T_bg)-1)

    # line intensity
    y = {}
    if sw_I:
        y = (amp - C_I*term_T)*term_tau # intensity in Jy/sr
    else:
        y = C_T2/np.log(1+1/term_tau/(amp/C_I-term_T)) 
            # Brightness temperature in K

    return y



def func_tau(x,model,T_bg,Switch_dict,*para): # --------------------------------
    '''
    To generate the spectrum given parameters Tex, Tau, v, sigma_v.


    Parameters
    ----------
    x: [km/s], 1D array, line-of-sight velocity in the LSR frame. The 
        independent variable for calculating the spectrum

    model: string, name of the transition, see variable Splatalog for available
        models

    T_bg: [K], scalar, background brightness temperature

    Switch_dict: dict of switches for controlling the output. Values are bools.

        'thin': whether to use the optically thin limit
        'intensity': whether to calculate intensity (T) or brightness T (F)
        'hyperfine': whether to output the spectra of hyperfine lines, only 
            valid when thin=True
        'collapse_hyperfine': whether to set the hyperfine velocity offsets to 
            0, used for demonstrating only one line, shouldn't use for 
            quantitative analyses

    *para: 1D array, parameters of the model, extend para if you have 
        multiple velocity components. Format: [tex, Tau1, v0_1, sigma_v1, 
        Tau2, v0_2, sigma_v2, ...], where the universal tex is in K, 
        Taui, v0_i, sigma_vi are mean opacity, 
        bulk velocity [km/s] and 1-sigma velocity dispersion [km/s] of the ith 
        velocity component.


    Returns
    -------
    y: the modeled spectrum.
        If sw_thin & sw_hyper, y is 2D array with shape (n_hyper+1, len(x)),
        where y[0] is the sum of all the hyperfine spectra. Otherwise, 
        y is 1D array with the same length as x's.
        If sw_I=True y is the intensity in Jy/sr, else y is the brightness 
        temperature in K.


    Notices
    -------
    See documentation spectra_rot.md for details of the models.
    '''


    para = np.array(para)
    tmp = Splatalog[model]

    # switches
    sw_thin = Switch_dict['thin']
    sw_I = Switch_dict['intensity']
    sw_hyper = Switch_dict['hyperfine']
    sw_collapse_hyper = Switch_dict['collapse_hyperfine']


    # component parameters, shape = (n_component,1)
    tex = para[0] # [K]
    tau0, v0, sigma_v = np.reshape(para[1:],(-1,3)).transpose()[...,np.newaxis]

    # structure parameters
    sw_ls = tmp['sw_ls']
    C_N = tmp['C_N']
    C_T1 = tmp['C_T1']
    C_T2 = tmp['C_T2']
    C_I = tmp['C_I']
    v = tmp['v']
    R = tmp['R']

    # whether collapse hyperfine lines
    if sw_collapse_hyper:
        v = np.zeros_like(v)

    # reshape
    C_N = C_N[...,np.newaxis,np.newaxis] # shape = (n_hyper, 1, 1)
    v = v[...,np.newaxis,np.newaxis] # shape = (n_hyper, 1, 1)
    R = R[...,np.newaxis,np.newaxis] # shape = (n_hyper, 1, 1)


    # Q_rot
    Q_rot = {}
    if sw_ls:
        C_Q1 = tmp['C_Q1'] 
        Q_rot = 1/C_Q1*tex*np.exp(C_Q1/3/tex)

    else:
        C_Q2 = tmp['C_Q2']
        C_Q3 = tmp['C_Q3']
        C_Q4 = tmp['C_Q4']
        Q_rot = C_Q2*tex**1.5*np.exp(C_Q3/tex)*(1+C_Q4/tex**2)

    # opacity
    TAU = 2*(np.log(2)/np.pi)**.5*R*tau0*np.exp(-(x-v0-v)**2/(2*sigma_v**2))
        # shape = (n_hyper, n_component, len(x))

    if sw_thin and sw_hyper:
        tau = np.sum(TAU,axis=1) # shape = (n_hyper, len(x))
        tau = np.concatenate((np.sum(tau,axis=0)[np.newaxis],tau))
             # shape = (n_hyper+1, len(x))
    else:
        tau = np.sum(TAU,axis=(0,1)) # shape = x.shape

    # terms
    term_tau = {}
    if sw_thin:
        term_tau = tau 
    else:
        term_tau = (1-np.exp(-tau))

    term_T = 1/(np.exp(C_T2/tex)-1) - 1/(np.exp(C_T2/T_bg)-1)

    # line intensity
    y = {}
    if sw_I:
        y = C_I*term_T*term_tau
            # intensity in Jy/sr
    else:
        y = C_T2/np.log(1+1/term_tau/term_T)
            # Brightness temperature in K

    return y




def fitter(f,x,y,model,para0,bounds,T_bg,Switch_dict): # ---------------
    '''
    To fit the spectrum with the model in func.

    Parameters
    ----------
    f: the model function
    x: [km/s], 1D array, los velocity of the spectrum
    y: 1D array, the spectrum
    model: string, name of the model
    para0: 1D array, initial guesses
    bounds: ([p1_min, p2_min, ...], [p1_max, p2_max, ...])
    T_bg, Switch_dict: cf. spectra_rot.func

    Returns
    --------
    success: {True|False}, whether fitting is successful
    para: 1D array, optimized parameters
    perr: 1D array, uncertainties of parameters
    '''

    try:
        para, pcov = curve_fit(lambda x,*para: f(x,model,T_bg,Switch_dict,
            *para),x,y,para0,bounds=bounds) 
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
    2. Tau is proportional only to the column density. Actually the parameters 
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








