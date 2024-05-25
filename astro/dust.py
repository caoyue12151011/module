"""
To analyze the properties of insterstellar dust.
"""
import os
import numpy as np
import matplotlib.pyplot as plt 
import astropy.units as u
import astropy.constants as ct
from matplotlib.widgets import Slider, Button

# change default matplotlib fonts (window layouts)
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12
plt.rcParams['axes.labelsize'] = 16
plt.rcParams['legend.fontsize'] = 13

# get path of this script
path = os.path.dirname(os.path.abspath(__file__))

# constants 
mu_H2 = 2.8  # molecular mass of ISM per hydrogen atom
kappa0 = 10*u.cm**2/u.g
lamda0 = 300  # [um]

# constants used in dust_SED()
C_tau = (mu_H2*ct.m_n*kappa0).to_value(u.cm**2)
C_T = (ct.h*ct.c/ct.k_B).to_value(u.K*u.um)
C_I = (2*ct.h*ct.c).to_value(u.um**3*u.Jy)


def dust_SED(lamda, T_dust, lgN, beta, Gamma, sw_thin=False):
    '''
    To calculate the SED of interstellar dust using the graybody model. 
    This function should be fast since it'll be used by pixel-by-pixel SED 
    fitting. N_H2 and I_nu (see below) are logarithmized since they span 
    orders of magnitude and may cause error in SED fitting.

    Inputs
    ------
    lamda: [um], scalar or 1darray, wavelength
    T_dust: [K], scalar, dust temperature
    lgN: = lg(N_H2/cm^-2), scalar, log H2 column density
    beta: scalar, opacity power index
    Gamma: scalar, gas-to-dust mass ratio
    sw_thin: bool, whether using the optically thin limit

    Returns
    -------
    lgI: = lg(I/(Jy/sr)), log intensity
    '''
    tau = C_tau * (lamda/lamda0)**-beta * 10**lgN / Gamma
    tau = tau.astype(float)
    B_nu = C_I / lamda**3 / (np.exp(C_T/T_dust/lamda) - 1)
    if sw_thin:
        lgI = np.log10(B_nu * tau)
    else:   
        lgI = np.log10(B_nu * (1-np.exp(-tau)))

    return lgI


def test():
    ''' dust_SED()
    # parameters 
    lamda = np.linspace(50, 1000, 300)  # [um]
    T = 20  # [K]
    lgN = 20
    beta = 2
    Gamma = 100

    # figure parameters 
    figsize = (8, 5)
    fx_p = .4  # fraction of the plotting panel
    fy_p = .15
    fw_p = .55
    fh_p = .8
    fx_s1 = .05
    fx_s2 = .15
    fx_s3 = .25
    fy_s = .1
    fw_s = .05
    fh_s = .8

    # spectrum
    I = 10**dust_SED(lamda, T, lgN, beta, Gamma)
    I_thin = 10**dust_SED(lamda, T, lgN, beta, Gamma, sw_thin=True)

    # plotting panel
    fig = plt.figure(figsize=figsize)
    ax = fig.add_axes([fx_p, fy_p, fw_p, fh_p])
    plt.yscale('log')
    l, = plt.plot(lamda, I, color='k', label='Real')
    l_thin, = plt.plot(lamda, I_thin, ls='--', color='k', 
                       label='Optically thin')
    plt.legend()
    plt.grid()
    plt.xlabel(r'$\lambda\rm\ (\mu m)$')
    plt.ylabel(r'$\rm Intensity\ (Jy\ sr^{-1})$')
    
    # sliders 
    sld1 = Slider(ax=fig.add_axes([fx_s1,fy_s,fw_s,fh_s]),
                  label=r'$T_{\rm dust}\rm\ (K)$',
                  valmin=0, valmax=200, valinit=T, orientation='vertical')
    sld2 = Slider(ax=fig.add_axes([fx_s2,fy_s,fw_s,fh_s]),
                  label=r'$\lg(N_{\rm H_2}/\rm cm^{-2})$',
                  valmin=5, valmax=30, valinit=lgN, orientation='vertical')
    sld3 = Slider(ax=fig.add_axes([fx_s3,fy_s,fw_s,fh_s]),
                  label=r'$\beta$',
                  valmin=-5, valmax=8, valinit=beta, orientation='vertical')

    # The function to be called anytime a slider's value changes
    def update(val):
        # new slider values
        T = sld1.val
        lgN = sld2.val
        beta = sld3.val 

        # new data
        I = 10**dust_SED(lamda, T, lgN, beta, Gamma)
        I_thin = 10**dust_SED(lamda, T, lgN, beta, Gamma, sw_thin=True)

        ymin = min(np.nanmin(I), np.nanmin(I_thin))
        ymax = max(np.nanmax(I), np.nanmax(I_thin))

        # upload plots
        l.set_ydata(I)
        l_thin.set_ydata(I_thin)
        ax.set_ylim(ymin, ymax)
        fig.canvas.draw_idle()

    # register the update function with each slider
    sld1.on_changed(update)
    sld2.on_changed(update)
    sld3.on_changed(update)
    plt.show()
    #'''




















