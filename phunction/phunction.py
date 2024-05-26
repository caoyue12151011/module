'''
Contents
--------
Blackbody
    SED_bb: blackbody SED
    SED2color: SED to RGB color
    T2color_bb: blackbody temperature to RGB color

Main-sequence star
    m2L_MS: mass to luminosity 
    L2m_MS: luminosity to mass 
    m2R_MS: mass to radius 
    m2T_MS: mass to surface temperature 
    m2color_MS: mass to RGB color 

Interstellar dust
    SED_dust
    
Water
    T2rho_H2O: temperature to density 
    T2surface_tension_H2O: temperature to surface tension 
    T2vapor_P_H2O: temperature to vapor pressure  
    T2mu_H2O: temperature to dynamic viscosity  
    vapor_P2T_H2O: vapor pressure to temperature  

Earth's atmosphere
    h2property_atm: altitude to properties

Earth
    r2rho_earth: radius to density 
    r2g_earth: radius to gravity acceleration

demo: demonstrate the relations
'''
import os
import colour
import numpy as np 
import astropy.units as u 
import astropy.constants as ct
import matplotlib.pyplot as plt 
from ambiance import Atmosphere as Atm 
from matplotlib.widgets import Slider, Button

# change default matplotlib fonts
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12
plt.rcParams['axes.labelsize'] = 16
plt.rcParams['legend.fontsize'] = 13

# get path of this script
path = os.path.dirname(os.path.abspath(__file__))

# load data 
Data = {}
name = 'earth_density_vs_radius'
r, rho = np.loadtxt(f'{path}/data/{name}.txt').transpose()
Data[name] = {'r': r, 'rho': rho}

name = 'H2O_density_vs_T'
T, rho = np.loadtxt(f'{path}/data/{name}.txt').transpose()
Data[name] = {'T': T, 'rho': rho}

name = 'H2O_surface_tension_vs_T'
T, sigma = np.loadtxt(f'{path}/data/{name}.txt').transpose()
Data[name] = {'T': T, 'sigma': sigma}

name = 'H2O_vapor_P_viscosity_vs_T'
T, P, mu = np.loadtxt(f'{path}/data/{name}.txt').transpose()
Data[name] = {'T': T, 'P': P, 'mu': mu}

# blackbody -------------------------------------------------------------------

# constants
C_T = (ct.h*ct.c/ct.k_B).to_value(u.K*u.um)
C_I = (2*ct.h*ct.c).to_value(u.um**3*u.Jy)

def SED_bb(lamda, T):
    """
    To effeciently calculate the blackbody SED. 

    Inputs
    ------
    lamda: [um], scalar or 1darray, wavelength
    T: [K], temperature

    Returns
    -------
    B_nu: [Jy/sr], intensity
    """
    return C_I / lamda**3 / (np.exp(C_T/T/lamda) - 1)


def SED2color(lamda, I_lamda):
    '''
    Derive the RGB color given the SED.

    Inputs
    ------
    lamda: [nm], 1darray, wavelength. lamda[0]<780, lamda[-1]>360, and the 
        interval must be 1, 5, 10 or 20nm
    I_lamda: [any unit], spectral energy per wavelength, i.e. dI/d(lamda)

    Returns
    -------
    RGB: [0~1], array of length 3. 
    '''
    # check inputs 
    if lamda[0] >= 780 or lamda[-1] <= 360:
        print('Error: wavelengths not in the legit range [360,780] nm.')

    # interpolate SED
    lamda0 = np.linspace(360, 780, 421)
    I_lamda0 = np.interp(lamda0, lamda, I_lamda)

    sp = colour.SpectralDistribution(I_lamda0, lamda0)
    RGB = colour.XYZ_to_sRGB(colour.sd_to_XYZ(sp))
    RGB /= max(RGB)
    RGB[RGB<0.] = 0.

    return RGB


def T2color_bb(T):
    '''
    Derive the RGB color array of a blackbody of temperature T. Nan inputs give
    white color.

    Inputs
    ------
    T: [K], scalar or 1d array

    Returns
    -------
    RGB: [0~1], 3 or (len(T),3) array
    '''
    is_scalar = not hasattr(T, "__len__")
    if is_scalar:
        T = np.array([T])
    
    SED = [colour.sd_blackbody(t) for t in T]
    RGB = np.array([colour.XYZ_to_sRGB(colour.sd_to_XYZ(sed)) for sed in SED])
    _max = np.max(RGB,axis=1)
    RGB /= _max[:,np.newaxis]
    RGB[RGB<0.] = 0.

    if is_scalar:
        RGB = RGB[0]
    return RGB

# main-sequence star ----------------------------------------------------------

def m2L_MS(m):
    '''
    To calculate the luminosity of a main sequence star given its mass. See
    https://en.wikipedia.org/wiki/Mass%E2%80%93luminosity_relation.

    Inputs
    ------
    m: [Msun], scalar or np.array

    Returns
    -------
    L: [Lsun], scalar or np.array
    '''
    return (.23*m**2.3 *(m<.43) + m**4 *(.43<=m)*(m<2) +
            1.4*m**3.5 *(2<=m)*(m<55) + 32000*m *(55<=m) )
    

def L2m_MS(L):
    # same as m2L_MS() but reversed
    return ((L/.23)**(1/2.3) *(L<.033) +
            L**.25 *(.033<=L)*(L<1.189) +
            (L/1.4)**(1/3.5) *(1.189<=L)*(L<1727418) +
            L/32000 *(1727418<=L))
    

def m2R_MS(m):
    '''
    To calculate the radius of a main sequence star given its mass. See 
    http://personal.psu.edu/rbc3/A534/lec18.pdf.

    Inputs
    ------
    m: [Msun], scalar or np.array

    Returns
    -------
    R: [Rsun], scalar or np.array
    '''
    return m**.8 *(m<1) + m**.57 *(1<=m)
    

def m2T_MS(m):
    '''
    To calculate the surface temperature of a main sequence star given its 
    mass.

    Inputs
    ------
    m: [Msun], scalar or np.array

    Returns
    -------
    T: [K], scalar or np.array
    '''
    L = m2L_MS(m)
    R = m2R_MS(m)
    T = ((L*u.L_sun/(4*np.pi*ct.sigma_sb*(R*u.R_sun)**2))**.25).to_value(u.K)

    return T


def m2color_MS(m):
    '''
    Derive the RGB color array of a main-sequence star given its mass, assuming
    blackbody radiation. Nan inputs give white color.

    Inputs
    ------
    m: [M_sun], scalar or 1darray

    Returns
    -------
    RGB: of shape 3 or (len(m),3)
    '''
    return T2color_bb(m2T_MS(m))

# Interstellar dust -----------------------------------------------------------

# constants 
mu_H2 = 2.8  # molecular mass of ISM per hydrogen atom
kappa0 = 10*u.cm**2/u.g
lamda0 = 300  # [um]

# derived constants
C_tau = (mu_H2*ct.m_n*kappa0).to_value(u.cm**2)

def SED_dust(lamda, T, lgN, beta, Gamma, sw_thin=False):
    '''
    To calculate the SED of interstellar dust using the graybody model. 
    This function should be fast since it'll be used for pixel-by-pixel SED 
    fitting. N_H2 and I_nu (see below) are logarithmized since they span 
    orders of magnitude and may cause error in the SED fitting.

    Inputs
    ------
    lamda: [um], scalar or 1darray, wavelength
    T: [K], scalar, dust temperature
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
    B_nu = SED_bb(lamda, T)
    if sw_thin:
        lgI = np.log10(B_nu * tau)
    else:   
        lgI = np.log10(B_nu * (1-np.exp(-tau)))

    return lgI

# Water -----------------------------------------------------------------------

def T2rho_H2O(T): 
    ''' Temperature to density at 1 atm.  

    Inputs
    ------
    T: [°C], scalar or array. [-100,100]°C, otherwise returns nan

    Returns
    -------
    rho: [kg/m^3], scalar or array 
    '''
    data = Data['H2O_density_vs_T']
    return np.interp(T, data['T'], data['rho'], np.nan, np.nan)


def T2surface_tension_H2O(T): 
    ''' Temperature to surface tension.  

    Inputs
    ------
    T: [°C], scalar or array. [0,374]°C, otherwise returns nan 

    Returns
    -------
    S: [mN/m^2], scalar or array 
    '''
    data = Data['H2O_surface_tension_vs_T']
    return np.interp(T, data['T'], data['sigma'], np.nan, np.nan)


def T2vapor_P_H2O(T): 
    ''' Temperature to vapor pressure.  

    Inputs
    ------
    T: [°C], scalar or array. [0.01,360]°C, otherwise returns nan 

    Returns
    -------
    P: [bar], scalar or array 
    '''
    data = Data['H2O_vapor_P_viscosity_vs_T']
    return np.interp(T, data['T'], data['P'], np.nan, np.nan)


def vapor_P2T_H2O(P):
    data = Data['H2O_vapor_P_viscosity_vs_T']
    return np.interp(P, data['P'], data['T'], np.nan, np.nan)


def T2mu_H2O(T): 
    ''' Temperature to dynamic viscosity.  

    Inputs
    ------
    T: [°C], scalar or array. [0.01,360]°C, otherwise returns nan 

    Returns
    -------
    mu: [mPa s], scalar or array 
    '''
    data = Data['H2O_vapor_P_viscosity_vs_T']
    return np.interp(T, data['T'], data['mu'], np.nan, np.nan)

# Earth's atmosphere ----------------------------------------------------------

def h2property_atm(h, prop):
    ''' Altitude to properties.

    h: [m], scalar or array 
    prop: string, returns
        'p': [bar], pressure
        'rho': [km/m^3]
        'T': [°C]
        'mu': [kg/m/s], dynamic viscosity
        'cs': [m/s], sound speed
    '''
    model = Atm(h, check_bounds=False)
    res = None 
    if prop == 'p':
        res = model.pressure/1e5
    elif prop == 'rho':
        res = model.density
    elif prop == 'T':
        res = model.temperature_in_celsius
    elif prop == 'mu':
        res = model.dynamic_viscosity
    elif prop == 'cs':
        res = model.speed_of_sound

    if not hasattr(h, '__len__'):
        res = res[0]
    return res

# Earth -----------------------------------------------------------------------

def r2rho_earth(r):
    ''' Density given radius from the Earth center.  

    Inputs
    ------
    r: [km], scalar or array. [0, 6371] km, otherwise returns nan 

    Returns
    -------
    rho: [g/cm^3], scalar or array 
    '''
    data = Data['earth_density_vs_radius']
    return np.interp(r, data['r'], data['rho'], np.nan, np.nan)


def r2g_earth(r):
    ''' Gravity acceletation given radius from the Earth center.  

    Inputs
    ------
    r: [km], scalar or array. [0, 6371] km, otherwise returns nan 

    Returns
    -------
    g: [m/s^2], scalar or array 
    '''
    is_scalar = not hasattr(r, '__len__')
    if is_scalar:
        r = np.array([r])

    R = np.linspace(0, r, 10000)  # (R, r)
    Rho = r2rho_earth(R)
    g = 4*np.pi*ct.G / r**2 * np.trapz(R**2*Rho, R, axis=0) * u.g/u.cm**3*u.km
    g = g.to_value(u.m/u.s**2)
    g[r==0] = 0  # fix the r=0 part 

    if is_scalar:
        g = g[0]
    return g


def demo():
    # Interstellar dust .......................................................

    #''' SED_dust
    # parameters 
    lamda = np.linspace(.36, .78, 300)  # [um]
    T = 6000  # [K]
    lgN = 15.9
    beta = 3
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
    I = 10**SED_dust(lamda, T, lgN, beta, Gamma)
    I_thin = 10**SED_dust(lamda, T, lgN, beta, Gamma, sw_thin=True)

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
                  valmin=0, valmax=60000, valinit=T, orientation='vertical')
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
        I = 10**SED_dust(lamda, T, lgN, beta, Gamma)
        I_thin = 10**SED_dust(lamda, T, lgN, beta, Gamma, sw_thin=True)

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
    plt.savefig(f'{path}/image/dust_SED.pdf')
    plt.show()
    #'''

    # water ...................................................................

    ''' T2rho_H2O
    T = np.linspace(-100,100,1000)
    rho = T2rho_H2O(T)

    plt.figure(figsize=(7,4.5))
    plt.plot(T,rho,color='k')
    plt.xlabel('Temperature (°C)')
    plt.ylabel('Density '+r'$\rm(g\ cm^{-3})$')
    plt.grid()
    plt.tight_layout()
    plt.savefig('%s/image/phunction/T2rho_H2O.pdf'%dire)
    plt.close()
    #'''

    ''' T2surface_tension_H2O
    T = np.linspace(0,400,1000)
    S = T2surface_tension_H2O(T)

    plt.figure(figsize=(7,4.5))
    plt.plot(T,S,color='k')
    plt.xlabel('Temperature (°C)')
    plt.ylabel('Surface tension '+r'$\rm(mN\ m^{-2})$')
    plt.grid()
    plt.tight_layout()
    plt.savefig('%s/image/phunction/T2surface_tension_H2O.pdf'%dire)
    plt.close()
    #'''

    ''' T2vapor_P_H2O
    T = np.linspace(0,400,1000)
    P = T2vapor_P_H2O(T)

    plt.figure(figsize=(7,4.5))
    plt.yscale('log')
    plt.plot(T,P,color='k')
    plt.xlabel('Temperature (°C)')
    plt.ylabel('Vapor pressure (bar)')
    plt.grid()
    plt.tight_layout()
    plt.savefig('%s/image/phunction/T2vapor_P_H2O.pdf'%dire)
    plt.close()
    #'''

    ''' T2mu_H2O
    T = np.linspace(0,400,1000)
    mu = T2mu_H2O(T)

    plt.figure(figsize=(7,4.5))
    plt.plot(T,mu,color='k')
    plt.xlabel('Temperature (°C)')
    plt.ylabel('Dynamic viscosity '+r'$\rm(mPa\ s)$')
    plt.grid()
    plt.tight_layout()
    plt.savefig('%s/image/phunction/T2mu_H2O.pdf'%dire)
    plt.close()
    #'''

    # Earth's atmosphere ......................................................

    '''
    # parameters 
    h = np.linspace(-6e3,85e3,1000) # [m]
    h_tp = 11019 # [m], tropopause
    h_sp = 47350 # stratopause
    h_mp = 86000 # mesopause
    xlabel = ['Pressure (bar)',r'$\rm Density\ (kg\ m^{-3})$',
        r'$\rm Temperature\ (^\circ C)$',
        r'$\rm Dynamic\ viscosity\ (kg\ m^{-1}s^{-1})$',
        r'$\rm Sound\ speed\ (m\ s^{-1})$']
    name = ['p','rho','T','mu','cs']
    xticks = np.arange(-5,83,5) # [km]
    funs = [h2P_atm,h2rho_atm,h2T_atm,
        h2mu_atm,h2cs_atm]

    xlabels = [str(i) for i in xticks]

    for i in range(len(name)):
        res = funs[i](h)

        plt.figure(figsize=(8,4))
        if name[i] in ['p','rho']:
            plt.yscale('log')

        plt.plot(h/1e3,res,color='k')
        ax = plt.gca()
        plt.text(.6,.9,'ICAO standard atmosphere (1993)',fontsize=12,
            transform=ax.transAxes)

        plt.axvline(0,color='k',ls='--',lw=1,alpha=.7)
        plt.axvline(h_tp/1e3,color='r',ls='--',lw=1,alpha=.7)
        plt.axvline(h_sp/1e3,color='b',ls='--',lw=1,alpha=.7)
        plt.axvline(h_mp/1e3,color='g',ls='--',lw=1,alpha=.7)

        _min = np.nanmin(res)
        plt.text(-3,_min,'See level',color='k',fontsize=12,rotation=90)
        plt.text(h_tp/1e3-3,_min,'Tropopause',color='r',fontsize=12,
        rotation=90)
        plt.text(h_sp/1e3-3,_min,'Stratopause',color='b',fontsize=12,
            rotation=90)
        plt.text(h_mp/1e3-3,_min,'Mesopause',color='g',fontsize=12,rotation=90)

        ax.set_xticks(xticks)
        ax.set_xticklabels(xlabels)
        plt.xlabel('Altitude (km)')
        plt.ylabel(xlabel[i])
        plt.grid()
        plt.tight_layout()
        plt.savefig('%s/image/phunction/%s_atm.pdf'%(dire,name[i]))
        plt.close()
    #'''

    # Earth ...................................................................

    ''' rho vs r
    r = np.linspace(0, 6371, 1000)
    rho = r2rho_earth(r)

    plt.figure(figsize=(7,4.5))
    plt.plot(r, rho, color='k')
    plt.axvline(1216, color='k', ls='--', lw=1)  # out bound of inner core
    plt.text(1216, .5, 'Inner core', fontsize=12, rotation='vertical')
    plt.axvline(3486, color='k', ls='--', lw=1)  # out bound of outer core
    plt.text(3486, .5, 'Outer core', fontsize=12, rotation='vertical')
    plt.axvline(5671, color='k', ls='--', lw=1)  # out bound of rigid mantle
    plt.text(5671, .5, 'Rigid mantle', fontsize=12, rotation='vertical')
    plt.axvline(6271, color='k', ls='--', lw=1)  # out bound of asthenosphere
    plt.text(6271, .5, 'Asthenosphere', fontsize=12, rotation='vertical')
    plt.xlim(xmin=0)
    plt.ylim(ymin=0)
    plt.xlabel(r'$r\rm\ (km)$')
    plt.ylabel(r'$\rho\rm\ (g\ cm^{-3})$')
    plt.grid()
    plt.tight_layout()
    plt.savefig(f'{path}/image/r2rho_earth.pdf')
    plt.close()
    #'''

    ''' g vs r
    r = np.linspace(0, 6371, 1000)
    g = r2g_earth(r)

    plt.figure(figsize=(7,4.5))
    plt.plot(r, g, color='k')
    plt.axvline(1216, color='k', ls='--', lw=1)  # out bound of inner core
    plt.text(1216, .5, 'Inner core', fontsize=12, rotation='vertical')
    plt.axvline(3486, color='k', ls='--', lw=1)  # out bound of outer core
    plt.text(3486, .5, 'Outer core', fontsize=12, rotation='vertical')
    plt.axvline(5671, color='k', ls='--', lw=1)  # out bound of rigid mantle
    plt.text(5671, .5, 'Rigid mantle', fontsize=12, rotation='vertical')
    plt.axvline(6271, color='k', ls='--', lw=1)  # out bound of asthenosphere
    plt.text(6271, .5, 'Asthenosphere', fontsize=12, rotation='vertical')
    plt.xlim(xmin=0)
    plt.ylim(ymin=0)
    plt.xlabel(r'$r\rm\ (km)$')
    plt.ylabel(r'$g\rm\ (m\ s^{-2})$')
    plt.grid()
    plt.tight_layout()
    plt.savefig(f'{path}/image/r2g_earth.pdf')
    plt.close()
    #'''

    # others ..................................................................

    ''' boiling point vs altitude

    h = np.linspace(-5,100,1000) # [km]
    c = 'floralwhite'

    P = h2P_atm(h*1e3) 
    T = vapor_P2T_H2O(P) # [°C]

    fig = plt.figure()
    ax = plt.gca()
    fig.patch.set_facecolor(c)
    ax.set_facecolor(c) 

    plt.plot(h,T,color='gray',lw=3)
    plt.grid()
    plt.xlabel('Altitude (km)')
    plt.ylabel('Boiling point (°C)')
    plt.tight_layout()
    plt.savefig('%s/image/phunction/h_vs_boiling_point.pdf'%dire)
    plt.close()
    #'''
