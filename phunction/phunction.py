'''
Contents
--------
Main-sequence star
    m2L_MS: mass to luminosity 
    L2m_MS: luminosity to mass 
    m2R_MS: mass to radius 
    m2T_MS: mass to surface temperature 
    SED2color: SED to RGB color
    T2color_bb: blackbody temperature to RGB color
    m2color_MS: mass to RGB color 

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
import matplotlib.pyplot as plt 
import astropy.units as u 
import astropy.constants as ct
from ambiance import Atmosphere as Atm 

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


def SED2color(SED):
    '''
    Derive the RGB color given the SED.

    Inputs
    ------
    SED: SED object or array of SED object. See the doc of colour.sd_to_XYZ. 
        For example, colour.sd_blackbody(T) generates an SED object.

    Returns
    -------
    RGB: array of shape 3 or (SED_#,3). 
    '''
    RGB = None
    if hasattr(SED,"__len__"):
        RGB = colour.XYZ_to_sRGB(colour.sd_to_XYZ(SED))
        RGB /= max(RGB)
    else:
        RGB = np.array([colour.XYZ_to_sRGB(colour.sd_to_XYZ(sed)) 
            for sed in SED])
        _max = np.max(RGB,axis=1)
        RGB /= _max[:,np.newaxis]

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
    RGB: 3 or (len(T),3) array
    '''
    if hasattr(T, "__len__"):
        return SED2color(colour.sd_blackbody(T))
    else:
        return SED2color([colour.sd_blackbody(t) for t in T])


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
