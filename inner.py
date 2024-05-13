import time
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.patches import Circle, Ellipse

def circle(x, y, x0, y0, r, boolean='union'):
    '''
    To check whether point (x,y) is in circles.

    Inputs
    ------
    x,y: scalar or np.ndarray of the same shape, coords of points
    x0,y0: scalar or 1darray, coords of centers of circles
    r: scalar or 1darray, radii
    boolean: 'intersection' or 'union', boolean operation on circles
        
    Returns
    -------
    res: bools of the same shape as x, True for in the region
    '''
    # make scalars as arrays
    if np.isscalar(x):
        x, y = np.array([x]), np.array([y])
    if np.isscalar(x0):
        x0, y0, r = np.array([x0]), np.array([y0]), np.array([r])
        
    # main procedure
    if boolean == 'union':
        res = np.zeros_like(x, bool)
        for i in range(len(x0)):
            res = res + (((x-x0[i])**2+(y-y0[i])**2)**.5 < r[i])
    elif boolean=='intersection':
        res = np.ones_like(x, bool)
        for i in range(len(x0)):
            res = res * (((x-x0[i])**2+(y-y0[i])**2)**.5 < r[i])
    
    if np.isscalar(x):
        res = res[0]    
    return res


def ellipse(x, y, x0, y0, a, b, alpha):
    '''
    To check if point (x,y) is in the eclipse.

    Inputs
    ------
    x,y: scalar or np.ndarray of the same shape, coord of the points
    x0,y0: scalar, coord of the center
    a, b: scalar, the major and minor axis of the eclipse
    alpha: [deg], scalar, the angle between the major axis and x-axis
    
    Returns
    -------
    res: bools of the same shape as x, True for in the region
    '''
    alpha *= np.pi/180
    c = np.cos(alpha)
    s = np.sin(alpha)
    x1 =  (x-x0)*c + (y-y0)*s
    y1 = -(x-x0)*s + (y-y0)*c
    return (x1/a)**2 + (y1/b)**2 < 1


def polygon(x, y, coords):
    '''
    To check if the point (x,y) is in the polygon.
    
    Inputs
    ------
    x, y: scalar or np.ndarray of the same shape, coord of the points
    coords: N*2 array, coords of the vertices of the polygon [x,y]

    Returns
    -------
    res: bools of the same shape as x, True for in the region
    '''
    if np.isscalar(x):
        return Path(coords).contains_point([x,y])
    else:
        shape = x.shape 
        x = x.flatten()
        y = y.flatten()
        points = np.transpose([x, y])
        return Path(coords).contains_points(points).reshape(shape)


def sector_ring(x, y, x0, y0, r1, r2, th1, th2):
    '''
    Check if point (x,y) is in the sector ring.

    Inputs
    ------
    x, y: scalar or np.ndarray of the same shape, coord of the points
    x0, y0: scalar, center of the sector ring 
    r1, r2: scalar, minor & major radius of the ring
    th1, th2: [deg], starting & ending angles of the sector, [0, 360], 
              anticlockwise from x-axis
    
    Returns
    -------
    res: bools of the same shape as x, True for in the region
    '''
    th1 *= np.pi/180
    th2 *= np.pi/180
    r = np.hypot(x-x0, y-y0)
    th = np.arctan2(y-y0, x-x0)  # [-pi, pi]
    th[th<0] = th[th<0] + 2*np.pi  # [0, 2*pi]

    m1 = r > r1
    m2 = r < r2
    m3 = th > th1
    m4 = th < th2

    if th1 <= th2:
        res = np.logical_and.reduce((m1,m2,m3,m4))
    else:
        res = np.logical_and.reduce((m1,m2,np.logical_or(m3,m4)))

    return res


def elliptic_sector_ring(x, y, x0, y0, a1, a2, b1, b2, theta, th1, th2):
    '''
    Check if point (x,y) is in the elliptic sector ring.

    Inputs
    ------
    x, y: scalar or np.ndarray of the same shape, coord of the points
    x0, y0: scalar, center of the ellipses
    a1, a2, b1, b2: scalar, major and minor axes of the ellipse
    theta: [deg], scalar, tilt angle, between a- & x-axis
    th1, th2: [deg], position angles of the sector, [0,360]

    Returns
    -------
    res: bools of the same shape as x, True for in the region
    '''
    theta *= np.pi/180
    th1 *= np.pi/180
    th2 *= np.pi/180

    c = np.cos(theta)
    s = np.sin(theta)

    x1 =  (x-x0)*c + (y-y0)*s
    y1 = -(x-x0)*s + (y-y0)*c
    r1 = (x1/a1)**2 + (y1/b1)**2
    r2 = (x1/a2)**2 + (y1/b2)**2
    th = np.arctan2(y-y0, x-x0)
    th[th<0] = th[th<0] + 2*np.pi

    m1 = r1 > 1
    m2 = r2 < 1
    m3 = th > th1
    m4 = th < th2

    if th1 <= th2:
        res = np.logical_and.reduce((m1,m2,m3,m4))
    else:
        res = np.logical_and.reduce((m1,m2,np.logical_or(m3,m4)))

    return res


def test():
    # parameters 
    x = np.linspace(-10, 10, 450)
    y = np.linspace(-10, 10, 500)
    X, Y = np.meshgrid(x, y)

    ''' circle ...............................................................
    x0 = np.random.uniform(-10, 10, 10)
    y0 = np.random.uniform(-10, 10, 10)
    r = np.full_like(x0, 3)

    t = time.time()
    res = circle(X, Y, x0, y0, r, 'union')
    t = time.time() - t
    print(f'Time: {1e3*t:.2f} ms.')

    plt.figure()
    ax = plt.gca()
    plt.imshow(res, origin='lower', extent=[x.min(),x.max(),y.min(),y.max()],
               alpha=.7)
    plt.colorbar()
    for i in range(len(x0)):
      ax.add_patch(Circle((x0[i],y0[i]), r[i], fc='none', ec='k'))
    plt.grid()
    plt.xlabel('x')
    plt.ylabel('y')
    plt.tight_layout()
    plt.savefig('/home/dev/Documents/coding/module/image/inner/circle.pdf')
    plt.close()
    #'''

    ''' ellipse ..............................................................
    x0 = np.random.uniform(-10, 10)
    y0 = np.random.uniform(-10, 10)
    a = np.random.uniform(0, 10)
    b = np.random.uniform(0, 10)
    alpha = np.random.uniform(0, np.pi)

    t = time.time()
    res = ellipse(X, Y, x0, y0, a, b, alpha)
    t = time.time() - t
    print(f'Time: {1e3*t:.2f} ms.')

    plt.figure()
    ax = plt.gca()
    plt.imshow(res, origin='lower', extent=[x.min(),x.max(),y.min(),y.max()],
               alpha=.7)
    plt.colorbar()
    ax.add_patch(Ellipse((x0,y0), 2*a, 2*b, alpha*180/np.pi, 
                         fc='none', ec='k'))
    plt.grid()
    plt.xlabel('x')
    plt.ylabel('y')
    plt.tight_layout()
    plt.savefig('/home/dev/Documents/coding/module/image/inner/ellipse.pdf')
    plt.close()
    #'''

    ''' polygon ..............................................................
    coords = np.random.uniform(-10, 10, (10,2))

    t = time.time()
    res = polygon(X, Y, coords)
    t = time.time() - t
    print(f'Time: {1e3*t:.2f} ms.')

    plt.figure()
    ax = plt.gca()
    plt.imshow(res, origin='lower', extent=[x.min(),x.max(),y.min(),y.max()],
               alpha=.7)
    plt.colorbar()
    plt.plot(coords[:,0], coords[:,1], color='k', lw=1)
    plt.grid()
    plt.xlabel('x')
    plt.ylabel('y')
    plt.tight_layout()
    plt.savefig('/home/dev/Documents/coding/module/image/inner/polygon.pdf')
    plt.close()
    #'''

    ''' sector_ring ..........................................................
    x0 = np.random.uniform(-5, 5)
    y0 = np.random.uniform(-5, 5)
    r1 = np.random.uniform(0, 3)
    r2 = np.random.uniform(3, 5)
    th1 = np.random.uniform(0, 360)
    th2 = np.random.uniform(0, 360)
    
    t = time.time()
    res = sector_ring(X, Y, x0, y0, r1, r2, th1, th2)
    t = time.time() - t
    print(f'Time: {1e3*t:.2f} ms.')

    text = f'Starting angle: {th1:.1f}°\nEnding angle: {th2:.1f}°'

    plt.figure()
    ax = plt.gca()
    plt.imshow(res, origin='lower', extent=[x.min(),x.max(),y.min(),y.max()],
               alpha=.7)
    plt.text(.1, .9, text, fontsize=12, transform=ax.transAxes)
    plt.colorbar()
    plt.grid()
    plt.xlabel('x')
    plt.ylabel('y')
    plt.tight_layout()
    plt.savefig('/home/dev/Documents/coding/module/image/inner/'
                'sector_ring.pdf')
    plt.close()
    #'''
    
    ''' elliptic_sector_ring .................................................
    x0 = np.random.uniform(-5, 5)
    y0 = np.random.uniform(-5, 5)
    a1 = np.random.uniform(0, 3)
    a2 = np.random.uniform(3, 5)
    b1 = np.random.uniform(0, 3)
    b2 = np.random.uniform(3, 5)
    theta = np.random.uniform(0, 180)
    th1 = np.random.uniform(0, 360)
    th2 = np.random.uniform(0, 360)
    
    t = time.time()
    res = elliptic_sector_ring(X, Y, x0, y0, a1, a2, b1, b2, theta, th1, th2)
    t = time.time() - t
    print(f'Time: {1e3*t:.2f} ms.')

    text = f'Starting angle: {th1:.1f}°\nEnding angle: {th2:.1f}°'

    plt.figure()
    ax = plt.gca()
    plt.imshow(res, origin='lower', extent=[x.min(),x.max(),y.min(),y.max()],
               alpha=.7)
    plt.text(.1, .9, text, fontsize=12, transform=ax.transAxes)
    plt.colorbar()
    plt.grid()
    plt.xlabel('x')
    plt.ylabel('y')
    plt.tight_layout()
    plt.savefig('/home/dev/Documents/coding/module/image/inner/'
                'elliptic_sector_ring.pdf')
    plt.close()
    #'''
    