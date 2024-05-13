import time
import numpy as np 
import matplotlib.pyplot as plt 
from numba import jit 

def cubic_solver(a, b, c, d):
    """
    Solve cubic equation ax^3 + bx^2 + cx + d = 0.
    a, b, c, d: scalar or np.ndarray of same shape, real numbers 
        (complex numbers not tested).
    Returns 3 roots (can be equal).
    Solving time ~2 µs if scalar. Sympy solving time ~40 ms.
    """
    is_scalar = not hasattr(a, '__len__')

    # np.float to float, otherwise returns nan in fractional power of 
    # negative numbers
    if is_scalar:
        a, b, c, d = float(a), float(b), float(c), float(d)

    w = complex(-.5, 3**.5/2)
    shift = b/(3*a)
    p = (3*a*c - b**2) / (3*a**2)
    q = (2*b**3 - 9*a*b*c + 27*a**2*d) / (27*a**3)
    delta = q**2/4 + p**3/27

    # to avoid loss of significance when p close to 0
    if is_scalar:
        if q < 0:  
            u = (-.5*q + delta**.5)**(1/3)
        else:
            u = (-.5*q - delta**.5)**(1/3)
    else:
        u = np.full(q.shape, np.nan, np.complex_)
        ind = q < 0
        u[ind] = (-.5*q[ind] + delta[ind].astype(np.complex_)**.5)**(1/3)
        ind = q >= 0
        u[ind] = (-.5*q[ind] - delta[ind].astype(np.complex_)**.5)**(1/3)

    v = -p/(3*u)
    x1 = u + v - shift
    x2 = w*(u + w*v) - shift
    x3 = w*(w*u + v) - shift

    return x1, x2, x3


def dichotomy(f,low,high,tol=1e-8,rel_tol=True,*args):
    '''
    To find the root of f(x)=0 using the dichotomy search. 

    Inputs
    ------
    f: the single-variable function 
    low,high: the initial interval for the search, low <= high and 
            f(low)*f(high) <=0.
    tol: tolerance 
    rel_tol: whether the tolerance is relative 
    *args: other arguments in f

    Returns
    -------
    x: the root
    '''
    # initial check 
    if f(low,*args)*f(high,*args)>0:
        raise ValueError('f(low)*f(high)>0 in the dichotomy search.')

    mid = (low+high)/2 
    if rel_tol:
        while abs((high-low)/mid)>tol:
            if f(low,*args)*f(mid,*args)<=0:
                high = mid 
            else:
                low = mid 
            mid = (low+high)/2 
    else:
        while abs(high-low)>tol:
            if f(low,*args)*f(mid,*args)<=0:
                high = mid 
            else:
                low = mid 
            mid = (low+high)/2 

    return mid


def point_set_match(Coord1,Coord2,r0,method):
    '''
    To find the point pairs (A,B) from two point sets Sa & Sb with some 
    method. The method can be one of the following

    'soulmate': (A,B) are called soulmates if
        1. A is the nearest neighbor of B among all the points in Sa;
        2. B is the nearest neighbor of A among all the points in Sb;
        3. The separation |A-B|<=r0.
        
        If two points are both the nearest neighbors, choose the one with 
        smaller index. Set r0 = +inf to eliminate the 3rd requirement.

    'friend': B is a friend of A if 
        1. B is the nearest neighbor of A among all the points in Sb;
        2. The separation |A-B|<=r0.

    Inputs
    ------
    Coord1: n1*k array, where n1 the point number of set Sa, k is the dimension
        of space. Cartesian coordinates of the points A in set Sa. E.g.
        [[x1,y1,z1,...], [x2,y2,z2,...], ...]
    Coord2: same as Coord1 but for set Sb
    r0: separation criterion 
    method: the method for pairing points, see above. When ='friend', 
        this function finds the friends of Coord1 among Coord2 

    Outputs
    -------
    Ind1: indexes of the paired points A in set Sa
    Ind2: indexes of the paired points B in set Sb
    '''

    Coord1 = np.expand_dims(Coord1,1) # n1*1*k
    Coord2 = np.array(Coord2) # n2*k
    n1 = len(Coord1)
    n2, k = Coord2.shape

    # distance matrix
    Dist = (Coord1-Coord2)**2 # n1*n2*k 
    Dist = np.sum(Dist,axis=-1)**.5 # n1*n2 

    # indexes of the nearest neighbor B of point A 
    Ind20 = np.argmin(Dist,axis=1) 

    # indexes of the nearest neighbor A of point B
    Ind10 = np.argmin(Dist,axis=0) 


    # indexes of point pairs
    Ind1, Ind2 = [], []

    if method=='friend':
        for i in range(n1):
            j = Ind20[i]
            if Dist[i,j]<=r0:
                Ind1.append(i)
                Ind2.append(j)

    elif method=='soulmate':
        for i in range(n1):
            for j in range(n2):
                if Ind20[i]==j and Ind10[j]==i and Dist[i,j]<=r0:
                    Ind1.append(i)
                    Ind2.append(j)

    Ind1 = np.array(Ind1) 
    Ind2 = np.array(Ind2) 


    return Ind1, Ind2


def test():
    #''' cubic_solver
    # parameters 
    N = 10000  # iteration number
    M = None  # array size, None for scalar

    t = []
    for i in range(N):
        a = np.random.uniform(-10, 10, M)
        b = np.random.uniform(-10, 10, M)
        c = np.random.uniform(-10, 10, M)
        d = np.random.uniform(-10, 10, M)
        tt = time.time()
        res = cubic_solver(a, b, c, d)
        t.append(time.time() - tt)
    print(f'Per iteration time: {1e6*np.mean(t)} µs.')
    #'''

    ''' dichotomy
    def f(x):
        return x 
    print(dichotomy(np.sin,2,np.pi+1e-6))
    #'''

    ''' point_set_match
    # parameters
    Coord1 = np.random.rand(10,2)
    Coord2 = np.random.rand(20,2)
    # Coord1 = np.array([[0,0],[.3,0],[.5,.7]])
    # Coord2 = np.array([[0,.5],[.2,.4]])
    r0 = 199

    # point pairs
    Ind1f, Ind2f = point_set_match(Coord1,Coord2,r0,'friend') 
    Ind1s, Ind2s = point_set_match(Coord1,Coord2,r0,'soulmate') 

    plt.figure(figsize=(7,6))
    plt.axis('equal')
    plt.scatter(Coord1[:,0],Coord1[:,1],s=10,label='Set Sa')
    plt.scatter(Coord2[:,0],Coord2[:,1],s=10,label='Set Sb')

    for i in range(len(Coord1)):
        plt.text(Coord1[i,0],Coord1[i,1],str(i),fontsize=14)
    for i in range(len(Coord2)):
        plt.text(Coord2[i,0],Coord2[i,1],str(i),fontsize=14)

    plt.plot([Coord1[:,0][Ind1f],Coord2[:,0][Ind2f]],
        [Coord1[:,1][Ind1f],Coord2[:,1][Ind2f]],color='b',lw=4,alpha=.3)
    plt.plot([Coord1[:,0][Ind1s],Coord2[:,0][Ind2s]],
        [Coord1[:,1][Ind1s],Coord2[:,1][Ind2s]],color='r',lw=1)

    plt.plot(0,0,color='b',lw=4,alpha=.3,label='friend')
    plt.plot(0,0,color='r',lw=1,label='soulmate')
    # plt.plot([0,r0],[0,0],color='k')
    plt.legend()
    plt.tight_layout()
    plt.show()
    #'''