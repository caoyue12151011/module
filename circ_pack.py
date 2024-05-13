from scipy.stats import loguniform 
import numpy as np 
import matplotlib.pyplot as plt 



def circ_pack(S,eta,f):
    '''
    To pack circles in a rectangle without overlapping.

    Inputs
    ------
    S: 1D array, areas of the circles
    eta: width/height ratio of the rectangle in which the circles are packed
    f: filing factor, =total circle area/rectangle area

    Returns
    -------
    Xc, Yc: 1D array, coordinates of the circle centers. The rectangle is placed
        at (0,0) and (w,h).
    w: width of the rectangle
    h: height of the rectangle
    '''


    N = len(S) # number of circles
    R = (S/np.pi)**.5 # circle radii


    # find size of the rectangle
    h = (np.pi/eta/f*np.sum(R**2))**.5
    w = eta*h 


    # sort circles from small to big
    Ind = np.argsort(R)
    Ind_r = np.argsort(Ind) # for undo sorting
    R_s = R[Ind]


    # randomly place circles from big to small
    Xc = np.array([])
    Yc = np.array([])
    for i in range(N):
        r = R_s[N-1-i]

        while True:
            xc = np.random.uniform(r,w-r)
            yc = np.random.uniform(r,h-r)

            # check overlapping 
            if np.all(((Xc-xc)**2+(Yc-yc)**2)**.5 > r+R_s[N-1:N-len(Xc)-1:-1]):
                Xc = np.append(Xc,xc)
                Yc = np.append(Yc,yc)
                break

        if i%int(N/100)==0:
            print('Circle packing %.1f%% done.'%(100*i/N))

    # undo sorting
    Xc = Xc[::-1][Ind_r]
    Yc = Yc[::-1][Ind_r]

    return Xc, Yc, w, h 



def test():

    # parameters
    S = loguniform.rvs(.01,10,size=3000)
    eta = 1.5 
    f = .5
    figx = 8 # [in]

    Xc, Yc, R_s, w, h = circ_pack(S,eta,f)


    plt.figure(figsize=(figx,figx/eta))
    ax = plt.gca()
    plt.axis('equal')
    plt.plot([0,w,w,0,0],[0,0,h,h,0],color='k')
    for i in range(len(S)):
        circle = plt.Circle((Xc[i],Yc[i]),R_s[i],facecolor='r',edgecolor='none',
            alpha=.5)
        ax.add_artist(circle)
    # plt.xticks([])
    # plt.yticks([])
    plt.tight_layout()
    plt.show()











