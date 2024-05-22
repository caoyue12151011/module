import os
import numpy as np 
import matplotlib.pyplot as plt 

# get path of this script
path = os.path.dirname(os.path.abspath(__file__))

def pdf2points(x, pdf, N):
    '''
    To generate random points given probability density function (PDF).
    The PDF is regarded as a continuous function (see below). 
    If you are dealing with discrete PDF (i.e. PMF), use np.random.choice
    instead. 

    Inputs
    ------
    x: 1darray, x-axis sample of the PDF
    pdf: 1darray, must > 0 and without nans, can be non-normalized
    N: scalar, number of points to be generated

    Returns
    -------
    x_sample: scalar or 1darray.
    '''
    # create CDF
    cdf = np.cumsum(pdf*np.gradient(x)) 
    cdf = (cdf-cdf[0]) / (cdf[-1]-cdf[0])

    # generate random points
    x_sample = np.interp(np.random.rand(N), cdf, x)
    if N==1:
        x_sample = x_sample[0]
    return x_sample


def test():
    x = np.linspace(-2,2,1000)
    pdf = np.exp(-x**2/(2*.7**2)) + 5*np.exp(-(x-1)**2/(2*.1**2))
    N = 5000

    pdf /= np.trapz(pdf, x)
    x_sample = pdf2points(x, pdf, N)

    plt.figure()
    plt.plot(x, pdf, color='k', lw=1)
    plt.hist(x_sample, histtype='stepfilled', bins=int(N**.5), density=True,
             fc='r', ec='none')
    plt.grid()
    plt.xlabel('x')
    plt.ylabel('Probability density')
    plt.tight_layout()
    plt.savefig(f'{path}/image/pdf2points.pdf')
    plt.close()