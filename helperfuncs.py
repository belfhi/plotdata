import numpy as np

def search(arr, t, eps=0.05, forward=True):
    tii = np.where((arr < t+eps) & (arr > t-eps))[0]
    if len(tii) > 0: 
        if forward:
            return tii[0]
        else:
            return tii[-1]
    else:
        return #search(arr, t, eps=2*eps)
twopi = 2*np.pi

def search_indx(arr, t, eps=0.01):
    try:
        tii = np.where( (arr<t+eps) & (arr>t-eps))[0][0]
        return tii
    except IndexError:
        return

def gaussian(x, a,b,c):
    """
    Simple function that returns the gaussian
    of an array x
    """
    return a*np.exp((-(x-b)**2)/c)

def powerlaw(x, a, b, x0=0):
    """
    power-law with exponent b and (optional)
    shift of x0
    """
    return a*(x - x0)**b

def linear(x, a, b):
    """
    linear function a*x + b
    """
    return a + b*x

def parabola(x, a, b, c):
    """
    parabolo function
    """
    return a*(x-b)**2 - c

clrs = ['seagreen','red','blue','gray','purple','green','orange','dodgerblue','tomato','black']
