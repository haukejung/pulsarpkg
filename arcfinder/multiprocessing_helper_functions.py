# these are all helper functions to help with multiprocessing on the
# automatic arc curvature detections.

import numpy as np


# returns the value at x of the 1D Gaussian distribution specified by mu and sigma.
def gaussian(x, mu, sig):
    return np.exp(-(x - mu)**2 / (2. * sig**2))


# user specifies a parabola starting at the origin with formula y=ax**2 , 
# and also a point with x=px and y=py, and this function returns
# the x-coordinate of the closest point on the parabola.
def closest_point_on_the_parabola(a,px,py):
    
    thingy1 = 2.*a*py
    thingy2 = np.sqrt(-3+0j)
    thingy3 = 2**(1/3.)
    thingy4 = 2**(2/3.)
    thingy = (-108.*a**4*px + np.sqrt(11664.*a**8*px**2 - 864.*a**6*(-1 + thingy1)**3 + 0j))**(1/3.)
    Aone = (thingy3*(-1. + thingy1))
    Atwo = thingy
    Athree = thingy
    Afour = (6.*thingy3*a**2)
    Bone = ((1. + thingy2)*(-1. + thingy1))
    Btwo = (thingy4*thingy)
    Bthree = ((1. - thingy2)*thingy)
    Bfour = (12.*thingy3*a**2)
    Cone = (1. - thingy2)*(-1 + thingy1)
    Ctwo = thingy4*thingy
    Cthree = (1. + thingy2)*thingy
    Cfour = 12.*thingy3*a**2
    
    A = -np.real(Aone/Atwo + Athree/Afour)
    B =  np.real(Bone/Btwo + Bthree/Bfour)
    C =  np.real(Cone/Ctwo + Cthree/Cfour)
    
    solns = [A,B,C]
    solns_temp = []
    for soln in solns:
        solns_temp.append(np.abs(soln-px))
    
    val, idx = min((val, idx) for (idx, val) in enumerate(solns_temp))
    return solns[idx]

# returns how far the specified point is from the closest point on the parabola.
def dist_from_parabola(a,px,py):
    X = closest_point_on_the_parabola(a,px,py)
    return (X,np.sqrt(a**2*X**4-2.*a*X**2*py+py**2+X**2-2.*px*X+px**2))


def weight_function(eta,px,py,sigma=[1,1]):
    x = closest_point_on_the_parabola(eta,px,py)
    y = eta*x**2
    dist_from_origin = np.sqrt(px**2+py**2)
    ret_y = gaussian(py-y,0,sigma[0])
    ret_x = gaussian(px-x,0,sigma[1])
    ret = dist_from_origin*np.sqrt(ret_x*ret_y)
    if ret < (1/np.e**3):
        ret = None
    return ret

def weight_function2(y,x,py,px,sigma):
    ret_y = gaussian(py-y,0,sigma[0])
    ret_x = gaussian(px-x,0,sigma[1])
    ret = np.sqrt(ret_x*ret_y)
    if ret < (1/np.e**3):
        ret = None
    return ret

def weight_function3(eta,orig_y,orig_x,py,px,sigma):
    x = closest_point_on_the_parabola(eta,px,py)
    y = eta*x**2
    dist_from_origin = np.sqrt(orig_x**2+orig_y**2)
    ret_y = gaussian(py-y,0,sigma[0])
    ret_x = gaussian(px-x,0,sigma[1])
    ret = dist_from_origin*np.sqrt(ret_x*ret_y)
    if ret < (1/np.e**3):
        ret = None
    return ret

def crunchy(eta,sec,hand=None,sigma=None):
    powers = []
    powers_norm = []

    y_axis = sec.get_y_axis()
    x_axis = sec.get_x_axis()
    if sigma==None:
        sigma = [np.absolute(y_axis[1]-y_axis[0]),
                 np.absolute(x_axis[1]-x_axis[0])]

    for yi in range(len(y_axis)):
        y = y_axis[yi]
        for xi in range(len(x_axis)):
            x = x_axis[xi]
            this_weight = weight_function(eta,x,y,sigma)
            if this_weight is None:
                powers.append(None)
                powers_norm.append(None)
            else:
                variance = 1/this_weight
                powers.append(sec.get([yi,xi])/variance)
                powers_norm.append(1/variance)
    p = np.nansum(list(filter(None,powers)))
    pn = np.nansum(list(filter(None,powers_norm)))
    #print("eta: " + str(eta))
    #print(p)
    #print(pn)
    #print("p/pn: " + str(p/pn))
    return (eta,p/pn)

def crunchy2(pt_and_sigma,sec,hand=None):
    pt,sigma = pt_and_sigma
    py,px = pt

    powers = []
    powers_norm = []

    y_axis = sec.get_y_axis()
    x_axis = sec.get_x_axis()
    px_y = np.absolute(y_axis[1]-y_axis[0])
    px_x = np.absolute(x_axis[1]-x_axis[0])

    if sigma==None:
        sigma = [px_y,px_x]
    if sigma[0]<px_y:
        sigma = [px_y,sigma[1]]
    if sigma[1]<px_x:
        sigma = [sigma[0],px_x]

    for yi in range(len(y_axis)):
        y = y_axis[yi]
        for xi in range(len(x_axis)):
            x = x_axis[xi]
            this_weight = weight_function2(y,x,py,px,sigma)
            if this_weight is None:
                powers.append(None)
                powers_norm.append(None)
            else:
                variance = 1/this_weight
                powers.append(sec.get([yi,xi])/variance)
                powers_norm.append(1/variance)
    p = np.nansum(list(filter(None,powers)))
    pn = np.nansum(list(filter(None,powers_norm)))
    return (pt,p/pn)

def crunchy3(offset, eta, sec, sigma=None):

    powers = []
    powers_norm = []

    y_axis = sec.get_y_axis()
    x_axis = sec.get_x_axis()
    px_y = np.absolute(y_axis[1]-y_axis[0])
    px_x = np.absolute(x_axis[1]-x_axis[0])

    if sigma==None:
        sigma = [px_y,px_x]
    if sigma[0]<px_y:
        sigma = [px_y,sigma[1]]
    if sigma[1]<px_x:
        sigma = [sigma[0],px_x]

    for yi in range(len(y_axis)):
        y = y_axis[yi]
        for xi in range(len(x_axis)):
            x = x_axis[xi]
            y_eff = y + eta*offset**2
            x_eff = x - offset
            this_weight = weight_function3(eta,y,x,y_eff,x_eff,sigma)
            if this_weight is None:
                powers.append(None)
                powers_norm.append(None)
            else:
                variance = 1/this_weight
                powers.append(sec.get([yi,xi])/variance)
                powers_norm.append(1/variance)
    p = np.nansum(list(filter(None,powers)))
    pn = np.nansum(list(filter(None,powers_norm)))
    return (offset,p/pn)
