import numpy as np

def integrate_newton(x,f,alg = "trap"):
    """
    Integration of sample points using Newton-Cotes rules

    Parameters:
    -----------
        x:
            Array like, same shape containing coordinates and values of the sample points as f
        f:
            Array like, same shape containing coordinates and values of the sample points as x
        alg:
            (Optional) str, default value of "trap", other option is "simp", deciding which inetrgration rule to use

    Raises:
    -------
        ValueError:
            If alg contains a str othan than "trap" and "simp" 
        TypeError:
            If alg is not a str, python will check
        ValueError:
            If the dimensions of x and f are incompatible
                Too many dimensions and/or different length

    Return:
    -------
        Integral estimate:
            Float

    """
    #Make sure alg is not case sensitive
    alg = alg.strip().lower()
    
    #Make x and t arrays or floats
    x = np.array(x, dtype = float)
    f = np.array(f, dtype = float)

    #Check that the shape and length of new x,f arrays are the same shape
    if len(x.shape) != 1 or len(f.shape) != 1:
        raise ValueError ("Array is more than a 1D array")
    if len(x) != len(f):
        raise ValueError (f"The dimensions of x, {len(x)} and the dimensions of f, {len(f)} are not compatible")
    
    #Checking alg's entered value will be through a function
    N = len(x)
    integral = 0.0
    print(x[1] - x[0])
#------------------------------------
    #Trap rule
    if alg == "trap":
        #np.sum(start(including):stop(excluding):every_#points)
        integral = ((x[1] - x[0]) / 2) * (f[0] + 2 * np.sum(f[1:-1]) + f[-1])
        
#---------------------------------------------
    #Simpson's rules
    elif alg == "simp":
        #If the number of intervals is odd -> all 1/3 rule
        if N % 2:
            #np.sum(start(including):stop(excluding):every_#points)
            integral = ((x[1] - x[0]) / 3) * (f[0] + 4 * np.sum(f[1 : -1 : 2]) + 2 * np.sum(f[2 : -1 : 2]) + f[-1])
        
        #If the number of intervals is even -> needs one round of 3/8 rule
        else:
            #If the number is greater than 4, we will use 1/3 rule, if it's 4 then we use 3/8 rule
            if N > 4:
                integral = ((x[1] - x[0]) / 3) * (f[0] + 4 * np.sum(f[1 : -1 : 2]) + 2 * np.sum(f[2 : -1 : 2]) + f[-1])
            #Once N is 4, we will then add on 3/8 rule to finish
            integral += ((x[1] - x[0]) / 8) * (f[-4] + 3 * f[-3] + 3 * f[-2] + f[-1])
        
    #If alg is non of the above "trap" or "simp"-> raise ValueError
    else:
        raise ValueError ("Invalid algorithm entered, please enter 'simp' or trap'")
    
    return integral

def integrate_gauss(f, lims, npts = 3):
    """
    Integration from sample points using Gauss-Legendre quadrature

    Parameters:
    -----------
        f:
            Reference to a callable object, function or class that implements the __call__() method
        lims:
            Object with a len of 2 that contains the upper and lower bounds (x=a and x=b)
        npts:
            (Optional) Int that gives the number of integration points to use
                Possible values are 1,2,3,4,5 default will be 3

    Raises:
    -------
        TypeError:
            If f isn't callable
        ValueError:
            If lims doesn't have a len == 2
        ValueError:
            If lims[0] or lims[1] are not convertable to float
        ValueError:
            If npts is not [1,2,3,4,5]
    
    Return:
    -------
        Integral estimate:
            Float
    """
    #Check that f is callable
    if not callable(f):
        raise TypeError(f"The function {f} is not callable")

    #Check that lims has a length of 2
    if len(lims) != 2:
        raise ValueError (f"The dimensions of the limit array is {len(lims)}, must be 2")
    
    #Check that npts is in [1,2,3,4,5]
    if npts not in [1,2,3,4,5]:
        raise ValueError (f"The number of points is not one of [1,2,3,4,5] please enter one of those numbers of points")

    #Change the lims array values into integer values to use in furture
    #Raise error if not possible
    try:
        a = float(lims[0])
        b = float(lims[1])
    except ValueError:
        print("The upper or lower limit bounds cannot be changed to integers")
    
    #If the number of points = 1
    if npts == 1:
        s_k = [0.0]
        c_k = [2.0]
    
    #If the number of points = 2
    elif npts == 2:
        s_k = [- 1 / np.sqrt(3), 1 / np.sqrt(3)]
        c_k = [1.0, 1.0]
    
    #If the number of points = 3
    elif npts == 3:
        s_k = [- np.sqrt(3/5), 0.0, np.sqrt(3/5)]
        c_k = [5/9, 8/9, 5/9]

    #If the number of points = 4
    elif npts == 4:
        s_k = [(-1*np.sqrt((3/7) + (2/7)*np.sqrt(6/5))), -1*np.sqrt((3/7) - (2/7)*np.sqrt(6/5)), np.sqrt((3/7) - (2/7)*np.sqrt(6/5)), np.sqrt((3/7) + (2/7)*np.sqrt(6/5))]
        c_k = [(18 - np.sqrt(30))/36, (18 + np.sqrt(30))/36, (18 + np.sqrt(30))/36, (18 - np.sqrt(30))/36]
    #If the number of points = 5
    elif npts == 5:
        s_k = [-(1/3)*np.sqrt(5+2*np.sqrt(10/7)), -(1/3)*np.sqrt(5-2*np.sqrt(10/7)), 0, (1/3)*np.sqrt(5-2*np.sqrt(10/7)),(1/3)*np.sqrt(5+2*np.sqrt(10/7))]
        c_k = [(322 - 13*np.sqrt(70))/900, (322 + 13*np.sqrt(70))/900, 128/225, (322 + 13*np.sqrt(70))/900, (322 - 13*np.sqrt(70))/900]

    #Make an empty array for the new changed/transfromed weights and points
    transpoints = [0]
    transweights = [0]
    
    #Put the transformed points into the function into an array for the function
    for i, j in enumerate(s_k):
        tran_s_k = 0.5*(a + b)+ 0.5 * (b - a) * s_k[i]
        f_k = f(tran_s_k)
        transpoints.append(f_k)
        print(transpoints)
    
    for i, j in enumerate(c_k):
        wk = 0.5* (b - a) * c_k[i]
        transweights.append(wk)
        print(transweights)

    #Calculating the integral (sum of the product of the transformed weights and the function of the transformed points)
    integral = np.sum(np.array(transweights) * np.array(transpoints))
    return integral

