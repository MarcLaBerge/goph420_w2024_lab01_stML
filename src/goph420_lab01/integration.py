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

    #Check that the shape and length of the new x,f arrays are the same shape
    if len(x.shape) != 1 or len(f.shape) != 1:
        raise ValueError ("Array is more than a 1D array")
    if len(x) != len(f):
        raise ValueError (f"The dimensions of x, {len(x)} and the dimensions of f, {len(f)} are not compatible")
    
    #Checking alg's entered value will be through a function
    N = len(x)
    integral = 0
    print(x[1] - x[0])
    #Start with trap
    if alg == "trap":
        #Possible error in the f[1] so if weird some to this. f[1:-1] starting at index 1 (not 0) 
        #and then adding everything until the last value
        for i in range (0,N-1):
            integral+=((x[i + 1] - x[i]) / 2) * (f[i + 1] + f[i])
        return integral

    #Simpson's rules
    elif alg == "simp":
        #If the number of intervals is odd -> all 1/3 rule
        if N % 2 !=  0:
            #np.sum(start(including):stop(excluding):every_#points)
            integral = ((x[1] - x[0]) / 3) * (f[0] + 4 * np.sum(f[1 : -4 : 2]) + 2 * np.sum(f[2 : -4 : 2] + f[N-4]))
            
        #If the number of intervals is even -> needs one round of 3/8 rule
        else:
            #If the number is greater than 4, we will use 1/3 rule, if it's 4 then we use 3/8 rule
            if N > 4:
                intgeral = ((x[1] - x[0]) / 3) * (f[0] + 4 * np.sum(f[1 : -4 : 2]) + 2 * np.sum(f[2 : -4 : 2] + f[N-4]))
            #Once N is 4, we will then add on 3/8 rule to finish
            integral += ((x[1] - x[0]) / 8) * (f[N-4] + 3 * f[N-3] + 3 * f[N-2] + f[N-1])
            

            #for i in range (0, N-4, 2):
                #integral += ((x[i+2]-x[i]) / 6) * (f[i]+ 4 * f[i+1] + f[i+2])
            #Simpson's 3/8 rule
                #I'm not sure if this is right
            #integral += (1 / 8) * ((x[N-1] - x[N-4]) / 3) * (f[N-4] + 3*f[N-3] + 3*f[N-2] + f[N-1])
            #return integral
        #If the amount of intervals is odd, needs more than 3 data points
        #elif N-1 % 2 != 0 and N >= 3:
            #for i in range (0, N-2, 2):
                #integral+=((x[i + 2] - x[i]) / 6) * (f[i] + 4 *f[i + 1]+f[i + 2])
        return integral
        
    #If alg is non of the above "trap" or "simp"-> raise ValueError
    else:
        raise ValueError ("Invalid algorithm entered, please enter 'simp' or trap'")


def integrate_gauss(f, lims, npts):
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
    try:
        callable(f)
    except ValueError:
        print("Object cannot be called")

    #Check that lims has a length of 2
    if len(lims) != 2:
        raise ValueError (f"The dimensions of the limit array is {len(lims)}, must be 2")
    
    #Check that npts is in [1,2,3,4,5]
    if npts not in [1,2,3,4,5]:
        raise ValueError (f"The number of points is not one of [1,2,3,4,5] please enter one of those numbers of points")

    #Change the lims array values into integer values to use in furture
    #Raise error if not possible
    try:
        a = int(lims[0])
        b = int(lims[1])
    except ValueError:
        print("The upper or lower limit bounds cannot be changed to integers")



    #Start the integral at 0
    integral = 0

    #Change the value of npts into an interger just incase
    npts = int(npts)
    

    #If the number of points = 1
    if npts == 1:
        weights = [2.0]
        s_k = [0.0]
    
    #If the number of points = 2
    elif npts == 2:
        weights = [1.0, 1.0]
        s_k = [-1/np.sqrt(3), 1/np.sqrt(3)]
    
    #If the number of points = 3
    elif npts == 3:
        weights = [5/9, 8/9, 5/9]
        s_k = [-1/np.sqrt(3/5), 0, np.sqrt(3/5)]

    #If the number of points = 4
    elif npts == 4:
        weights = [(18 - np.sqrt(30))/36, (18 + np.sqrt(30))/36, (18 + np.sqrt(30))/36, (18 - np.sqrt(30))/36]
        s_k = [(-1*np.sqrt((3/7) + (2/7)*np.sqrt(6/5))), -1*np.sqrt((3/7) - (2/7)*np.sqrt(6/5)), np.sqrt((3/7) - (2/7)*np.sqrt(6/5)), np.sqrt((3/7) + (2/7)*np.sqrt(6/5))]

    #If the number of points = 5
    elif npts == 5:
        weights = [(322 - 13*np.sqrt(70))/900, (322 + 13*np.sqrt(70))/900, 128/225, (322 + 13*np.sqrt(70))/900, (322 - 13*np.sqrt(70))/900]
        s_k = [-(1/3)*np.sqrt(5+2*np.sqrt(10/7)), -(1/3)*np.sqrt(5-2*np.sqrt(10/7)), 0, (1/3)*np.sqrt(5-2*np.sqrt(10/7)),(1/3)*np.sqrt(5+2*np.sqrt(10/7))]

    #Make an empty array for the new changed weights and points
    transpoints = np.zeros(npts)
    transweights = np.zeros(npts)

    #Finding the changes weights and s_k
    for i in range (0, npts):
        transweights[i] = (((b - a) / 2) * weights[i])
        transpoints[i] = (0.5 * (a + b) + 0.5 * (b - a) * s_k[i])

    #Calculating the integral, need to put the new points into the function
    for i in range (0, npts):
        integral += transweights[i] * f(transpoints[i])
    return integral


