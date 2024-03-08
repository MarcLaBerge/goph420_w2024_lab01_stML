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
    
    #Make x and t arrays
    x = np.array(x)
    f = np.array(f)

    #Check that the shape and length of the new x,f arrays are the same shape
    if len(x.shape) != 1 or len(f.shape) != 1:
        raise ValueError ("Array is more than a 1D array")
    if len(x) != len(f):
        raise ValueError (f"The dimensions of x, {len(x)} and the dimensions of f, {len(f)} are not compatible")
    
    #Checking alg's entered value will be through a function
    N = len(x)
    integral = 0
    #Start with trap
    if alg == "trap":
        #Possible error in the f[1] so if weird some to this
        integral = (x[1] - x[0])/2 + (f[0] + 2*np.sum(f[1:-1]) + f[1])

        return integral

    #Simpson's rules
    elif alg == "simp":
        #If the amount of intervals is even, need more than 3 data points
        if N-1 % 2 == 0 and N >= 3:
            #Simpson's 1/3 rule
            for i in range (0,N-4,2):
                integral += ((x[N+2]-x[N]) / 6) * (f[N]+ 4 * f[N+1] + f[N+2])
            #Simpson's 3/8 rule
                #I'm not sure if this is right
            integral += (1 / 8) * ((x[N-1] - x[N-4]) / 3) * (f[N-4] + 3*f[N-3] + 3*f[N-2] + f[N-1])
            return integral
        #If the amount of intervals is odd, needs more than 3 data points
        elif N-1 % 2 != 0 and N >= 3:
            integral += (1 / 3)*((x[i+2] - x[i]) / 2) * (f[i]+ 4 * f[i+1] + f[i+2])
            return integral
        
    #If alg is non of the above -> raise ValueError
    else:
        raise ValueError ("Invalid algorithm entered, please enter 'simp' or trap'")



    

   

      
    
    return

def integrate_gauss(f, lims, npts):
    """
    Integration from sample points using Gauss-Legendre quadrature

    Parameters:
    -----------
        f:
            Reference to a callable object
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
    pass
alg = 'str'
integrate_newton(1,1,alg)


