

def integration_newton(x,f,alg):
    """
    Integration of sample points using Newton-Cotes rules
    Parameters:
    -----------
        x:
            Array like, same shape containing coordinatss and values of the sample points as f
        f:
            Array like, same shape containing coordinates and values of the sample points as x
        alg:
            (Optional) str, default value of "trap", other option is "simp", deciding which inetrgration rule to use

    Raises:
    -------
        ValueError:
            If alg contains a str othan than "trap" and "str" 
        TypeError:
            If alg is not a str
        ValueError:
            If the dimensions of x and f are incompatible
                Too many dimensions and/or different length

    Return:
    -------
        Integral estimate:
            Float

    """





