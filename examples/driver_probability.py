#Will contain code that generates analysis results
#and generates any plots for Section 5.2
import numpy as np
import matplotlib.pyplot as plt
from goph420_lab01.integration import integrate_gauss

def main():
    #initialize given values, standard deviation and the annual mean
    std_dev = 0.5
    mean = 1.5

    #Function from equation 17
    f = lambda x: (1/np.sqrt(2 * np.pi)) * (np.exp((-1 / 2)*((x - mean) / (std_dev) ** 2)))

    #Number of points
    npts = [1,2,3,4,5]

    #Creating arrays for the convergent integrals as well as the relative errors
    converged_integrals = []
    relative_error = []

    #Initialization past integral
    past_integral = 0

    #Calculate the probablitiy of a M > 4.0 event, lims upper bound is 4 while the mean of 0.5 is the lower bound
    for i in npts:
        new_integral = integrate_gauss(f, lims = [1.5, 4], npts = i)
        eps_a = np.abs((new_integral - past_integral) / new_integral)
        #Add the convergence of the integral to the array
        converged_integrals.append(new_integral)
        #Add the relative error to the empty array
        relative_error.append(eps_a)
        #Replace the old with the new integral, then move on to the next
        past_integral = new_integral

    plt.loglog(npts, relative_error, 'teal')
    plt.ylabel("Relative Error")
    plt.xlabel("Number of points used in the integration")
    plt.title("Probability of a Seimic Event with a magnitude greater than 4.0")
    plt.savefig("data/probability_seismic_event_convergence.png")
    plt.close("all")

    #-----------------------------------------------------------------------------
    #initialize given values, standard deviation and the annual mean
    std_dev = 0.05
    mean = 10.28

    #Function from equation 17
    f = lambda x: (1/np.sqrt(2 * np.pi)) * (np.exp((-1 / 2)*((x - mean) / (std_dev) ** 2)))
    
    #Creating arrays for the convergent integrals as well as the relative errors
    convergent_integrals = []
    relative_error = []

    #Initialization past integral
    past_integral = 0

    #Calculate the probablitiy of a M > 4.0 event, lims upper bound is 4 while the mean of 0.5 is the lower bound
    for i in npts:
        new_integral = integrate_gauss(f, lims = [10.25, 10.35], npts = i)
        eps_a = np.abs((new_integral - past_integral) / new_integral)
        #Add the convergence of the integral to the array
        convergent_integrals.append(new_integral)
        #Add the relative error to the empty array
        relative_error.append(eps_a)
        #Replace the old with the new integral, then move on to the next
        past_integral = new_integral

    plt.loglog(npts, relative_error, 'teal')
    plt.ylabel("Relative Error")
    plt.xlabel("Number of points used in the integration")
    plt.title("Probability that 10.25m <= Length <= 10.35m")
    plt.savefig("data/probability_length_interval_convergence.png")
    plt.close("all")

if __name__ == "__main__":
    main()