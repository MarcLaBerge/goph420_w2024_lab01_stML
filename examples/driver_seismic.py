#Will contain code that generates analysis results
#and generates any plots for section 5.1
import numpy as np
import matplotlib.pyplot as plt
from goph420_lab01.integration import integrate_newton

def main():
    #Loading velocity data from file
    wave_data = np.loadtxt('data/s_wave_data.txt')
    #Seperate data into time and velocity
    t = wave_data[:,0]
    v = wave_data[:,1]
    #X axis will be (t)ime (seconds), Y axis will be (v)elocity (mm/s)

    
    #Finding the end of the period of time that is reasonable
    v_abs_max = np.max(np.abs(v))
    v_end = 0.005*v_abs_max
    #i = np.where(np.abs(v) > v_end)
    #int_limit = int(max(t[i]))
    for i, _ in enumerate(v):
        if np.abs(v[i]) > v_end:
            int_limit = i
    
    
    #Plotting the raw data from s_wave_data.txt
    plt.plot(t, v, 'r-', label = 'S-Wave Arrivals', linewidth = 0.5)
    plt.vlines(x=[0], ymin=[-0.3], ymax=[0.3],colors='green',ls='--',label = 'Period start time')
    plt.vlines(x=t[int_limit], ymin=[-0.3], ymax=[0.3], colors='teal', ls='--', label = 'Period end time')
    plt.ylabel('Velocity (v) [mm/s]')
    plt.xlabel('Time (t) [second]')
    plt.title('Collected data for S-Wave Arrivals')
    plt.legend()
    plt.savefig('figures/s_wave_data_figure.png')
    plt.close('all')

    #---------------------------------------------------------------------------
    #integration using trap, simpson's 1/3 and 3/8 rules
    #   determine the period
    T = t[int_limit]
    #Slice the t and v values to stop at the int_limit (So that integration won't go past)
    t = t[:int_limit]
    v = v[:int_limit]

    #Estimating the integral with different sampling intervals, I think that bigger intervals won't have enough data points
    intervals = [1,2,4,8,16,32]
    #Stepsize 
    delta_a = [0.01, 0.02, 0.04,0.08,0.16,0.32]

    #Create empty array for the integrated values of the rules to be added to
    int_trap = []
    int_simp = []
    

    #Implementing the functions with imputs through an interval step size as well as the kind of algorithm to use
    for step in intervals:
        #Entering the required values for the newton function, ::step, the entire sliced array, going through it at intervals
        integrated_trap = integrate_newton(t[::step], v[::step], alg = "trap") #Starting with the trap rule
        integrated_simp = integrate_newton(t[::step], v[::step], alg = "simp") #Simpson's rule

        #Add the integration values into the empty arrays from earlier
        int_trap.append(integrated_trap)
        int_simp.append(integrated_simp)

    #Need to find the relative error for figure
        #From Equation 15
    eps_trap = np.abs(np.diff(int_trap)) / int_trap[:-1]
    eps_simp = np.abs(np.diff(int_simp)) / int_simp[:-1]


    #Plot the figure -> curve of the convergence of each integration rule
    plt.loglog(delta_a[1:], eps_trap, label = 'Trapezoid rule')
    plt.loglog(delta_a[1:], eps_simp, label = "Simpson's rule")
    plt.ylabel("Approximate Relative Error [eps_s]")
    plt.xlabel("Sampling Interval [deltax]")
    plt.legend()
    plt.savefig("figures/trap_simp_convergence.png")
    plt.close("all")
    
    #------------------------------------------------------------------------------


if __name__ == '__main__':
    main()
