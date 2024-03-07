#Will contain code that generates analysis results
#and generates any plots for section 5.1
import numpy as np
import matplotlib.pyplot as plt

def main():
    #Loading velocity data from file
    wave_data = np.loadtxt('data/s_wave_data.txt')
    #Seperate data into time and velocity
    t = wave_data[:,0]
    v = wave_data[:,1]
    #X axis will be (t)ime (seconds), Y axis will be (v)elocity (mm/s)

    
    #Finding the end of the period of time that is reasonable
    v_abs_max = max(abs(v))
    v_min = 0.005*v_abs_max
    i = np.where(v > v_min)
    end_time = max(t[i])
    
    #print("Largest velocity",v_abs_max)
    #print("Period cut off velocity",end_time)
    
    

    plt.plot(t, v, 'r-', label = 'S-Wave Arrivals', linewidth = 0.5)
    plt.vlines(x=[0], ymin=[-0.3], ymax=[0.3],colors='green',ls='--',label = 'Period start time')
    plt.vlines(x=[end_time], ymin=[-0.3], ymax=[0.3], colors='teal', ls='--', label = 'Period end time')
    plt.ylabel('Velocity (v) [mm/s]')
    plt.xlabel('Time (t) [second]')
    plt.title('Collected data for S-Wave Arrivals')
    plt.legend()
    plt.savefig('data/s_wave_data_figure.png')
    plt.close('all')

if __name__ == '__main__':
    main()
