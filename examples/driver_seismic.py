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

    plt.plot(t, v, 'k-', label = 'S-Wave Arrivals', linewidth = 0.5)
    plt.ylabel('Velocity (v) [mm/s]')
    plt.xlabel('Time (t) [second]')
    plt.title('Collected data for S-Wave Arrivals')
    plt.legend()
    plt.savefig('data/s_wave_data_figure.png')
    plt.close('all')

if __name__ == '__main__':
    main()
