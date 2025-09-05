import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import glob
import os
import re
import sys
from scipy.optimize import curve_fit

time_arr = []
peak_pos = []

# whether or not to plot the analytic wave function of coherent state
show_analytic_psi=True

def s(t, hbar, sigma_x):
    return 1.0+0.25*(hbar*t/(sigma_x*sigma_x))**2

def psi_real_imag(x, t, x0, v0, hbar, sigma_x):
    const = np.pow(2.0*np.pi*sigma_x**2, -0.25)
    fact = np.sqrt(1.0+0.5j*hbar*t/(sigma_x*sigma_x))
    sigma_t2 = sigma_x*sigma_x*s(t, hbar, sigma_x)
    exponent_real = -0.25*(x-x0-v0*t)**2/sigma_t2
    exponent_imag = (0.125*((x-x0)/(sigma_x*sigma_x))**2*hbar*t + v0*(x-x0)/hbar - (v0**2*t)/(2.0*hbar))/s(t, hbar, sigma_x)

    result = const*np.exp(exponent_real + 1j*exponent_imag)/fact

    return result.real, result.imag

def psi_squared_free(x, t, x0, v0, hbar, sigma_x):
    sigma_t = sigma_x*np.sqrt(1.0 + 0.25*(hbar*t)**2/(sigma_x**4))
    exponent = -0.5*((x-x0-v0*t)/sigma_t)**2

    return np.exp(exponent)/sigma_t/np.sqrt(2.0*np.pi)

def gaussian_func(x, A, mu,  sigma):
    exponent = -0.5*(x-mu)*(x-mu)/(sigma*sigma)

    return A*np.exp(exponent)

def read_header(filename):
    with open(filename, "r") as f:
        nmesh = int(f.readline().strip())
        hbar = float(f.readline().strip())
        tnow = float(f.readline().strip())
    
    return nmesh, hbar, tnow

def extract_index(filename):
    match = re.search(r'(\d+)(?=\.\w+$)', filename)
    return int(match.group(1)) if match else -1

def update(frame):
    nmesh, hbar, tnow = read_header(file_list[frame])
    data = np.loadtxt(file_list[frame], skiprows=3)
    x, y = data[:, 0], data[:, column]
    
    line.set_data(x, y)
    ax.relim()
    ax.autoscale_view()  # auto scale the vertical axis
    ax.set_title(f"$N$ = {nmesh} / $\hbar$ = {hbar:.2e} / $t$ = {tnow:.2f}",fontsize=10)

    return line,

def inspect_data_range(file_list):
    all_x = []
    all_y = []
    
    for f in file_list:
        data = np.loadtxt(f, skiprows=3)
        all_x.extend(data[:, 0])
        all_y.extend(data[:, column])
    
    xmin, xmax = min(all_x), max(all_x)
    ymin, ymax = min(all_y), max(all_y)

    xrange = xmax - xmin
    yrange = ymax - ymin
    
    margin_ratio = 0.05

    xlim = (xmin - margin_ratio * xrange, xmax + margin_ratio * xrange)
    ylim = (ymin - margin_ratio * yrange, ymax + margin_ratio * yrange)

    return xlim, ylim

def init():
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    return line,

if __name__  == "__main__":

    if len(sys.argv) != 3:
       print("Usage: python "+sys.argv[0]+" \"model_??.dat\" column")
       print("file list is given using the wildcard and must be double quated")
       sys.exit(1)
       
    model_name = sys.argv[1]

    file_list = sorted(glob.glob(sys.argv[1]), key=extract_index)

    column = int(sys.argv[2])
    xlim, ylim = inspect_data_range(file_list)


    fig = plt.figure(dpi=200)
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'

    ax = fig.add_axes((0.1, 0.1, 0.8, 0.8))
    ax.grid(axis='x',which='major', color='#e9e9e9')
    ax.grid(axis='y',which='major', color='#e9e9e9')

    line, = ax.plot([],[],'b-')

    ani = animation.FuncAnimation(
        fig, update, frames=len(file_list), init_func=init, blit=False, interval=500, repeat=False)

    #    plt.show()
    plt.close(fig)

    # eliminate extension
    name, _ = os.path.splitext(sys.argv[1])
    # eliminate wildcard characters
    prefix = re.sub(r"[?*]", "", name)
    
    for i, filename in enumerate(file_list):
        nmesh, hbar, tnow = read_header(filename)
        data = np.loadtxt(filename, skiprows=3)
        x = data[:,0]
        y = data[:,column]

        xx = np.linspace(-1,1,1024)
        psi_squared = psi_squared_free(xx, tnow, -0.5, 1.0, hbar, 0.05)
        psi_real, psi_imag = psi_real_imag(xx, tnow, -0.5, 1.0, hbar, 0.05)

        if column==1:
            fit_, cov_ = curve_fit(gaussian_func, x, y, p0=[1.0, x[np.argmax(y)], 0.2])
            amp_fit, peak_pos_fit, sigam_fit = fit_
            time_arr.append(tnow)
            peak_pos.append(peak_pos_fit)

        fig = plt.figure(dpi=200)
        plt.rcParams['xtick.direction'] = 'in'
        plt.rcParams['ytick.direction'] = 'in'
        
        ax = fig.add_axes((0.1, 0.1, 0.8, 0.8))

        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_xlabel("x")
        if column==1:
            ax.set_ylabel("$|\Psi(x)|^2$")
        if column==3:
            ax.set_ylabel("$Re[\Psi(x)]$")
        if column==4:
            ax.set_ylabel("$Im[\Psi(x)]$")
            
        ax.grid(axis='x',which='major', color='#e9e9e9')
        ax.grid(axis='y',which='major', color='#e9e9e9')

        ax.set_title(f"$N$ = {nmesh} / $\hbar$ = {hbar:.2e} / $t$ = {tnow:.2f}",fontsize=10)

        ax.scatter(x,y,marker="s", s=4);

        if show_analytic_psi:
            if column==1:
                ax.plot(xx, psi_squared, color="red",lw=0.5)
            if column==3:
                ax.plot(xx, psi_real, color="red", lw=0.5)
            if column==4:
                ax.plot(xx, psi_imag, color="red", lw=0.5)

        output_base, _= os.path.splitext(filename)
        output_filename = output_base+".png"

        plt.savefig(output_filename, bbox_inches="tight")
        plt.close(fig)

    if column==1:
        with open(prefix+"peak.dat", "w") as f:
            for i in range(len(time_arr)):
                print(time_arr[i], peak_pos[i], file=f)
