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
    ax.set_title(f"Frame {frame} - {file_list[frame]}")
    ax.text(0.05, 0.95,
            f"$N$ = {nmesh}\n$\hbar$ = {hbar:.2e}\n$t$ = {tnow:.2f}",
            transform=ax.transAxes,
            va="top", ha="left",
            fontsize=12, bbox=dict(boxstyle="round", facecolor="white"))
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

#    file_list = sorted(glob.glob('model_??.dat'), key=extract_index)
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

        if column==1:
            fit_, cov_ = curve_fit(gaussian_func, x, y, p0=[1.0, x[np.argmax(y)], 0.2])
            amp_fit, peak_pos_fit, sigam_fit = fit_
            #        print(tnow, peak_pos_fit)
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

        ax.text(0.05, 0.95,
                f"$N$ = {nmesh}\n$\hbar$ = {hbar:.2e}\n$t$ = {tnow:.2f}",
                transform=ax.transAxes,
                va="top", ha="left",
                fontsize=12, bbox=dict(boxstyle="round", facecolor="white"))

        ax.scatter(x,y,marker="s",s=16);

        output_base, _= os.path.splitext(filename)
        output_filename = output_base+".png"

        plt.savefig(output_filename, bbox_inches="tight")
        plt.close(fig)

    if column==1:
        with open(prefix+"peak.dat", "w") as f:
            for i in range(len(time_arr)):
                print(time_arr[i], peak_pos[i], file=f)
