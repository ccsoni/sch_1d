import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import glob
import os
import re
import sys

column = 1

def extract_index(filename):
    match = re.search(r'(\d+)(?=\.\w+$)', filename)
    return int(match.group(1)) if match else -1

def update(frame):
    data = np.loadtxt(file_list[frame])
    x, y = data[:, 0], data[:, column]
    line.set_data(x, y)
    ax.relim()
    ax.autoscale_view()  # auto scale the vertical axis
    ax.set_title(f"Frame {frame} - {file_list[frame]}")
    return line,

def inspect_data_range(file_list):
    all_x = []
    all_y = []
    
    for f in file_list:
        data = np.loadtxt(f)
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
        fig, update, frames=len(file_list), init_func=init, blit=False, interval=500
    )

    plt.show()
    plt.close(fig)

    for i, filename in enumerate(file_list):
        data = np.loadtxt(filename)
        x = data[:,0]
        y = data[:,column]

        fig = plt.figure(dpi=200)
        plt.rcParams['xtick.direction'] = 'in'
        plt.rcParams['ytick.direction'] = 'in'
        
        ax = fig.add_axes((0.1, 0.1, 0.8, 0.8))

        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.grid(axis='x',which='major', color='#e9e9e9')
        ax.grid(axis='y',which='major', color='#e9e9e9')

        ax.plot(x,y);

        output_base, _= os.path.splitext(filename)
        output_filename = output_base+".png"

        plt.savefig(output_filename, bbox_inches="tight")
        plt.close(fig)


        
