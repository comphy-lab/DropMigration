# Author: Vatsal Sanjay
# vatsalsanjay@gmail.com
# Physics of Fluids
# Last updated: 22-Oct-2023

import numpy as np
import os
import subprocess as sp
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.collections import LineCollection
from matplotlib.ticker import StrMethodFormatter
import pandas as pd
import concurrent.futures
import sys
import argparse


matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['text.usetex'] = True

def execute_process(exe):
    p = sp.Popen(exe, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout, stderr = p.communicate()
    return stderr.decode("utf-8").split("\n")

def get_segs(place):
    temp2 = execute_process(["./getFacets", place])
    segs = []
    skip = False
    temp2 = list(filter(None, temp2))
    for n1 in range(len(temp2)):
        temp3 = temp2[n1].split(" ")
        if not skip:
            temp4 = temp2[n1+1].split(" ")
            x1, y1 = map(float, [temp3[0], temp3[1]])
            x2, y2 = map(float, [temp4[0], temp4[1]])
            segs.extend([((x1, y1),(x2, y2))])
            skip = True
        else:
            skip = False
    return segs

def get_field_values(place, xmin, xmax, ymin, ymax, ny):
    temp2 = list(filter(None, execute_process(["./getData", place, str(xmin), str(ymin), str(xmax), str(ymax), str(ny)])))
    data = np.array([line.split() for line in temp2], dtype=float)
    nx = data.shape[0] // ny
    X = data[:,0].reshape((nx, ny)).transpose()
    Y = data[:,1].reshape((nx, ny)).transpose()
    T = data[:,2].reshape((nx, ny)).transpose()
    return X, Y, T, nx

def plot_graphics(t, name, xmin, xmax, ymin, ymax, segs, T):

    fig, ax = plt.subplots()
    fig.set_size_inches(19.20, 10.80)
    rect = matplotlib.patches.Rectangle((xmin, ymin), xmax-xmin, ymax-ymin, linewidth=2, edgecolor='k', facecolor='none')
    ax.add_patch(rect)
    line_segments = LineCollection(segs, linewidths=4, colors='green', linestyle='solid')
    ax.add_collection(line_segments)

    ax.imshow(T, cmap="coolwarm", interpolation='Bilinear', origin='lower', extent= [xmin, xmax, ymin, ymax], vmax = 4.0, vmin = 0.0)
        
    ax.set_aspect('equal')
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)

    ax.set_title(r"$t = %5.4f$" % t, fontsize=20)
    ax.axis('off')        

    plt.savefig(name, bbox_inches='tight', dpi=300)
    plt.close()

def process_file(ti, folder, tSnap, xmin, xmax, ymin, ymax, ny, lw=2):
    t = tSnap*ti
    place = f"intermediate/snapshot-{t:5.4f}"
    name = f"{folder}/{int(t*1000):08d}.png"

    if not os.path.exists(place):
        print(f"{place} File not found!")
        return
    elif os.path.exists(name):
        print(f"{name} Image present!")
        return

    segs = get_segs(place)
    X, Y, T, nz = get_field_values(place, xmin, xmax, ymin, ymax, ny)
    
    plot_graphics(t, name, xmin, xmax, ymin, ymax, segs, T)

    return nz


def main():
    # Set up command-line argument parser
    parser = argparse.ArgumentParser(description='Make videos of standing waves')
    parser.add_argument('--num_workers', type=int, default=8, help='Number of workers')
    parser.add_argument('--tSnap', type=float, default=0.1, help='Snapshot time interval')
    parser.add_argument('--L0', type=float, default=1., help='Length of the domain')

    args = parser.parse_args()

    nGFS = 300
    num_workers = args.num_workers
    tSnap = args.tSnap
    L0 = args.L0

    ny = 512
    xmin, xmax, ymin, ymax = -L0/2., L0/2., -L0/2., L0/2.
    lw = 2
    folder = 'Video'
    os.makedirs(folder,exist_ok=True)

    with concurrent.futures.ProcessPoolExecutor(max_workers=num_workers) as executor:
        futures = [executor.submit(process_file, ti, folder, tSnap, xmin, xmax, ymin, ymax, ny, lw=lw) for ti in range(nGFS)]
        for future in concurrent.futures.as_completed(futures):
            exception = future.exception()
            if exception is not None:
                print(f"Exception occurred: {exception}")
            else:
                nx = future.result()
                print(f"Done {futures.index(future)+1}/{nGFS}! with nx = {nx}")

if __name__ == "__main__":
    main()
