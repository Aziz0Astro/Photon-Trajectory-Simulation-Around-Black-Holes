#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 28 19:59:46 2024
@author: azizabdulaziz
"""

# Import necessary libraries
import numpy as np
import matplotlib.pyplot as plt
from coloredLine import colored_line_between_pts  # Custom function for colored lines
from SingleRayPath import RayPathUntilAccDisk  # Custom function to compute ray paths
from scipy.interpolate import make_interp_spline  # For smooth curve fitting
from matplotlib.collections import LineCollection  # For efficient line plotting

# Function to find the index of the first negative number in an array
def giveFirstNegativeNumberIndex(arr):
    index = 0
    for i in range(len(arr)):
        if arr[i] < 0:
            index = i
            break
    return index

# Function to compute multiple ray paths around a Schwarzschild black hole
def multipleRayPathsAroundBH(r0, nrays, M=1., phi0=0, t0=0):
    """
    Computes photon paths around a Schwarzschild black hole, given initial
    conditions. The paths are plotted with time-dependent color coding.
    
    Parameters:
        r0: Initial radial distance
        nrays: Number of rays to simulate
        M: Black hole mass (default: 1)
        phi0: Initial angular position (default: 0)
        t0: Initial time (default: 0)
    
    Returns:
        hitDiskY: Array of Y-coordinates where rays hit the accretion disk
        hitDiskTimes: Array of times at which rays hit the accretion disk
    """
    # Compute the initial angle psi based on r0
    if np.isclose(r0, 3 * M):
        psi = np.pi / 2.
    elif r0 < 3 * M:
        psi = abs(np.arctan(np.sqrt((r0 / (2. * M) - 1) / (r0 / (6. * M) + 1.)) * (1 / (r0 / (3. * M) - 1.))))
    elif r0 > 3 * M:
        psi = np.pi - abs(np.arctan(np.sqrt((r0 / (2. * M) - 1) / (r0 / (6. * M) + 1.)) * (1 / (r0 / (3. * M) - 1.))))

    # Initialize variables for ray tracing
    delta_deg_max = 179  # Maximum angular offset
    lines = []  # Store ray path data
    hitDiskTimes = []  # Times when rays hit the accretion disk
    hitDiskY = []  # Y-coordinates where rays hit the accretion disk

    # Generate ray paths while decreasing delta_deg_max
    while delta_deg_max > 0:
        r, phi, t = RayPathUntilAccDisk(r0=r0, phi0=phi0, t0=t0, Delta0=np.deg2rad(delta_deg_max), r_max=r0 + 100)
        x = r * np.cos(phi)
        y = r * np.sin(phi)
        delta_deg_max -= 0.5

        # If ray hits the accretion disk, record time and position
        if x[-1] < 0.001:
            hitDiskTimes.append(t[-1])
            hitDiskY.append(y[-1])

        # Store ray path data for plotting
        lines.append((x, y, t))
        lines.append((x, -y, t))  # Add symmetrical ray path

        # Terminate if certain conditions are met
        if y[-1] >= 50 and x[-1] < 0.001:
            break

    # Normalize time values for consistent color mapping
    t_all = np.concatenate([t for _, _, t in lines])
    norm = plt.Normalize(t_all.min(), t_all.max())

    # Plot all ray paths with color-coded time
    fig, ax = plt.subplots()
    for x, y, t in lines:
        points = np.array([x, y]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)
        lc = LineCollection(segments, cmap='jet', norm=norm)
        lc.set_array(t)
        ax.add_collection(lc)

    # Add additional plot features
    ax.autoscale()
    ax.set_title('Photon Paths Around a Schwarzschild Blackhole')
    cbar = plt.colorbar(lc, ax=ax)
    cbar.set_label('Time')
    ax.add_artist(plt.Circle((0, 0), 2, fill=False))  # Add black hole's event horizon
    plt.annotate(f'r₀ = {r0}', xy=(1, 0), xycoords='axes fraction', fontsize=8, ha='right', va='bottom')
    plt.ylim(-r0 - 5, r0 + 50)
    plt.xlim(-r0 - 5, r0 + 5)
    plt.grid()
    plt.show()

    # Return sorted results for further processing
    return np.sort(np.array(hitDiskY)), np.sort(np.array(hitDiskTimes))

# Function to plot time delays of rays reaching the accretion disk
def plotTimeDelays(Y, Delays, r0, nrays):
    """
    Plots time delays for rays reaching the accretion disk.
    Parameters:
        Y: List of Y-coordinates for ray impacts
        Delays: List of time delays for rays
        r0: List of initial radial distances
        nrays: Number of rays simulated
    """
    for i in range(len(Y)):
        x_smooth = np.linspace(Y[i].min(), Y[i].max(), 300)
        spl = make_interp_spline(Y[i], Delays[i])
        y_smooth = spl(x_smooth)
        plt.scatter(Y[i], Delays[i], label='_nolegend_')
        plt.plot(x_smooth, y_smooth)
    plt.legend([f"r₀ = {r}" for r in r0])
    plt.grid()
    plt.title("Time It Takes Rays to Reach the Accretion Disk")
    plt.xlabel("Y")
    plt.ylabel("Delay")
    plt.show()

# Main computation and plotting for r₀ ranging from 2.2 to 9.2
r0 = 2.2
all_delays = []
while r0 < 9.2:
    Y, Delays = multipleRayPathsAroundBH(r0, 100)
    r0 += 1
    x_smooth = np.linspace(Y.min(), Y.max(), 300)
    spl = make_interp_spline(Y, Delays)
    y_smooth = spl(x_smooth)
    mask = y_smooth <= 50  # Filter time delays for plotting
    x_filtered = x_smooth[mask]
    y_filtered = y_smooth[mask]
    spl = make_interp_spline(x_filtered, y_filtered)
    all_delays.append(spl(np.linspace(x_filtered.min(), x_filtered.max(), 300)))

# Display all time delay results as a heatmap
plt.imshow(np.array(all_delays).T, cmap='rainbow', interpolation='nearest', aspect='auto', origin='lower', extent=[2.2, 9.2, 2, 50])
plt.colorbar().ax.set_ylabel('Time Delay', rotation=270)
plt.xlabel("r₀")
plt.ylabel("Y")
plt.title("Time It Takes to Reach Accretion Disk")
plt.show()
