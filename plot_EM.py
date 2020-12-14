#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 11 12:53:48 2018

@author: yannis
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.ticker as ticker

def addPotentialColorMap(myAxis):
    """
    Create a colormap for the potential on myAxis
    Add labels on X,Y axes
    """
    plt.xlabel('X(m)', fontsize=14)
    plt.ylabel('Y(m)', fontsize=14)
    return myAxis.pcolor(xGrid, yGrid, potential,
                       norm=mcolors.SymLogNorm( linthresh=1E8, linscale=1,
                                              vmin=potential.min(), vmax=potential.max() ),
                       cmap='jet')

def fixedPotential_cBar(myFig,fraction="5%"):
    """
    This function fixes the colorbar of the potential colormap
    """
    
    ## I create another axis object for the colorbar 
    ## in order to use later tbe plt.tight_layout() function
    divider = make_axes_locatable(plt.gca())
    cax = divider.append_axes("right", fraction, pad="3%")
    cbar = myFig.colorbar(pColorMap, cax=cax)
    cbar.ax.set_ylabel('Potential (Volts)', fontsize=14)
    
    ## NOTE: the fact that symmetric log-scale colorbar tick labels 
    ## near zero don't appear correctly is a known pyplot issue
    ## see: https://matplotlib.org/users/colormapnorms.html
    ## here the large threshold partially alleviates the problem 
    ## but still we will do some modifications (remove minor ticks, fix v==0 tick) 
    
    ## fixing the ticks of the colorbar
    loc = ticker.SymmetricalLogLocator(linthresh=1E8,base=10)
    cbar.ax.yaxis.set_major_locator(loc)
    cbar.set_ticks(cbar.ax.yaxis.get_major_locator().tick_values(potential.min(), potential.max()))

    return

e0 =  8.85E-12


x = np.linspace(0, 1, 101)
y = np.linspace(0, 1, 101)
xGrid, yGrid = np.meshgrid(x, y)

# source coordinates:
# making sure that they don't coincide with grid 
# points to avoid divisions by zero

q1 = -1
posX1 = 0.4555
posY1 = 0.500

q2 = 1
posX2 = 0.5555
posY2 = 0.500

potential1 = q1/(4 * np.pi * e0 * np.sqrt( (xGrid-posX1)**2 + (yGrid-posY1)**2 ))
potential2 = q2/(4 * np.pi * e0 * np.sqrt( (xGrid-posX2)**2 + (yGrid-posY2)**2 ))

potential = potential1 + potential2

deltaX = x[1]-x[0]
deltaY = y[1]-y[0]

electricX = -(potential[1:-1,2:] - potential[1:-1,:-2])/ (2*deltaX)
electricY = -(potential[2:,1:-1] - potential[:-2,1:-1])/ (2*deltaY)

electricNorm = np.sqrt(electricX**2 + electricY**2)
electricAngles = np.arctan2(electricY,electricX)

#### Plots

## Plot the electric potential
fig1 ,ax1= plt.subplots(num='Potential', figsize=(6, 6))
ax1.set_title('Electric Potential V', fontsize=14)
ax1.set_aspect(1)

## setting log scale for potential values 
pColorMap = addPotentialColorMap(ax1)

fixedPotential_cBar(fig1)

plt.tight_layout()
plt.show()


fig2 , ax2 = plt.subplots(1,1, figsize=(12, 6), num='Electric Field Magnitude and Angle')
ax2.set_title('Magnitude of the Electric Field', fontsize=14)

## Plot Magnitude Field
ax2.set_aspect(1)
eColorMap = ax2.pcolor(xGrid[1:-1,1:-1], yGrid[1:-1,1:-1], electricNorm,
                   norm=mcolors.LogNorm(vmin=electricNorm.min(), vmax=electricNorm.max()),
                   cmap='gnuplot2')
plt.xlabel('X(m)', fontsize=14)
plt.ylabel('Y(m)', fontsize=14)


divider = make_axes_locatable(plt.gca())
cax = divider.append_axes("right", "5%", pad="3%")
cbar = fig2.colorbar(eColorMap, cax=cax)
cbar.ax.set_ylabel('E (Volts/m)', fontsize=14)


## Plot Angle Field
ax2 = divider.append_axes("right", "100%", pad="40%")
ax2.set_title('Direction Map of the Electric Field', fontsize=14) 
ax2.set_aspect(1)
aColorMap = ax2.pcolor(xGrid[1:-1,1:-1], yGrid[1:-1,1:-1], electricAngles, cmap='hsv')
plt.xlabel('X(m)', fontsize=14)
plt.ylabel('Y(m)', fontsize=14)


cax = divider.append_axes("right", "5%", pad="3%")
cbar = fig2.colorbar(aColorMap, cax=cax)
cbar.ax.set_ylabel('Angles (Radians)', fontsize=14)

plt.tight_layout()
plt.show()


## Mask the fields very close to the charges to get rid of very large gradients
maskRadius = 0.04
mask = np.zeros(xGrid[1:-1,1:-1].shape, dtype=bool)
maskPoints1 = (xGrid[1:-1,1:-1]-posX1)**2 + (yGrid[1:-1,1:-1]-posY1)**2 <= maskRadius**2
maskPoints2 = (xGrid[1:-1,1:-1]-posX2)**2 + (yGrid[1:-1,1:-1]-posY2)**2 <= maskRadius**2
mask[maskPoints1] , mask[maskPoints2]= True , True
electricXMasked = np.ma.array(electricX, mask=mask)
electricYMasked = np.ma.array(electricY, mask=mask)


##Vector Plot
fig3 , ax3 = plt.subplots(1,1, figsize=(12, 6),num='Electric Vector Field')
ax3.set_title('Electric Field Vector Plot', fontsize=14) 

pColorMap = addPotentialColorMap(ax3)

plt.quiver(xGrid[1:-1,1:-1], yGrid[1:-1,1:-1], electricXMasked, electricYMasked,scale=1E13, scale_units='inches',units='width', pivot='mid', alpha=.5)
plt.xlim(min(posX1,posX2)-0.1,max(posX1,posX2)+0.105)
plt.ylim(min(posY1,posY2)-0.1,max(posY1,posY2)+0.105)

fixedPotential_cBar(fig3,fraction='2%')

plt.tight_layout()
plt.show()


##Stream Plot
fig4 , ax4 = plt.subplots(1,1, figsize=(6, 6),num='Electric Field StreamPlot')
ax4.set_title('Electric Field Stream Plot', fontsize=14) 
ax4.set_aspect(1)
pColorMap = addPotentialColorMap(ax4)

plt.streamplot(xGrid[1:-1,1:-1], yGrid[1:-1,1:-1], electricXMasked, electricYMasked,color='k')

fixedPotential_cBar(fig4)

plt.tight_layout()
plt.show()

