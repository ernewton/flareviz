# -*- coding: utf-8 -*-
"""
Created on Mon Aug 15 09:57:29 2016

@author: ellie
"""
import seaborn as sns
from scipy.interpolate import interp1d
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rc('xtick', labelsize=16) 
mpl.rc('ytick', labelsize=16) 

myfile = 'data/ERN_3120_flare_616'

time_d, flux, eflux = np.genfromtxt(myfile, unpack=True)
toff=4
time = (time_d - time_d[0])*24+toff
flux = np.power(10, -1*flux/2.5) + np.median(flux[time<-1])
eflux = flux*eflux*np.log(10.)/2.5

order = np.argsort(time)
time = time[order]
flux = flux[order]
eflux = eflux[order]

thalf = 0.7
xnew = np.arange(-10.8561,10,0.001)
ynew = np.zeros_like(xnew)*0
rise = (1. + 1.941*(xnew/thalf)\
        -0.175*(xnew/thalf)**2.\
        -2.246*(xnew/thalf)**3.\
        -1.125*(xnew/thalf)**4.)*(np.max(flux)-1)
ynew[(rise>0)*(xnew<0)] = rise[(rise>0)*(xnew<0)]
ynew[xnew>=0] = 1.1*(np.max(flux)-1)*np.exp(-0.965*xnew[xnew>=0]/thalf)


plt.figure(figsize=(8,6))
colors = sns.color_palette()

err = 0.002
i = 1
plt.fill_between(xnew+toff, ynew+1-err*2, ynew+1+err*2, alpha=0.2, color=colors[i])
plt.fill_between(xnew+toff, ynew+1-err, ynew+1+err, alpha=0.2, color=colors[i])
plt.plot(xnew+toff, ynew+1, color=colors[i], zorder=1)

plt.errorbar(time, flux, eflux, elinewidth=6,
             linestyle='None', color=colors[3], zorder=2)

plt.xlabel('Time in hours', size=16)
plt.ylabel('Relative brightness', size=16,  labelpad=-10)
plt.xlim(-1,10.2)
plt.ylim(0.978,1.07)
plt.yticks([1,1.02,1.04,1.06], ['Normal','+2%','+4%','+6%'])


plt.annotate(s='', xy=(-4+toff,1), xytext=(-4+toff,1.04), 
             arrowprops=dict(arrowstyle="simple, head_width=1, tail_width=0.4", 
            color=colors[4]))
plt.annotate(s='', xy=(-4+toff,1.065), xytext=(-4+toff,1.02), 
             arrowprops=dict(arrowstyle="simple, head_width=1, tail_width=0.4", 
            color=colors[4]))
plt.text(-3.7+toff, 1.033, '6-7% increase in brightness', rotation=90, size=14,
         verticalalignment='center', color=colors[4])
 
yp = 0.99
plt.annotate(s='', xy=(-0.5+toff,yp), xytext=(1+toff,yp), 
             arrowprops=dict(arrowstyle="simple, head_width=1, tail_width=0.4", 
            color=colors[4]))
plt.annotate(s='', xy=(2+toff,yp), xytext=(toff,yp), 
             arrowprops=dict(arrowstyle="simple, head_width=1, tail_width=0.4", 
            color=colors[4]))
plt.text(0.71+toff, 0.983, '~2 hour duration', size=14,
         horizontalalignment='center', color=colors[4])