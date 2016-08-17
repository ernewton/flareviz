# -*- coding: utf-8 -*-
"""
Created on Wed Aug 17 09:48:21 2016

@author: enewton
"""

import pickle
import matplotlib.pylab as plt
import matplotlib as mpl
import numpy as np
from astropy.time import Time
from zachopy import oned as zponed
import math


def short_period(lc, xlim=None, days=False, save=False):

  lc.hjd_off = 0.
  try:
    print xlim-1
  except:
    xlim = [lc.hjd[0], lc.hjd[-1]]      

  if days:
    cc = 1.
    xlab = 'Time (days)'
  else:
    cc = 24
    xlab=  'Time (hours)'
  
  fig = plt.figure(figsize=(6,4.5))
  ax = fig.add_subplot(111)
  mpl.rc('xtick',labelsize = 16)
  mpl.rc('ytick',labelsize = 16)

  time = lc.phase_fold(binned=True) # reference to time
  yy = lc.flux_reduce(binned=True)
  yy = np.power(10, -1*yy/2.5)-1

  plt.scatter(time*lc.period*cc, yy, color='#004080', s=40, alpha=0.5)
  plt.scatter((time+1)*lc.period*cc, yy, color='gray', s=40, alpha=0.2)
  plt.scatter((time-1)*lc.period*cc, yy, color='gray', s=40, alpha=0.2)
  tmp = np.abs(yy) < (5*np.nanstd(np.abs(yy)))
  lim = 5*np.nanstd(np.abs(yy[tmp]))
  print (5*np.nanstd(np.abs(yy))), lim
  plt.ylim(-lim, lim)
  plt.xlim(-1*lc.period*0.4*cc, (lc.period)*1.4*cc)     
  ax.set_ylabel('Rel. brightness (mag)', fontsize=18)
  ax.set_xlabel(xlab, fontsize=18)
  
  y = [0.02, 0.01, 0.0, -0.01, -0.02]
  plt.yticks(y)
  plt.tight_layout()
  if save:
    plt.savefig('short_phased.jpg')

  
def long_period(lc, xlim=None, alpha=0.6, binwidth=None, color='#004080', binned=False, save=False):
  

  lc.hjd_off = 0.

  fig = plt.figure(figsize=(12,5))
  ax = fig.add_subplot(111)
  mpl.rc('xtick',labelsize = 16)
  mpl.rc('ytick',labelsize = 16)

  if binned:
    xx = Time(lc.hjd_bin+2400000.5,format='jd')
    xxx= lc.hjd_bin
    yy = lc.flux_reduce(binned=True)
  else:
    xx = Time(lc.hjd+2400000.5,format='jd')
    xxx= lc.hjd
    yy = lc.flux_reduce()
  yy = np.power(10, -1*yy/2.5)-1
  print yy
  ax.plot_date(xx.plot_date,yy, zorder=0,
               markeredgecolor='none', color=color, alpha=alpha)

  if binwidth is not None:
    xa, y, error = zponed.binto(x=xxx[np.isfinite(yy)], y=yy[np.isfinite(yy)], 
                                         sem=True, binwidth=binwidth)
    x = Time(xa+2400000.5, format='jd')
    plt.errorbar(x.plot_date, y, error, zorder=1, 
                 elinewidth=3, capthick=3, capsize=3, linewidth=0, color='#1f78b4')
  y = [0.03, 0.02, 0.01, 0.0, -0.01, -0.02, -0.03]
  ylab = ['3%', 0.02, 0.01, 0.0, -0.01, -0.02, -0.03]
  plt.yticks(y, ylab)
  
  tmp = np.abs(yy) < (5*np.nanstd(np.abs(yy)))
  lim = 5*np.nanstd(np.abs(yy[tmp]))
  plt.ylim(-lim, lim)
  print np.nanstd(np.abs(yy))    
  plt.gcf().autofmt_xdate()
  
  ax.set_ylabel('Relative brightness (%)', fontsize=18)
  ax.set_xlabel('Date', fontsize=18)
  
  plt.tight_layout()
  if save:
    plt.savefig('long_stretch.jpg')


def star(lspm, load=True, south=False, model=False, plotlong=None):
  
  # southern stars have a different naming convention
  if south:
    prefix = "2massj"
    suffix = "daily"
    ls = lspm
  else:
    prefix = "lspm"
    suffix = "lc"
    ls = str(int(lspm))
  
  # load the data from pickle?
  if load:
    lc = pickle.load( open("data/"+prefix+ls+".pkl", "rb") )
  # else must read it in!
  else:
    f = "data/"+prefix+ls+"_"+suffix+".fits"
    lc = melc.LightCurve(f, id=lspm, south=south)
    buf = lc.prep_period()
    fit = melc.fit_period(buf, pmax=200)
    lc.update_model(fit)
    pickle.dump( lc, open( "data/"+prefix+ls+".pkl", "wb" ) )
  plt.close()
  plt.close()

  if (lc.period > 10) or plotlong==True:
    long_period(lc)
    if model:
        step = (lc.period/50)
        axrange = np.arange(plt.xlim()[0],plt.xlim()[1], step)
        y = lc.amp*np.sin(2*math.pi*(axrange)/lc.period + lc.phase)
        plt.plot(axrange,y,c='k', zorder=1)
  else:
    short_period(lc)

  print "Period is: ", lc.period, "days"
  print "   ... or  ", lc.period*24, "hours"


def add_model(lspm, south=False):
  # southern stars have a different naming convention
  if south:
    prefix = "2massj"
    suffix = "daily"
    ls = lspm
  else:
    prefix = "lspm"
    suffix = "lc"
    ls = str(int(lspm))
  
  lc = pickle.load( open("data/"+prefix+ls+".pkl", "rb") )
  step = (lc.period /50)
  axrange = np.arange(plt.xlim()[0],plt.xlim()[1], step)
  y = lc.amp*np.sin(2*math.pi*(axrange)/lc.period + lc.phase)
  plt.plot(axrange,y,c='k', zorder=1)

plt.close()