# -*- coding: utf-8 -*-
"""
Created on Fri Feb  3 15:18:39 2023

@author: Mohammed
"""

import numpy as np 
import math 

def Gauss_distr(X, Y = None, sigma=0.1,): 
    ##  g(x)=1/(sigma*sqrt(2*pi))*{exp(-(x-b)^2/(2*sigma^2))}
    if Y==None: 
        Y = np.ones(X.shape) 
    norm = (sigma*np.sqrt(2*np.pi))**-1 
    
    if np.ndim(X) ==1: 
        x=np.arange(X[0], X[-1],0.05)
        y = np.zeros((len(x)))
        #Dressing with gaussian function
        for J in X:
            yy =norm * np.exp(-1*(x-J)**2/(2*sigma**2))
            y +=yy 
        return x,y
    
    elif np.ndim(X) ==2: 
        x=np.arange(X.min(), X.max(),0.05)
        yy = np.zeros((len(x), X.shape[1])) 
        for i in range(X.shape[1]): 
            y = np.zeros((len(x)))
        #Dressing with gaussian function
            for J in X[:,i]:
                yy =norm * np.exp(-1*(x-J)**2/((2*np.pi*sigma)**0.5))
                y +=yy 
            yy[:,i] = y  
        return x, yy 

def Dress_abs(X,OS, unit=None, linewidth=0.4,): 
    
    """ 
    Absorption calculated using following equation. 
    Ei = {sqrt(pi)*e^2*N / (1000*ln(10).c^2*me)}* {fi/sigma} * exp[-{(x-xi) / sigma}^2 ] 
    electron charge: e = 4.803204×10-10 esu (cm3/2.g1/2.s-1) 
    electron mass: me = 9.10938×10-28 g
    Avogrado number: N = 6.022*10**23 mol-1 
    speed of light: c=29979245800.0 cm·s-1 
    sigma in eV unit, and converted to nm unit 
    """ 
    e = 4.803204*10**-10
    me = 9.10938*10**-28 
    c = 29979245800.0  # cm·s-1 
    N = 6.022*10**23  
    
    sigma_cm =  10**-7*1240/linewidth #linewidth in ev unit 
         
    norm = np.sqrt(np.pi)*N*e**2 / (1000*np.log(10)*me*(1/sigma_cm)*c**2)
    
    ev2nm = lambda a : 1240/a 
    
    if unit: 
        if unit.lower() == 'ev': 
            x =np.arange(min(X[0],0.05), max(X[-1],12),0.005) 
        elif unit.lower() == 'nm': 
            x = np.logspace(np.log2(max(3000,X[0])), np.log2(min(100,X[-1])), num=2000 , base = 2, endpoint=False)
    else: 
        if X.max() > 50:
            unit='nm'
            print('Absorption spectra is dressing as nm unit')
            x = np.logspace(np.log2(max(3000,X[0])), np.log2(min(100,X[-1])), num=2000 , base = 2, endpoint=False)
        else: 
            unit = 'ev'
            x =np.arange(min(X[0],0.1), max(X[-1],12),0.005)

    Y = np.zeros((len(x)))
    #Dressing with gaussian function
    if unit == 'ev': 
        for j in range(len(X)):
            y =norm * OS[j] *  np.exp(-1*((x-X[j])/linewidth)**2) 
            Y +=y 
        return x,Y 
    elif unit=='nm': 
        for j in range(len(X)):
            y =norm * OS[j] *  np.exp(-1*((1/x-1/X[j])/(1/ev2nm(linewidth)))**2)   
            Y +=y 
        return x,Y
    
def unitcell2cartbasis(a,b,c,alpha=90,beta=90,gamma=90,volume=False): 
    ''' Conversion of latice vector of a unit cell and angles to parallelepped vectors in a cartesian real space 
        alpha is angle between b & c 
        beta is angle between a & c and 
        gamma is angle between a & b '''
    cos = lambda r : math.cos(math.radians(r))
    sin = lambda r : math.sin(math.radians(r)) 
    
    cc = (cos(alpha)-cos(beta)*cos(gamma) )/ sin(gamma) 
    va = np.array([a,0,0]) 
    vb = np.array([b*cos(gamma), b*sin(gamma),0]) 
    vc = np.array([c*cos(beta), c*cc, c*(1-(cos(beta))**2-cc*2)**0.5]) 
    if volume==True: 
        v = a*b*c*(1+2*cos(alpha)*cos(beta)*cos(gamma)-(cos(alpha))**2-(cos(beta))**2-(cos(gamma))**2)**0.5    
        return np.array([va,vb,vc]), v 
    else: 
        return np.array([va,vb,vc]) 
    
def Rzyx(alpha,beta,gamma): 
    ''' Rotation matrix, in X Y and Z orders, angles are in degrees ''' 
    cos = lambda r : math.cos(math.radians(r))
    sin = lambda r : math.sin(math.radians(r))
   
    Rzyx = np.array([[cos(alpha)*cos(beta) , 
                  cos(alpha)*sin(beta)*sin(gamma)-sin(alpha)*cos(gamma),
                  cos(alpha)*sin(beta)*cos(gamma)+ sin(alpha)*sin(gamma)],
                [sin(alpha)*cos(beta) , 
                 sin(alpha)*sin(beta)*sin(gamma)+cos(alpha)*cos(gamma),
                 sin(alpha)*sin(beta)*cos(gamma) - cos(alpha)*sin(gamma)],
                [-sin(beta),cos(beta)*sin(gamma),cos(beta)*cos(gamma)]]) 
    return Rzyx 
    

