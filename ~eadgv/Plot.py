# -*- coding: utf-8 -*-
"""
Created on Mon Dec 19 11:25:53 2022

@author: Mohammed A Jabed Jabed.abu@gmail.com 
"""


#sns.set(context='notebook', style='white',
#        color_codes=True, rc={'legend.frameon': False})
import matplotlib.pyplot as plt
import math 
import numpy as np 

plt.rcParams["font.family"] = 'Arial'
plt.rcParams.update({'font.size': 20})
Color = plt.rcParams['axes.prop_cycle'].by_key()['color']
#Color = ['red', 'skyblue', 'green','gray','orange']

def plot_line(x,y, xlim=None, ylim=None, ls=None, lw=None, colors = None, xticks=None, yticks=None, 
              xticklabels=None, yticklabels=None,xlabel=None, ylabel=None, figsize=(7,5), name=None):  

    argnan = lambda a,b : a if a != None else b 
    
    fig=plt.figure(figsize=figsize) 
    ax= fig.add_subplot(1,1,1) 

    ax.plot(x,y,color=argnan(colors,'b'),lw=argnan(lw,2),ls=argnan(ls,'-')) 
    if xlim: 
        ax.set_xlim(xlim[0],xlim[1]) 
    if ylim: 
        try: 
            ax.set_ylim(ylim[0],ylim[1])
        except: 
            ax.set_ylim(bottom=0) 
    if yticks: 
        ax.set_yticks(yticks) 
    if xticks: 
        ax.set_xticks(xticks) 
    if yticklabels: 
        ax.set_yticklabels(yticklabels) 
    if xticklabels: 
        ax.set_xticklabels(xticklabels) 
    if xlabel: 
        ax.set_xlabel(xlabel) 
    if ylabel: 
        ax.set_ylabel(ylabel) 
    if name: 
        plt.tight_layout()
        plt.savefig(name,dpi=330)         
    else: 
        return ax 


def plot_waterfall(x,z,xlim=[-3,3],ylim=None, gamma = 45, theta=45, alpha=0.4,
                   figsize=(7,5),name=None, ls=None, lw=None, color = None, xlabel=None, ylabel=None ): 
    print('This is experimental work, outcome may not produce expected quality')
    n_prec =  0 if math.floor(math.log10(z.max() - z.min())) > 1 else 2-  math.floor(math.log10(z.max() - z.min()))
    z_shift = z.min() 
    if z_shift<0: 
        z = z-z_shift  
    gamma= math.radians(gamma) 
    theta = math.radians(theta) 
    
    Rx = np.array([[1, 0,0],
                   [0, math.cos(theta), math.sin(theta)],
                   [0,-math.sin(theta),math.cos(theta)]])   
    Rz = np.array([[math.cos(gamma), math.sin(gamma), 0],
                   [-math.sin(gamma),math.cos(gamma),0],
                   [0,0,1]])  
    # Rotation matrix around x and z axis 
    Rxz = np.array([[math.cos(gamma),math.cos(theta)*math.sin(gamma), math.sin(theta)*math.sin(gamma)],
                    [-math.sin(gamma), math.cos(theta)*math.cos(gamma),math.sin(theta)*math.cos(gamma)],
                    [0,-math.sin(theta),math.cos(theta)]])


    fig = plt.figure(figsize=(7,5)) 
    ax = fig.add_subplot(1,1,1)  
    
    x_range = xlim[1] - xlim[0] 
    if ylim:
        yz_ratio = (ylim[1]-ylim[0])*0.2
    else: 
        yz_ratio = (z.max()-z.min())*0.2
    
    if z.shape[1]*x_range*0.02 > x_range*0.25: 
        yx_ratio = 0.25/z.shape[1] 
    else: 
        yx_ratio = 0.02      
    print(yx_ratio)
    for i in np.arange(0,z.shape[1],1): 
        ay = [0,i*x_range*yx_ratio,0]
        ay_trans = np.matmul(Rxz,ay) 
        X = x + np.abs(ay_trans[0]) 
        Y = z[:,i] + np.abs(ay_trans[2]) * yz_ratio
        ax.fill_between(X, Y, 0, color='white',edgecolor='w', zorder=-i, alpha=alpha)
        if color: 
            ax.plot(X,Y,zorder=-i,color=color ) 
        else: 
            ax.plot(X,Y,zorder=-i)
    #   print(ay_trans)
        
    ax.set_ylim(bottom = z.min())
    ylim = ax.get_ylim()

    X2 = xlim[0] + np.abs(ay_trans[0] ) 
    Y2 = np.abs(ay_trans[2]) * yz_ratio 

    ax.set_xlim(xlim[0],X2+x_range)
    axlw = ax.spines['bottom'].get_linewidth() 

    ax.plot([xlim[0],X2], [z.min(),z.min()+Y2],color='k', zorder=51)  
    ax.plot([xlim[1],X2+x_range], [z.min(),z.min()+Y2],color='k', zorder=51)  

    ax.plot([xlim[0],X2], [ylim[1],ylim[1]+Y2],  lw = axlw, color='k', zorder=51)  
    ax.plot([xlim[1],X2+x_range], [ylim[1],ylim[1]+Y2], lw = axlw, color='k', zorder=51)   

    ax.vlines(x=X2,ymin=Y2, ymax =ylim[1] + Y2,  lw = axlw, color='k', zorder = -2)  
    ax.vlines(x=xlim[1] , ymin=z.min(), ymax =ylim[1],  lw = axlw, color='k', zorder = 10)  
    #print( X2)
    ax.hlines(y=Y2, xmin = X2 ,xmax = x_range + X2,  lw = axlw, color = 'k', zorder = -5)  
    ax.hlines(y=ylim[1], xmin= xlim[0], xmax = xlim[1],  lw = axlw, color='k', zorder = -10)  

    ax.set_ylim(0,ylim[1]+Y2)
 
    ax.spines['bottom'].set_bounds(xlim[0], xlim[1])
    ax.spines['top'].set_bounds(X2, X2+x_range) #xlim[1]+np.abs(ay_trans[0]))
    ax.spines['right'].set_bounds(Y2, ylim[1] + Y2 )  
    ax.spines['left'].set_bounds(z.min(), ylim[1] )#+ Y2 )   

    new_ticks = ax.get_yticks()+z_shift 
    print(new_ticks)
    ax.set_yticklabels([f"{j:.{n_prec}f}" for j in new_ticks])
    
    if xlabel: 
        ax.set_xlabel(xlabel) 
    if ylabel: 
        ax.set_ylabel(ylabel) 
    if name: 
        #plt.tight_layout()
        plt.savefig(name,dpi=630)         
    else: 
        return ax 
        