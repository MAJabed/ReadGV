# -*- coding: utf-8 -*-
"""
Created on Wed Mar 22 16:04:58 2023

@author: Mohammed A Jabed 
"""
import numpy as np 


class read_poscar(): 
    def __init__ (self, poscar='POSCAR'): 
        self.poscar = poscar  
        self.coordinates = self.poscar_data() 
        self.comments, self.scaling_factor, self.pbc, self.ions_symbol, self.ions_n, self.selective, self.dir_cart, self.ions_coord = self.poscar_data()  
        
    def poscar_data(self): 
        poscar = open(self.poscar, 'r') 
        
        comments = poscar.readline().strip('\n')  
        scaling_factor  = float(poscar.readline()) 
        pbc = np.loadtxt([poscar.readline() for i in range(3)],dtype=float) 
        ions_symbol = poscar.readline().split()  
        ions_n = np.array(poscar.readline().split(),dtype=int)   
        line = poscar.readline().strip() 
        if line[0].lower=='s': 
            selective = True
            dir_cart = poscar.readline().lstrip()[0].lower()  
        else: 
            selective = False 
            dir_cart = line.lstrip()[0].lower()
        ions_coord = [poscar.readline() for i in range(sum(ions_n))] 
        ions_coord = np.loadtxt(ions_coord)[:,:3] 
        return comments, scaling_factor, pbc, ions_symbol, ions_n, selective, dir_cart, ions_coord 
