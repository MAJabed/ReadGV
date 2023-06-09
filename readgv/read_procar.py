# -*- coding: utf-8 -*-
"""
Created on Wed Mar 12 16:04:58 2023

@author: Mohammed A Jabed 
"""


import numpy as np 
import re 
import warnings

class read_procar: 
    def __init__ (self, filename='PROCAR'): 
        self.filename = filename 
        self.procar = open(self.filename).readlines() 
        
    def nkpoints(self): 
        nkpoints = int(re.findall(r'of k-points:(.*) # of bands:', self.procar[1])[0]) 
        return  nkpoints 
    
    def nbands(self) :
        nbands = int(re.findall(r'# of bands:(.*)# of ions:',self.procar[1])[0])
        return nbands 
    
    def nions(self): 
        nions = int(re.findall(r'of ions:(.*)\n', self.procar[1])[0]) 
        return nions 
    
    def energy(self, n=1): 
        ene_occ  = self.procar[5+(n-1)*(self.nions+5)]
        energy  = float(re.findall(r'energy (.*) # occ.',ene_occ)[0])
        return energy 
    
    def occupation(self, n=1): 
        ene_occ  = self.procar[5+(n-1)*(self.nions+5)]
        occ  = float(re.findall(r'# occ. (.*)\n',ene_occ)[0])
        return occ  
        
    def read_band(self, n): 
        nkpoints = self.nkpoints() 
        nbands = self.nbands() 
        nions = self.nions()  
        ene_occ  = self.procar[5+(n-1)*(nions+5)]
        
        energy  = float(re.findall(r'energy (.*) # occ.',ene_occ)[0])
        occ  = float(re.findall(r'# occ. (.*)\n',ene_occ)[0]) 
        
        proj_bands =  self.procar[5+(n-1)*(nions+5)+3 : 5+(n-1)*(nions+5)+3+nions] 
        proj_bands = np.loadtxt(proj_bands, dtype=float) 
        
        proj_tot = self.procar[ 5+(n-1)*(nions+5)+3+nions]
          
        proj_tot = np.asarray(proj_tot.split()[1:]).astype(float) 
        return proj_bands, proj_tot, energy, occ
