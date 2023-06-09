# -*- coding: utf-8 -*-
"""
Created on Mon Mar 27 15:29:31 2023

@author: Mohammed A Jabed jabed.abu@gmail.com
"""

import re 
import numpy as np 

# default atomic weights, elements 0 to 109
defatw = [ 0.00000000,  1.00782504,  4.00260325,  7.01600450,  9.01218250,
          11.00930530, 12.00000000, 14.00307401, 15.99491464, 18.99840325,
          19.99243910, 22.98976970, 23.98504500, 26.98154130, 27.97692840,
          30.97376340, 31.97207180, 34.96885273, 39.96238310, 38.96370790,
          39.96259070, 44.95591360, 47.94794670, 50.94396250, 51.94050970,
          54.93804630, 55.93493930, 58.93319780, 57.93534710, 62.92959920,
          63.92914540, 68.92558090, 73.92117880, 74.92159550, 79.91652050,
          78.91833610, 83.91150640, 84.91170000, 87.90560000, 88.90540000,
          89.90430000, 92.90600000, 97.90550000, 98.90630000,101.90370000,
         102.90480000,105.90320000,106.90509000,113.90360000,114.90410000,
         117.90180000,120.90380000,129.90670000,126.90040000,131.90420000,
         132.90542900,137.90500000,138.90610000,139.90530000,140.90740000,
         141.90750000,144.91270000,151.91950000,152.92090000,157.92410000,
         158.92500000,163.92880000,164.93030000,165.93040000,168.93440000,
         173.93900000,174.94090000,179.94680000,180.94800000,183.95100000,
         186.95600000,189.95860000,192.96330000,194.96480000,196.96660000,
         201.97060000,204.97450000,207.97660000,208.98040000,208.98250000,
         210.98750000,222.01750000,223.01980000,226.02540000,227.02780000,
         232.03820000,231.03590000,238.05080000,237.04800000,242.05870000,
         243.06140000,246.06740000,247.07020000,249.07480000,252.08290000,
         252.08270000,255.09060000,259.10100000,262.10970000,261.10870000,
         262.11410000,266.12190000,264.12470000,  0.00000000,268.13880000]


class read_gaus_log: 
    def __init__(self, file='test.log'): 
        self.file = file 
        f  = open(self.file)  
        ff = f.read() 
        basis_line  = re.findall(r'NBsUse=.*',ff)[0]
        n_basis = int(re.findall(r'NBsUse=\s*(\d*)',basis_line)[0])
        self.no_basis = n_basis 
        self.n_scf = len(re.findall(r'SCF Done:',ff))
        alpha_beta  = re.findall('.*alpha ele.*', ff)[0]
        self.n_alpha_elec = int(re.findall(r'(.*)alpha ', alpha_beta)[0])
        self.n_beta_elec = int(re.findall(r'alpha electrons(.*)beta electrons', alpha_beta)[0])
        self.frag_name = np.array(list(set(re.findall(r'\(Fragment=(.*)\)',ff))))
        self.n_frag = len(self.frag_name)
        f.close() 
        
    def fragment_line(self,line): 
        frag_contri = {} 
        fragmets = re.findall(r'(Fr[\s\S]*?)=',line) 
        for i in range(len(fragmets)-1): 
            fr_str = re.compile('%s=([\s\S]*)%s'%(fragmets[i],fragmets[i+1])) 
            frag_i = int(fragmets[i].replace('Fr',''))
            atom_contr_frag_i  = float(re.findall(fr_str,line)[0]) 
            frag_contri[frag_i] = atom_contr_frag_i 
        # for the last fragment listed 
        fr_str = re.compile('%s=([\s\S]*)'%fragmets[-1]) 
        frag_i = int(fragmets[-1].replace('Fr',''))
        atom_contr_frag_i  = float(re.findall(fr_str,line)[0]) 
        frag_contri[frag_i] = atom_contr_frag_i
        return frag_contri 

    def get_MO_energies(self): 
        file = self.file 
        f = open(file) 
        MO_energies = np.zeros((2,self.no_basis, 2)) 
        
        line  = f.readline() 
        scf_count=0 
        while line:
            if 'SCF Done:' in line: 
                scf_count+=1 
                if scf_count == self.n_scf:
                    line = f.readline() 
                    while line:
                        if 'Population analysis using the SCF Density' in line: 
                            line = f.readline() 
                            qa = 0
                            qb = 0
                            while line: 
                                if ('Alpha' in line)  & ('eigenvalues' in line) : 
                                    line_temp =line.split('--')[1]
                                    ll = [line_temp[1+10*i:10*(i+1)+1] for i in range(int(len(line_temp)/10))]  # line.split('--')[1].split() 
                                    for i in ll: 
                                        MO_energies[0,qa,0] = qa  
                                        MO_energies[0,qa,1] = float(i)  
                                        qa+=1 
                                    line = f.readline() 
                                elif  ('Beta' in line)  & ('eigenvalues' in line) :
                                    line_temp =line.split('--')[1]
                                    ll = [line_temp[1+10*i:10*(i+1)+1] for i in range(int(len(line_temp)/10))]  # line.split('--')[1].split()
                                    for i in ll: 
                                        MO_energies[1,qb,0] = qb 
                                        MO_energies[1,qb,1] = float(i) 
                                        qb+=1 
                                    line = f.readline() 
                                elif qb>= self.no_basis :
                                    f.close() 
                                    return MO_energies
                                else:
                                    line = f.readline() 
                        else:
                            line = f.readline()  
            else: 
                line = f.readline() 

    def get_fragment_contribution_MOs(self): 
        file = self.file 
        f = open(file) 
        Atomic_contribution = np.zeros((2,self.no_basis, 2+self.frag_name.astype(int).max()))  
        line  = f.readline() 
        scf_count=0 
        while line: 
            if 'SCF Done:' in line: 
                scf_count+=1 
                if scf_count == self.n_scf: 
                    line = f.readline() 
                    while line:
                        if 'Atomic contributions to Alpha molecular orbitals:' in line:  
                            line = f.readline() 
                            while line: 
                                if 'Alpha ' in line: 
                                    n_mo = int(re.findall(r'(.*)OE',line )[0][10:]) # int(re.findall(r'occ|vir([\s\S]*)OE',line)[0]) 
                                    en = float(re.findall(r'OE=([\s\S]*)is',line)[0]) 
                                    Atomic_contribution[0,n_mo-1,0] = n_mo
                                    Atomic_contribution[0,n_mo-1,1] = en
                                
                                    line = f.readline()
                                    frag_contri = self.fragment_line(line) 
                                    
                                    for i in frag_contri.keys(): 
                                       # print(n_mo, i, frag_contri[i]) 
                                        Atomic_contribution[0,n_mo-1,i+2-1] = frag_contri[i] 
                                    line = f.readline() 
                                    if 'Fr' in line: 
                                        frag_contri = self.fragment_line(line)
                                        for i in frag_contri.keys(): 
                                            Atomic_contribution[0,n_mo-1,i+2-1] = frag_contri[i] 
                                        line = f.readline() 
                                elif len(line.strip())==0 : 
                                    break 
                                else: 
                                    print('Something is wrong ')
                                    break  

                        elif 'Atomic contributions to Beta  molecular orbitals:' in line:   
                            line = f.readline() 
                            while line: 
                                if 'Beta ' in line:
                                    n_mo = int(re.findall(r'(.*)OE',line )[0][10:]) # int(re.findall(r'occ|vir([\s\S]*)OE',line)[0]) 
                                    en = float(re.findall(r'OE=([\s\S]*)is',line)[0]) 
                                    Atomic_contribution[1,n_mo-1,0] = n_mo
                                    Atomic_contribution[1,n_mo-1,1] = en
                                
                                    line = f.readline()
                                    frag_contri = self.fragment_line(line) 
                                    
                                    for i in frag_contri.keys(): 
                                       # print(n_mo, i, frag_contri[i]) 
                                        Atomic_contribution[1,n_mo-1,i+2-1] = frag_contri[i] 
                                    line = f.readline() 
                                    if 'Fr' in line: 
                                        frag_contri = self.fragment_line(line)
                                        for i in frag_contri.keys(): 
                                            Atomic_contribution[1,n_mo-1,i+2-1] = frag_contri[i] 
                                        line = f.readline() 
                                elif len(line.strip())==0 : 
                                    return Atomic_contribution
                                else: 
                                    print('Beta Something is wrong ')
                                    break  
                        else: line = f.readline()
                    else: line = f.readline() 
                else: line = f.readline() 
            else: line = f.readline() 
        f.close() 
    
    def SCF_energy(self): 
        file = self.file 
        f = open(file) 
        SCF_energy = [] 
        line  = f.readline() 
        while line: 
            if 'SCF Done:' in line: 
                E = float(line.split()[4]) 
                SCF_energy.append(E) 
        return np.array(SCF_energy) 
    



