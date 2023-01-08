
import math, os.path, sys, time
from datetime import datetime, timedelta
from glob import glob
from argparse import ArgumentParser
import numpy as np 
import os 

try : 
    from read_outcar import * 
    from Plot import * 
    from fileconversion import * 
except: 
    from .read_outcar import * 
    from .Plot import * 
    from .fileconversion import * 
    
class Warning: 
    def file_exist( name): 
        if  not  os.path.isfile(name): 
            print('warning! file not found - %s' %name) 
            exit()
        else : pass 
            


def main(): 
    parser = ArgumentParser() 
    parser.add_argument("-job", dest="job", action="store", default=None, type=str.lower,
                        help="task to perform")
    parser.add_argument("--file", dest="fname", default=None, type=str, 
                        help = 'file name to convert'),
  
    parser.add_argument("--fout", dest="fout", default=None, type=str, 
                        help = 'output file name'),
    parser.add_argument("--ftype", dest="ftype", default="xyz", type=str, 
                        help = 'output file type')

    parser.add_argument("--select", dest="select", default=False, type=lambda x: True if x=='True' else False,  
                        help = 'output file type')

    parser.add_argument("--direct", dest="direct", default=True, type=lambda x: False if x=='False' else True, 
                        help = 'output file type')

    
    (options, args) = parser.parse_known_args() 
    print(options)
    print( args)
    
    if options.job == 'vasp2gaus': 
        if options.fname is None: 
            options.fname = input('What is the POSCAR file name? :').strip() 
        fname = options.fname
        
        if options.fout is None: 
            options.fout = os.path.splitext(options.fname)[0]
        name_out = os.path.splitext(options.fout)[0]  

        if options.ftype =='gaus': 
            Gaus_out = True 
            if os.path.isfile('%s.com'%name_out):
                print ('%s.com is exist'%name_out) 
                name_out_ = input('What should be the output file name?:').replace(" ", "").rstrip('.com')  
                if not len(name_out_)==0: 
                    name_out= name_out_     
                if os.path.isfile('%s.com'%name_out):  # overwritting the file if no name is given 
#                    name_out = name_out+'_'+''.join(random.choices(string.ascii_lowercase, k=5)) 
                    os.remove('%s.com'%name_out) 
        if options.ftype =='xyz': 
            Gaus_out = False
            if os.path.isfile('%s.xyz'%name_out):
                print ('%s.xyz is exist'%name_out) 
                name_out_ = input('What should be the output file name?:').replace(" ", "").rstrip('.xyz')  
                if not len(name_out_)==0: 
                    name_out= name_out_     
                if os.path.isfile('%s.xyz'%name_out): 
                    os.remove('%s.xyz'%name_out) 
                    print('Overwriting the file %s.xyz'%name_out) 
   
        contcar2com(fname, name_out,Gaus_out=Gaus_out) 
    
    elif options.job == 'gaus2poscar': 
        fname = options.fname 
        
        Warning.file_exist(fname)  
        
        file = open(fname,'r') 
        lines = file.read() 
        if not 'tv' in lines.lower() :
            print('Warning! Periodic cell vectors are not found in file %s' %fname)     
        file.close() 
        
        name_out = os.path.splitext(options.fout)[0] 
        
        
        com2poscar(fname, name_out=name_out, select=options.select, direct = options.direct)  


    

if __name__ == "__main__":
    main() 
    
    
    