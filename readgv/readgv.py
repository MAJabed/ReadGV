
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
    parser.add_argument("--bands", dest="bands", default=100, type=int, 
                        help = 'No of states to read, up-below homo')
    parser.add_argument("--frames", dest="frames", default='final', type=str, 
                        help = 'No of MD frame to read, e.g. final, initial, all, -n, n')
    parser.add_argument("--nthkpoint", dest="nthkpoint", default=1, type=int, 
                        help = 'Read nth KPOINT data')
    parser.add_argument("--coord", dest="coord", default='cart', type=str, 
                        help = 'Output coordinate types')

    
    (options, args) = parser.parse_known_args() 
    #print(options)
   #print( args)
    
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

    elif options.job == 'read_bands':
        fname = options.fname 
        nbands = options.bands 
        frames = options.frames
        nthkpoints = options.nthkpoint -1 
        
        #Translate frame argument in int variables 
        get_frame =  lambda a, k, nf : a[k,-1,:,:] if nf == 'final' \
        else ( a[k,0,:,:] if nf == 'initial' \
              else (a[k,:,:,:] if  nf ==  'all' 
                    else (a[k,:int(nf),:,:] if nf.isdigit() 
                          else (a[k,int(nf):,:,:] if nf.lstrip('-').isdigit()
                              else a ))))   
        
        outcar = read_outcar(filename=fname)
        Bands = outcar.band_energies()   # Output dimension is K points x frames x bands x 3 
        
        Bands = get_frame(Bands,nthkpoints, frames) 
        Fermi_1 = outcar.fermi_energies()[0,0] 
        nhomo  = int(Bands[0,:,0][Bands[0,:,1] < Fermi_1][-1])  
        nbands = nbands if nbands < nhomo else nhomo-1 
        
        Bands = Bands[:,nhomo-nbands-1:nhomo+nbands-1,:] 
              
        # Data write in a file options.fout, else print
        if options.fout != None: 
            if Bands.ndim <3: #it indicates one frame only 
                with open(options.fout,'w') as f: 
                    np.savetxt(f,Bands, fmt='%-4.i   %-10.6f    %-10.6f')
                f.close() 
            elif Bands.ndim ==3: 
                with open(options.fout,'w') as f: 
                    for j in range(Bands.shape[0]): 
                        f.write('Frame %i \n'%(j+1)) 
                        np.savetxt(f,Bands[j,:,:], fmt='%-4.i   %-10.6f    %-10.6f') 
                        f.write('\n')
            else: 
                print('Can not find a way to write in a file') 
                print(Bands) 
        else: 
            print(Bands) 

    elif options.job == 'get_geom':
        fname = options.fname 
        fout = options.fout   
        coord = options.coord 
        frames = options.frames
        
        get_frame =  lambda a, n : a[-1,:,:] if n == 'final' \
        else ( a[0,:,:] if n == 'initial' \
              else (a[:,:,:] if  n ==  'all' 
                    else (a[:int(n),:,:] if n.isdigit() 
                          else (a[int(n):,:,:] if n.lstrip('-').isdigit()
                              else a )))) 
        
        outcar = read_outcar(filename=fname)
        Trajectory = outcar.md_traj() 
        Trajectory = get_frame(Trajectory,options.frames) 
        tot_ions =  outcar.tot_ions  
        ion_list = outcar.ion_list 
        if fout != None: 
            traj_out = open(fout,'w') 
            for i in range(Trajectory.shape[0]): 
                traj_out.write('{} \nFrame {} \n'.format(tot_ions,i+1))
                for j in range(tot_ions): 
                    STR = ' {:<6}{:>12}{:>12}{:>12}'.format(ion_list[j] , Trajectory[i,j,0], Trajectory[i,j,1], Trajectory[i,j,2]) 
                    traj_out.write(STR+'\n')
                traj_out.write('\n') 
            traj_out.close()
        else: 
            print(Trajectory)
            

                    
        #    print(a.outcar2incar())
        #    a.md_traj(out='Temmp_md.dat')  
        # F
    

if __name__ == "__main__":
    main() 
    
    
    