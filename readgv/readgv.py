
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
    from mathstats import * 
    
except: 
    from .read_outcar import * 
    from .Plot import * 
    from .fileconversion import * 
    from .mathstats import * 
    
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
                        help = 'input file name'),
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

    parser.add_argument("--plot", dest="plot",  type=lambda x: True if x=='True' else x, default=False, 
                        help="Plot output")
    parser.add_argument("--sigma", dest='sigma',default=0.1, type=float,
                        help = 'Linewidth for the Gaussian broadening, default 0.1 eV')
    parser.add_argument("--shiftfermi", dest="shiftfermi", action='store_true', default=False, 
                        help = 'Plot shift to Fermi energies')
    parser.add_argument("--xlim", dest="xlim",  nargs='+', type=float, default=None, 
                        help = 'Lower and upper limit of the X axis')
    parser.add_argument("--alpha", dest="alpha",  type=float, default=0.4, 
                        help = 'Set transpaency of a matplotlib graph') 
    parser.add_argument("--theta", dest="theta", type=float, default=45, 
                        help = 'X-axis rotation for 3D view of waterfall graphs')    
    parser.add_argument("--gamma", dest="gamma", type=float, default=20, 
                        help = 'Z-axis rotation for 3D view of waterfall graphs')    
    parser.add_argument("--unit", dest="unit", type=str.lower, default='ev', 
                            help = 'Absorption spectra in eV or nm unit?, default in eV unit') 
    parser.add_argument("--figsize", dest="figsize",  nargs='+', type=int, default=(7,5), 
                            help = 'Plot figure size, integer list, default (7,5)')    
    
    
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
        Bands_all = outcar.band_energies()   # Output dimension is K points x frames x bands x 3 
        
        Bands = get_frame(Bands_all,nthkpoints, frames) 
        
        Fermi = outcar.fermi_energies()
        Fermi_frames = Fermi[-1,0] if frames =='final' else (Fermi[0,0] if frames=='initial' 
                                              else(Fermi[:,0] if frames=='all'
                                                   else (Fermi[int(frames),0] if frames.isdigit() 
                                                         else(Fermi[int(frames),0] if frames.lstrip('-').isdigit() 
                                                              else Fermi[:,0] ))))
        print(Fermi[0], 'Fermi zerooooooo')
        nhomo  = int(Bands[0,:,0][Bands[0,:,1] < Fermi_frames[0]][-1])  
        nbands = nbands if nbands < nhomo else nhomo-1 
        
        if options.plot == True:  
            print(Bands_all.shape)
            x = Bands_all[nthkpoints,-1,:,1] .astype(float)
           
            if options.shiftfermi: 
                x = x - outcar.fermi_energies()[-1,0]  
            if options.xlim: 
                xlim=options.xlim 
                #print(xlim)
                x = x[(x>xlim[0]) & (x<xlim[1]) ] 
            sigma = options.sigma
            x,y = Gauss_distr(x, sigma=sigma )  
            if not options.fout : 
                name = fname
            else: 
                name ='DOS_%s'% options.fout 
                
            plot_line(x,y, name = name, xlim=xlim, ylim=(0,),ls=None, lw=None, colors ='g', xticks=None, yticks=None, 
                          xticklabels=None, yticklabels=None,xlabel=r'E-$E_{fermi}$', ylabel=None)  
            np.savetxt('%s.dat'%name, np.c_[x,y],fmt='%0.4f  %0.4f' )    
        # Data write in a file options.fout, else print
    
        elif options.plot=='dos_all': 
            if not options.fout : 
                name = fname
            else: 
                name = options.fout
            if options.xlim: 
                xlim=options.xlim 
            else: 
                xlim=[Fermi_frames[0]-3,Fermi_frames[0]+3]   
            Bands_for_dos = Bands [:,:,1] 

            if options.shiftfermi: 
                Bands_for_dos =Bands_for_dos - Fermi_frames[:,None]
                xlim = xlim - Fermi_frames[0]
            print(xlim)
           # print(Bands_for_dos[0,:] >  ) 
            if options.xlim: 
                 xlim=options.xlim 
                 #print(xlim)
                 Bands_for_dos = Bands_for_dos.T[(Bands_for_dos[0,:] >xlim[0]) & (Bands_for_dos[0,:]<xlim[1]) ]
 #                Bands_for_dos = Bands_for_dos
            sigma = options.sigma            
            x,y = Gauss_distr(Bands_for_dos, sigma=sigma) 

            plot_waterfall(x,y,xlim=xlim,ylim=None, gamma = options.gamma, theta=options.theta, alpha=options.alpha,
                           color=None,name='DOS_all_%s'%name)
            np.savetxt('bands_all_%s.dat'%name,Bands_for_dos,fmt='%0.5f') 
            np.savetxt('Dress_dos_all_%s.dat'%name, np.c_[x,y], fmt='%0.4f') 

        elif options.fout != None: 
            Bands = Bands[:,nhomo-nbands-1:nhomo+nbands-1,:] 
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
            Bands = Bands[:,nhomo-nbands-1:nhomo+nbands-1,:]
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
            
    elif options.job == 'abs_gaus': 
        if options.fname: 
            fname  = options.fname 
        else: 
            print(parser.print_help()) 
            print('Provide gaussian TDDFT output file name') 
            sys.exit () 

        if options.fout: 
            fout = options.fout 
        else: 
            fout = os.path.splitext(os.path.basename(options.fname))[0]
        
        unit = options.unit 
        
        flog = open(fname)
        lines = flog.read() 
        Excited_state = re.findall(r'Excited State.*',lines) 
        Excited_state =np.array([[i.split()[2].strip(':'),i.split()[4],i.split()[6],i.split()[8].split('=')[-1]] for i in Excited_state]).astype(float) 

        if options.unit=='ev': 
            X,OS = Excited_state[:,1],Excited_state[:,-1] 
        if options.unit=='nm': 
            X,OS = Excited_state[:,2],Excited_state[:,-1] 

        XX,YY = Dress_abs(X,OS,linewidth=options.sigma, unit=options.unit) 
        np.savetxt('Absorption_%s_lw%s.dat'%(fout,options.sigma), np.c_[XX,YY], fmt='%0.4f  %0.4f' if unit=='ev' else '%0.2f  %0.4f') 

        if options.xlim: 
            print('User provided xlim', options.xlim) 
            x1,x2 = options.xlim 
        else:
            if unit=='nm': 
                x1 = 250
                x2 = 50*(XX[YY>YY.max()*0.001][0]//50+1) 
            else: 
                x1 = 0.0
                x2 = 1+int(XX[YY>YY.max()*0.001][-1]) 
        if not  all([abs(x1-x2)<50, unit=='ev']) | all([abs(x1-x2)>50, unit=='nm']):
            print('\n**Unit of absorption and xlim are not compatible**\n') 
            
        y2 = YY[(XX>x1) & (XX<=x2)].max()*1.2   

        ylabel = r'Absorption, $Lmol^{-1}cm^{-1}$'
        xlabel = 'Energy, %s' %('eV' if unit=='ev' else unit) 

        plt.rcParams.update({'font.size': 20})
        
        plot_line(XX,YY, xlim=[x1,x2], ylim=[0,y2], xlabel=xlabel, ylabel=ylabel, figsize=options.figsize, 
                  name='Absorption_%s_lw%s.png'%(fout,options.sigma))  

        
        
                    
        #    print(a.outcar2incar())
        #    a.md_traj(out='Temmp_md.dat')  
        # F
    

if __name__ == "__main__":
    main() 
    
    
    