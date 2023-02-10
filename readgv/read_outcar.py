
import numpy as np 
import re 
import sys 
import warnings


class read_outcar:
    def __init__(self, filename='OUTCAR') :  
        self.filename = filename
    
        self.charge = {}
        self.energies = None
        self.last_energy = None
        self.forces = None
        self.fermi = None
        self.bands = None
        self.positions = None
        self.kpoints = None
        self.array_sizes = {}
        self.species = None
        self.final_data = {}
        self.iteration_data = []
        self.outcar = open(self.filename, 'r').read() 
        
        
        self.poscar = None 
        self.paw_ions = []      # Ions from the potential function order in Outcar 
        self.no_ions = []  # Number of ions per ion type 
        self.paw_rwigs_rad = [] #RWIGS from postcar informatio in outcar 
        self.ibrion = None 
        self.incar = {} 
        #self.outcar_distance_warning()
        self.outcar_parser() 
        
        
    def outcar_parser(self): 
        f = self.outcar 
        #Reading poscar name 
        poscar = re.findall('POSCAR.*=.*',f)[-1].split('=') [-1].strip() 
        self.poscar = poscar  
        
        # reading the IONS from POTCAR information in outcar file 
        paw_ = re.findall('TITEL.*',f)
        PAW_ = [i.split('=')[1].split()[0] for i in paw_] 
        PAW_ions =  [i.split('=')[1].split()[1] for i in paw_] 
        self.paw_ions = PAW_ions 
        
        no_ions = [int(i) for i in re.findall('ions per type.*',f)[0].split('=')[-1].split()]
        self.no_ions = no_ions 
        
        paw_rwigs_rad = [float(i.split('=')[-1].split()[0]) for i in re.findall('RWIGS.*RWIGS.*',f) ] 
        self.paw_rwigs_rad = paw_rwigs_rad          
        
        ibrion = int(re.findall(r'IBRION.*',f)[0].split('=')[1].split()[0]) 
        self.ibrion = ibrion
        
        tot_ions = int(re.findall(r' number of ions     NIONS.*',f)[0].split('=')[1].split()[0])
        self.tot_ions = tot_ions 
        
        nkpoints = int(re.findall(r'NKPTS .*',f)[0].split('=')[1].split()[0]) 
        self.nkpoints = nkpoints 
        
        self.n_iteration = int(re.findall(r"Iteration   .*", f)[-1][10:16]) 
        self.nbands = int(re.findall(r"NBANDS.*", f)[0].split('=')[1].split()[0] ) 
        ion_list = [] 
        for i,j in enumerate(no_ions):
            ion_list = ion_list+[PAW_ions[i]]*j  
        self.ion_list = ion_list   
        
    def outcar2incar(self): 
        f = self.outcar 
        Tags = ['PREC','ISTART','ICHARG','ISPIN','LSORBIT','METAGGA','ENCUT','ENINI','NELM','EDIFF',
                'LREAL','EDIFFG','NSW','IBRION','ISIF','POTIM','TEIN','TEBEG','TEEND','SMASS','RWIGS',
                'EMIN','EFERMI','ISMEAR','IALGO','LWAVE','LCHARG','LORBIT','LMONO','LDIPOL','GGA','AEXX','NBANDS'] 
        Incar_var= {} 
        for var in Tags: 
            tag_value_  = re.findall(fr"{var}.*", f)[0].split('=')[1].split()[0] 
            Incar_var[var] = tag_value_
        self.INCAR = Incar_var 
        return Incar_var        
            
    def outcar_distance_warning(self) :                
        # Checking bond distances if unusual bond length in input geometry 
        # Getting the ions coordinates from outcar 
        f =  open(self.filename,'r')
        line = f.readline() 
        ion_pos = []
        ion_distance = [] 
        ion_distancewith = [] 
        while line: 
            if 'nearest neighbor table' in line: 
                line = f.readline() 
                while line:
                    try: 
                        kk,ll = line.split('-') 
                        n_ion = int(kk.split()[0]) 
                                                 
                        ion_pos.append([float(i) for i in kk.split()[-3:]]) 
                        k = int(len(ll.split())/2) 
                        ion_distance.append([float(ll.split()[1+i*2]) for i in range(k)])
                        ion_distancewith.append([int(ll.split()[i*2]) for i in range(k)])
                          #  print(line)
                        line = f.readline() 
                        if n_ion == sum(self.no_ions) : 
                            line = False  
                            break 
                    except ValueError as error: 
                        print(error.args[0])
                        if ('base 10' in error.args[0]) or ('to float' in error.args[0]): 
                            print(line)
                            print(error)
                            warnings.warn("Ions distance check: Distance can't convert to float in above line in OUTCAR file!" )
                            line = False 
                        else: 
                            line = f.readline() 
                            
                    except IndexError:  
                        line = f.readline() 
                break 
            else : 
                line = f.readline()            
        
        Corr_geom=True 
        for i,j in enumerate(ion_pos): 
            if any([l<0.10 for l in ion_distance[i]]): 
                print('Following ions are found in too close distance, \n check the input geometry') 
                print('Atom', j, ' and ', ion_distancewith[i], ' distance' , ion_distance[i] ) 
                Corr_geom=False 
                f.close() 
                sys.exit() 
        if Corr_geom==True: 
            print('No error found in initial coordinates')
        #f.close() 
    
    def input_geom(self, out=False,fout=None, **kwargs):  
        if fout==None: 
            fout = '%s_ini.xyz'%self.poscar
        f =  self.outcar   
        var =  "position of ions in fractional coordinates.*\n((.*\n){%i})"%self.tot_ions 
        frac_coordinates = re.findall(fr"{var}", f)[0][0].split('\n')  
        frac_coordinates = [i for i in frac_coordinates if len(i.strip()) != 0] 
                
        var =  "position of ions in cartesian coordinates.*\n((.*\n){%i})"%self.tot_ions 
        cart_coordinates = re.findall(fr"{var}", f)[0][0].split('\n')  
        cart_coordinates = [i for i in cart_coordinates if len(i.strip()) != 0]    

        if (len(kwargs)!=0 ) & ('coord' not in kwargs): 
            print('use kwargs coord=farc or cart to print coordinates')  
            return None 
        elif len(kwargs)==0: 
            kwargs['coord'] = 'cart' 
        if  kwargs['coord'].lower()[:4] == 'frac' :  
            CC = frac_coordinates 
        else : 
            CC = cart_coordinates 
        if len(CC) != self.tot_ions: 
            print("Fractional coordinated can't read properly, \n please check the output file fomrat")
            return None  
        else : 
            CC = ['{:<3}'.format(self.ion_list[i])+CC[i] for i in range(self.tot_ions)] 
            if out is True:  # Destination: user output, Flase : return a variable for other function 
                ff = open(fout,'w')
                ff.write(' %s\n %s \n'%(len(CC),self.poscar)) 
                for line in CC: 
                    ff.write('{:<3}{:>16}{:>16}{:>16} \n'.format(line.split()[0],line.split()[1],line.split()[2],line.split()[3]))
#                    ff.write(line+'\n') 
            else: 
#            print('Following is the input geometry in fractional coordinates') 
                return CC
        
    def md_traj(self,out=False, fout=None):
        if fout==None: 
            fout = '%s_traj.xyz'%self.poscar 
        f = self.outcar 
        n_iteration = self.n_iteration  # = int(re.findall(r"Iteration   .*", f)[-1][10:16])
        Trajectory = np.zeros((n_iteration, self.tot_ions,3))
        # Regex itreation is quite slow, need to change it 
        var =  "POSITION                                       TOTAL-FORCE.*\n.*\n((.*\n){%i})"%self.tot_ions
        geom_i = re.findall(fr"{var}", f)
#        geom_i = re.findall(r'POSITION                                       TOTAL-FORCE.*\n.*\n((.*\n){124})',f) 
        for j,i in enumerate(geom_i): 
            geom = np.loadtxt(i[0].split('\n')[:-1])[:,:3] 
            Trajectory[j,:,:] = geom          
        if len(geom_i) != n_iteration: 
            print('All ionic iteractions aren\'t finished completely')
        # Writing the file in xyz format 
        #print(out)
        if out != False: 
            traj_out = open(fout,'w') 
            for i in range(n_iteration): 
                traj_out.write('{} \nFrame {} \n'.format(self.tot_ions,i+1))
                for j in range(self.tot_ions): 
                    STR = ' {:<6}{:>12}{:>12}{:>12}'.format(self.ion_list[j] , Trajectory[i,j,0], Trajectory[i,j,1], Trajectory[i,j,2]) 
                    traj_out.write(STR+'\n')
                traj_out.write('\n') 
            traj_out.close()  
        else: 
            return Trajectory 
    
    def kpoint_dimension(self): 
        f=self.outcar 
        var = "k-points in units of 2pi.*\n((?:.*\n){%i})"%self.nkpoints 
        kpoints_dim_2pi = np.loadtxt( re.findall(fr"{var}", f).split('\n'),dtype=float)
        var = "k-points in reciprocal lattice and weights.*\n((?:.*\n){%i})"%self.nkpoints 
        kpoints_dim_recipr = np.loadtxt( re.findall(fr"{var}", f).split('\n'),dtype=float)  
        return kpoints_dim_2pi, kpoints_dim_recipr 
        
    def fermi_energies(self): 
        f=self.outcar
        fermi_all = re.findall(r'E-fermi : .*',f)  
        fermi_alphabeta = np.array([[float(i.split(':')[1].split()[0]),float(i.split(':')[3])] for i in fermi_all])
        return fermi_alphabeta 
    
    def band_energies(self) :
        f=self.outcar 
        States = np.zeros((self.nkpoints, self.n_iteration, self.nbands,3))
        
        for k in range(self.nkpoints): 
            var =  "k-point     %i.*\n.*band No.  band energies     occupation .*\n((?:.*\n){%i})"%(k+1,self.nbands) 
            bands  = re.findall(fr"{var}", f) 
            bands = np.array([np.loadtxt(i.split('\n')[:self.nbands] ) for i in bands ] )
            States[k,:,:,:] = bands 
#            bands = np.array([np.loadtxt(i[0].split('\n')[:nbands] ) for i in bands ] ) 
      #  j.split() for i in bands for j in i[0].split('\n')[:nbands]  
        return States 
       
    def free_energy(self): 
        f=self.outcar 
        free_en = np.array([i.split() for i in re.findall(r'free  energy.*', f)])[:,-2].astype(float) 
        return free_en 

class read_procar: 
    def __init__(self, filename='PROCAR') :  
        self.filename = filename
    

            
#if __name__ == '__main__': 
#    f = 'Datafiles/OUTCAR' 
#    a = Output(filename=f)
#    print(a.outcar2incar())
#    a.md_traj(out='Temmp_md.dat')  
#  'NKPTS', 'NBANDS', 'NEDOS', 'NIONS', 'NGX', 'NGY', 'NGZ', 'NGXF', 'NGYF', 'NGZF', 'ISPIN']: 