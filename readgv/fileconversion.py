# -*- coding: utf-8 -*-
"""
Created on Fri Dec 23 21:33:16 2022

@author: Mohammed a jabed, jabed.abu@gmail.com
"""
import numpy as np 
import math 
import os 

def angle(a,b):
    import numpy as np 
    a=np.asarray(a) 
    b=np.asarray(b) 
    unitV1 = a/np.linalg.norm(a)
    unitV2 = b/np.linalg.norm(b)
    return  np.arccos(np.dot(unitV1, unitV2))



def contcar2com(name, name_out,Gaus_out=False):  # name is the contcar file name 
    f_cont = open(name,'r') 
    # reading file, and cell normalizing factor  
    line_cont = f_cont.readlines()  
    sys_name  = line_cont[0] #compname  
    sigma = float(line_cont[1].split()[0]) 

    #Checking file format, cell dimension check 
    if not all([len(i.split())==3 for i in line_cont[2:5]]): 
        print(line_cont[2:5])
        print('File format does\'t match with CONTCAR file format, script is quiting')
        exit() 

        #Reading unit cell dimention 
    xyz = np.loadtxt(line_cont[2:5], dtype=float)
    print('The unit cell size is: ')
    print(xyz)

    Elements = line_cont[5].split() 	#Elements list 
    El_num = [int(i) for i in line_cont[6].split()]    #number of atoms of each elements 
    if line_cont[7].strip()[0].lower()=='s': 
        selec_dynm = True
        if line_cont[8].strip()[0].lower()=='d': 
            FracCoord=True 
        elif line_cont[8].strip()[0].lower()=='c':
            FracCoord=False
    elif line_cont[7].strip()[0].lower()=='d': 
        FracCoord=True
        selec_dynm = False
    elif line_cont[7].strip()[0].lower()=='c':
        FracCoord=False
        selec_dynm = False
    else: 
        "Can\'t read cartesian or direct coordinate correctly, line number 8 and 9. \nScript is quiting now!! :( "
        exit() 

    if selec_dynm is True:  
        Cord=np.loadtxt(line_cont[9:9+sum(El_num)], dtype=str)[:,:3].astype(float)
        AtomFrez = (np.loadtxt(line_cont[9:9+sum(El_num)], dtype=str)[:,3:]=='F').all(axis=1)
        Freez =  AtomFrez*-1
    else:
        try: 
            Cord = np.loadtxt(line_cont[8:8+sum(El_num)], dtype=str).astype(float) 
        except ValueError:  
            print('Trying to read Direct coordinate format, found string in line, \nScript is quiting') 
            exit() 

 
    ions_list = [] # Atom_symbol =[] 
    for i in range(len(Elements)): 
        ions_list = ions_list + [Elements[i]]*El_num[i] 

    #Convertion matices 
    a1,a2,a3 = xyz[0,:]
    b1,b2,b3 = xyz[1,:]
    c1,c2,c3 = xyz[2,:] 
    #Length of the Parallelepiped edge 
    a,b,c= np.linalg.norm(xyz[0,:]), np.linalg.norm(xyz[1,:]), np.linalg.norm(xyz[2,:])
    # angel of Parallelepiped 
    alpha = angle([b1,b2,b3],[c1,c2,c3])  
    beta  = angle([a1,a2,a3],[c1,c2,c3]) 
    gamma  = angle([a1,a2,a3],[b1,b2,b3]) 
    #Simplified term. look on the Wikipedia for more 
    n = (math.cos(alpha)-math.cos(gamma)*math.cos(beta))/math.sin(gamma) 
    Trans_matrix_inv = np.array([[a,  b*math.cos(gamma), c*math.cos(beta),],
                                 [0,  b*math.sin(gamma), c*n],
                                 [0,  0,                 c*((math.sin(beta))**2-n**2)**0.5]])  

    if FracCoord is True:
        Cart_Cord = np.matmul(Trans_matrix_inv,Cord.T).T
    elif FracCoord is False:
        Cart_Cord = Cord 

    if selec_dynm is True: 
        AA = np.c_[np.asarray(ions_list),Freez]  
        Cord_Str = np.c_[AA,np.asarray(Cart_Cord).round(4)] 
       # print(np.asarray(Cart_Cord).round(4))
    elif selec_dynm is False: 
        Cord_Str = np.c_[np.asarray(ions_list),np.asarray(Cart_Cord).round(4)] 
      #  Cord_Str = np.c_[np.asarray(ions_list),np.asarray(np.array_str(np.asarray(Cart_Cord),precision=4))] 
   
    print(name_out)
    if Gaus_out == True: 
        Header = '''%%mem=20GB
%%nprocshared=40
%%chk=%s.chk 

#p opt pbe1pbe/gen nosymm pseudo=read scf=maxcycles=10000

%s

0 1
''' %(name_out,sys_name.rstrip('\n')) 

        with open('%s.com'%name_out, "w") as fout:
            fout.write(Header)
   #         print(Cord_Str)
            if selec_dynm is True: 
                np.savetxt(fout, Cord_Str,fmt='%-5.10s %-5.10s %-10.10s  %-10.10s  %-10.10s') 
                
            elif selec_dynm is False: 
                np.savetxt(fout, Cord_Str,fmt='%-5.10s  %-10.6s  %-10.10s  %-10.10s')
               
            PBC = np.c_[np.asarray(['Tv','Tv','Tv']),xyz.astype('<U10')] 
            np.savetxt(fout, PBC,fmt='%-5.10s  %-010.8s  %-10.8s  %-10.10s') 
             
        print('Gaussian input file is written in the file name %s.com' %name_out)  
        fout.close() 


    else :
        with open('%s.xyz'%name_out.rstrip('.xyz'), "a") as fout: 
            #print(len(Cord_Str))
            fout.write('%i \n' %len(Cord_Str))
            fout.write('%s \n' %sys_name.rstrip('\n')) 
            np.savetxt(fout,np.c_[Cord_Str[:,0],Cord_Str[:,-3:]],
                       fmt='%-5.10s  %-10.10s  %-10.10s  %-10.10s') 
        print('xyz file is written in the file name %s.xyz' %name_out) 
        fout.close() 


def com2poscar(name,name_out=True, select=False, direct=True): 
    def angle(a,b):
        a=np.asarray(a) 
        b=np.asarray(b) 
        unitV1 = a/np.linalg.norm(a)
        unitV2 = b/np.linalg.norm(b)
        return  np.arccos(np.dot(unitV1, unitV2))
    
    file = open(name) 
    line = file.readline() 
    topheader = [] 
    commandkeys = [] 
     
    while line : 
        print('%' in line)
        if '%' in line: 
#            print(line)
            topheader.append(line)
            line = file.readline() 
#            print('command')
        else: 
            break 
    if len(line.strip())==0: 
        line = file.readline() 
    while line : 
        if '#' in line: 
            commandkeys.append(line)
            line=file.readline() 
#            print('keys')
        else: 
            break 
    if len(line.strip())==0: 
        line=file.readline() 
        
    comp_name = line.strip('\n') 
    line = file.readline() 
#    print('name')
        
    if len(line.strip())==0: 
        line=file.readline()  
    Coord = []
    TVs= [] 
    while line: 
        if (len(line.split())%2 ==0 ) & all([i.isdigit() for i in line.split()]) :
#            print('charg_multi')
            chrg_multi = line.strip('\n') 
            line = file.readline() 
            while line :
                if not 'tv' in line.lower(): 
                    Coord.append(line.strip('\n') )
                else: 
                    TVs.append(line.strip('\n').split()) 
                line = file.readline()
                if len(line.strip())==0: 
                    break
            break 
        
    # Warning 
    if len(TVs) <3: 
        print('Periodic cell dimensions not found.') 
        return None  
    
    freeztag = lambda a : False if a==-1 else True 
    ionsFrag = lambda a : a.split('(')[0].strip() if 'frag' in a.lower() else a.strip() 
    
    coordinates_xyz = np.empty((len(Coord),5), dtype=object) 
    for j,coord_ in enumerate(Coord):
        ll = coord_.split() 
        coordinates_xyz[j,0]  = ionsFrag( ll[0])   
        coordinates_xyz[j,2:5] = ll[-3:] 
        if len(ll) ==4: 
            coordinates_xyz[j,1] = True  
        elif len(ll) ==5: 
            coordinates_xyz[j,1] = freeztag(int(ll[1]))  
   
    coordinates_xyz = coordinates_xyz[coordinates_xyz[:,0].argsort()] 
    
    
    
    #working on Cartesian Fractional coordinates conversion 
    print(np.asarray(TVs)) 
    TVV = np.asarray(TVs)[:,1:].astype(float) 
    a1,a2,a3 = TVV[0,:]
    b1,b2,b3 = TVV[1,:]
    c1,c2,c3 = TVV[2,:] 
    #Length of the Parallelepiped edge 
    a,b,c= np.linalg.norm(TVV[0,:]), np.linalg.norm(TVV[1,:]), np.linalg.norm(TVV[2,:])
    # angel of Parallelepiped 
    alpha = angle([b1,b2,b3],[c1,c2,c3])  
    beta  = angle([a1,a2,a3],[c1,c2,c3]) 
    gamma  = angle([a1,a2,a3],[b1,b2,b3]) 
    #Simplified term. look on the Wikipedia for more 
    n = (math.cos(alpha)-math.cos(gamma)*math.cos(beta))/math.sin(gamma) 
    Trans_matrix_inv = np.array([[a,  b*math.cos(gamma), c*math.cos(beta),],
                                 [0,  b*math.sin(gamma), c*n],
                                 [0,  0,                 c*((math.sin(beta))**2-n**2)**0.5]]) 
    Trans_matrix = np.linalg.inv(Trans_matrix_inv)
#    delta = a*b*c*(1-math.cos(alpha)**2-math.cos(beta)**2-math.cos(gamma)**2+
 #                  2*math.cos(alpha)*math.cos(beta)*math.cos(gamma))**0.5 
    #Convertion matices 
#    Trans_matrix = np.array([[1/a, -1/(math.tan(gamma)*a), b*c*(math.cos(alpha)*math.cos(gamma)-math.cos(beta))/(delta*math.sin(gamma)) ],
#                          [0, 1/(b*math.sin(gamma)), a*c*(math.cos(beta)*math.cos(gamma) - math.cos(alpha))/(delta*math.sin(gamma))],
#                          [0,0, a*b*math.sin(gamma)/delta]
 #                         ])
     
    frac_coord = np.matmul(Trans_matrix,coordinates_xyz[:,2:].astype(float).T).T 
    ions_name = list(set(coordinates_xyz[:,0]))
    ions_name.sort() 
    print(ions_name)
    ions_count = [] 
    for j in ions_name: 
        ions_count.append((coordinates_xyz[:,0] == j).sum())

    if name_out==None:
        print('Parallelepiped vectors are - ') 
        print(TVV)
        print('Ions are', ions_name) 
        print('Ions counts are: ',ions_count) 
        return frac_coord
    else: 
        name_ = os.path.splitext(name)[0] 
        name_ = os.path.basename(name_)
        if name_out == True:            
            fout = open('{}.POSCAR'.format(name_),'w') 
        else: 
            fout = open('{}'.format(name_out),'w')  
        
        fout.write(name_+'\n')
        fout.write(' 1.000000  \n') 
        for j in range(3): 
        	fout.write('  %0.6f   %0.6f    %0.6f  \n' %(TVV[j,0],TVV[j,1],TVV[j,2])) 
        
        # Ions types 
        for i in ions_name: 
        	fout.write(' %s '%i) 
        fout.write('\n') 
        #Number of each Ions type 
        for i in ions_count: 
        	fout.write(' %i ' %i) 
        fout.write('\n') 
        
        #Writing the Coordinates
        if select: 
            fout.write('select \n') 
        if direct: 
            fout.write('Direct \n') 
        else: 
            fout.write('Cartesian \n') 
        
        for j in range(coordinates_xyz.shape[0]): 
            if direct == False: 
                new_xyz = np.asarray(coordinates_xyz[j,2:],dtype=float ) 
            else: 
                new_xyz = frac_coord[j,:] 
                
            if select:
                if coordinates_xyz[j,1] == True:
                    sel_str = '  T      T       T \n'
                elif coordinates_xyz[j,1] == False: 
                    sel_str = '  F      F       F  \n'
            else:
                sel_str = '  \n'
#            fout.write('%0.6f   %0.6f   %0.6f   %s' %(new_xyz[0],new_xyz[1],new_xyz[2],sel_str)) 
            fout.write('{:> 12.7f}   {:> 12.7f}   {:> 12.7f}   {}'.format(new_xyz[0],new_xyz[1],new_xyz[2],sel_str))
        print('POSCAR file written in the file %s' %name_out)
        fout.close() 
            



