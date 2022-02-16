
from telnetlib import SE
import numpy as np 
import math 
from sys import argv 
import os 

# This file will read Gaussian input file format, .com or .gif 
# Convert it too POSCAR file format 
# Input file should containe PBC information (Tv) 

def angle(a,b):
 import numpy as np 
 import math 
 a=np.asarray(a) 
 b=np.asarray(b) 
 unitV1 = a/np.linalg.norm(a)
 unitV2 = b/np.linalg.norm(b)
 return  np.arccos(np.dot(unitV1, unitV2))

TrueFalse = lambda a: True if a.lower()=='true' else False 
remFragm = lambda a : a.split('(')[0] if 'fragment' in a.lower() else a 

name = argv[1] 
f_outname = '%s.POSCAR'%name.split('.')[0] 
file = open(name) 

line = file.readline() 
Select=True
Direct = True 
for i in argv: 
	if i.lower().startswith('select'): 
		Select = TrueFalse(i.split('=')[-1]) 
	if i.lower().startswith('car'): 
		print(i)
		Direct = [False if TrueFalse(i.split('=')[-1]) else True][0]
	if i.lower().startswith('dir'): 
		Direct = TrueFalse(i.split('=')[-1])
	if i.lower().startswith('out'): 
		f_outname = i.split('=')[-1] 

print(Direct)
def Line2array(line): 
	global A, Select 	
	ll = line.split()
	if 'tv' in line.split()[0].lower(): 
		return ['Tv']+line.split()[-3:] 
	elif len(line.split())==4: 
		if (A != 'y' ) & (Select == True): 
			AA = input('All Ions is writting as "T T T " , Yes or No?:') 
			if AA.lower().strip()[0] == 'y':
				A = 'y'
				Select = True 
			else: 
				Select =False
		if A == 'y':
			Ion = remFragm(ll[0]) 
			return [Ion,'0',ll[1],ll[2],ll[3]] 
		else: 
			if Select == False: 
				return [remFragm(ll[0])]+ll[-3:]
			else: 
				print('Script is exiting, provide a correctly formated *.com file') 
				exit() 

	elif len(line.split())==5: 
		if line.split()[1].lstrip('-+').isdigit():
			Ion=remFragm(ll[0])  
			print([Ion]+ll[1:] )
			return [Ion]+ll[1:] 
		else: 
			if 'fragment' in ll[1].lower(): 
				A = input('All Ions are written as "T T T " , Yes or No?:') 
				if A.lower().strip()[0] == 'y': 
					Ion = remFragm(ll[0]) 
					return [Ion,'0',ll[1],ll[2],ll[3]] 
				else: 
					print('Script is exiting, Check the column no. 2 and 3 in the *.com file') 
					exit() 
	else: 
		print('Please provide a gaussian input file format as follows - \n \
[Ions	0/-1	X	Y	Z] \n \
			or \n \
[Ions	X	Y	Z] \n')  
		print('Script is terminating')
		exit() 


A=''
Coord = [] 
TV = [] 
while line: 
	if len(line.strip()) != 0: 
		if line.strip()[0] == '\%': 
			line=file.readline() 
			pass 
		elif line.strip()[0] == '#': 
			line=file.readline() 
			if len(line.strip()) ==0: 
				line = file.readline() 
				mod_name = line				                
			pass
		elif (len(line.split()) ==2 ) & all([i.isdigit() for i in line.split()]) :
			line=file.readline() 
			while line: 
				if not len(line.strip()) < 1: 
					if not 'tv' in line.lower().split()[0]: 
						#print(line)
						Coord.append( Line2array(line)) 
						#print(Line2array(line))
						try: 
							line = file.readline() 
						except: 
							break 
					else :  
						#print(line)
						TV.append(Line2array(line))  
						try: 
							line = file.readline() 
						except: 
							break 
				else: 
					file.close() 
					break 
		else: 
			try: 
				line=file.readline() 
			except: 
				break 
	else:
		try: 
			line=file.readline() 
		except: 
			break 

#Checking if PBC information is available 
if len(TV) ==0: 
    print('Can\'t find the PBC lattice dimention in the filr %s ' %name )
    print('No file is written')
    exit() 
    
#Writting the output file with Fractional coordinates  

if os.path.isfile(f_outname):
    name_=input('{} file exist, what should be the output file name? '.format(f_outname)).strip() 
    if len(name_.strip()) > 1: 
        f_outname = name_
		
f_out=open(f_outname,'w')  
f_out.write(mod_name) 
f_out.write(' 1.000000  \n') 

#Writing the PBC vector 
#print(TV)
for i in TV: 
	f_out.write('  %0.6f   %0.6f    %0.6f  \n' %(float(i[1]),float(i[2]),float(i[3]))) 


Coord = np.asarray(Coord) 
# Ions types 
Ele = np.sort(np.unique(Coord[:,0]))
for i in Ele: 
	f_out.write(' %s '%i) 
f_out.write('\n') 

#Number of each Ions type 
for i in Ele: 
	A = Coord[Coord[:,0]==i]
	f_out.write(' %s ' %A.shape[0])

#Whether selective coordinates or not 
print('Select: %s'%Select)	 
f_out.write('\n')
if Select == True: 
	f_out.write('select \n') 

#Direct coordinates 
if Direct == True: 
	f_out.write('Direct \n') 
else :
	f_out.write('Cartesian \n') 

#working on Cartesian Fractional coordinates conversion 
TVV = np.asarray(TV)[:,1:].astype(float) 
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
delta = a*b*c*(1-math.cos(alpha)**2-math.cos(beta)**2-math.cos(gamma)**2+2*math.cos(alpha)*math.cos(beta)*math.cos(gamma))**0.5 

#Convertion matices 
Conv_matr = np.array([[1/a, -1/(math.tan(gamma)*a), b*c*(math.cos(alpha)*math.cos(gamma)-math.cos(beta))/(delta*math.sin(gamma)) ], 
[0, 1/(b*math.sin(gamma)), a*c*(math.cos(beta)*math.cos(gamma) - math.cos(alpha))/(delta*math.sin(gamma))], 
[0,0, a*b*math.sin(gamma)/delta]
])
 
#Writing the Coordinates
for element in Ele: 
    A = Coord[Coord[:,0] == element] 
    print("Number of {:>2} ions is {:>4}".format(element,A.shape[0]))
    for j in range(A.shape[0]): 
        if Select == True: 
            if A[j,1] == '0': 
                sel_str = '  T      T       T \n' 
            elif A[j,1] == '-1': 
                sel_str = '  F      F       F  \n'
            else: 
                print('Working on coordinates line %s' %A[j,:]) 
                print('Selective coordinates are not found freezing coordinates in Gaussview format file %s' %name) 
                print('program is terminating')
        else: 
            sel_str = '  \n'
        xyz = np.asarray(A[j,-3:],dtype=float) 
        if Direct == False: 
            new_xyz = xyz 
        else: 
            new_xyz = np.matmul(xyz,Conv_matr) 
        f_out.write('%0.6f   %0.6f   %0.6f   %s' %(new_xyz[0],new_xyz[1],new_xyz[2],sel_str)) 
f_out.write('\n')
print('POSCAR file written in the file %s' %f_outname)