

import numpy as np 
import math 
from sys import argv 
import os 

# This script can read Gaussian input file format, .com or .gif and 
# Convert it POSCAR file format to use in VASP  
# Input file should containe PBC information (Tv) 
# Following Variable can read - 
# Gaussian input file name: argument[1]  
# Cartesian or Farctional coodinates: cartesian=True/False or Direct=True/False ; Default: Direct=True  (optional)
# Select optimization Select=True/False; Default: True if -1/0 flag in inputfile 
# Output file name out=filename  (optional) 
##___________________________________________________________________
# Example syntax: COM2POSCAR.py Example1.com Select=True Direct=True  
#--------------------------------------------------------------------

def angle(a,b):
 import numpy as np 
 a=np.asarray(a) 
 b=np.asarray(b) 
 unitV1 = a/np.linalg.norm(a)
 unitV2 = b/np.linalg.norm(b)
 return  np.arccos(np.dot(unitV1, unitV2))

TrueFalse = lambda a: True if a.lower()=='true' else False 
remFragm = lambda a : a.split('(')[0] if 'fragment' in a.lower() else a 
FloatCheck = lambda a : all([i.replace('.','').isdigit() for i in a]) 

name = argv[1] 
if not os.path.isfile(name): 
	print("No input file found with name %s" %name )
	exit() 


f_outname = '%s.POSCAR'%name.split('\\')[-1].split('/')[-1].split('.')[0] 
file = open(name) 

line = file.readline() 
Select=False
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
Freeztag='' 
def Line2array(line): 
	
	global Freeztag 	
	ll = line.split()
	if 'tv' in line.split()[0].lower(): 
		return ['Tv']+line.split()[-3:] 
	elif len(ll)==4:
		Freeztag = False 
		if 'frag' in line.lower(): 
			return [ll[0].split('(')[0] ,'0', ll[1],ll[2],ll[3]]	
		else : 
			return [ll[0],'0',ll[1],ll[2],ll[3]]

	elif len(ll) == 5: 
		if ll[1].lstrip('-+').isdigit():
			Freeztag = True 
			return [ll[0].split('(')[0] , ll[1],ll[2],ll[3],ll[4]] 
			
		elif 'fragment' in ll[1].lower(): 
			if FloatCheck(ll[-3:]): 
				Freeztag = False 
				return(ll[0],'0',ll[-3],ll[-2],ll[-1])
			else: 
				print(line)
				print('Last three column should be X,Y and Z coodinates') 
				print('Script is terminating') 
				exit() 			
	elif len(ll) == 6: 
		if ll[1].lstrip('-+').isdigit():
			Freeztag = True
			if FloatCheck(ll[-3:]): 
				return [ll[0].split('(')[0] , ll[1],ll[-3],ll[-2],ll[-1]] 
		elif ll[2].lstrip('-+').isdigit():
			Freeztag = True
			if FloatCheck(ll[-3:]): 
				return [ll[0].split('(')[0] , ll[2],ll[-3],ll[-2],ll[-1]] 
		elif FloatCheck(ll[-3:]):
			Freeztag = False 
			return [ll[0].split('(')[0] , '0',ll[-3],ll[-2],ll[-1]]
		else: 
			print(line)
			print('Last three column should be X,Y and Z coodinates') 
			print('Script is terminating') 
			exit() 		
		
	else:
		print(line)
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
		elif (len(line.split())%2 ==0 ) & all([i.isdigit() for i in line.split()]) :
			line=file.readline() 
			while line: 
				if not len(line.strip()) < 1: 
					if not 'tv' in line.lower().split()[0]: 
						
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

#Checking Select tag 
if (Freeztag == False  ) & (Select == True): 
	print('No Freezing ion labels found.')  
	AA = input('All Ions are considering as "T T T " , Yes or No?:') 
	if AA.strip()[0].lower() != 'y': 
		print('Script is terminating, comeback with right file format.')
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

print(Coord)
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
				exit() 
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

