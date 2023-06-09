#!//home/mohammed.jabed/bin/anaconda3/bin/python


######################################################################################
## Script to convert VASP-5.0 CONTCAR/POSCAR file to xyz file format (default) 		$$
## Or, gaussian input file (.com) 													$$
## The first argument is a input file or CONTCAR is default (if available) 			$$
## To get gaussian input file format, add argument 'com' 							$$
## To specify output file name, optional atgument:  out=filename					%%
## Example run: python CONTCAR2XYZ.py POSCAR com out=filename  						%%
##  			Mohammed A Jabed (jabed.abu@gmail.com) 								%%
##$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

import numpy as np
from sys import argv, platform  
import sys 
import os , select   

#Reading contcar file, if exist, else read the CONTCAR type file name given as a frist argument 
try: 
	name = argv[1] 
	if not os.path.isfile(name):  
		print('File not found: %s' %name) 
		name=input("Input file name to convert:  ").strip() 
except: 
	if 'lin' in platform: 
		print("What is the input file name to convert (30s waiting time):")
		i, o, e = select.select( [sys.stdin], [], [], 30 )
		if len(i.stript())==0: 
			print('No File name is given.') 
			if os.path.isfile('CONTCAR'):
				print('CONTCAR file is found, using the CONTCAR file')
				name = 'CONTCAR' 
		else:
			name = i.strip() 
	elif 'win' in platform: 
		name=input("Input file name to convert:  ").strip()  
		if len(name.strip())==0: 
			if os.path.isfile('CONTCAR'):
				print('CONTCAR file is found, using the CONTCAR file')
				name = 'CONTCAR' 
			
if os.path.isfile(name):
	f_CONT=open(name,'r') 
	name_out=name.split('.')[0]
else: 
	print('File not found: %s' %name) 
	exit() 

Gaus_out=False 
for i in argv:
	if i == "com":  
		Gaus_out=True 
	if i.startswith('out='): 
		name_out = i.split('=')[-1].rstrip('.com').rstrip('.xyz')


# reading file, and cell normalizing factor  
line_CONT = f_CONT.readlines()  
compname = line_CONT[0] 
sigma = float(line_CONT[1].split()[0]) 
#print(line_CONT[-5:])
#Checking file format 
if not all([len(i.split())==3 for i in line_CONT[2:5]]): 
	print(line_CONT[2:5])
	print('File format does\'t match with CONTCAR file format, script is quiting')
	exit() 
#Reading unit cell dimention 
XYZ = np.loadtxt(line_CONT[2:5], dtype=float)
print('The unit cell size is: ')
print(XYZ)

Elements = line_CONT[5].split() 	#Elements list 
El_num = [int(i) for i in line_CONT[6].split()]    #number of atoms of each elements 
if line_CONT[7].strip()[0].lower()=='s': 
	selec_dynm = True
	if line_CONT[8].strip()[0].lower()=='d': 
		FracCoord=True 
	elif line_CONT[8].strip()[0].lower()=='c':
		FracCoord=False
elif line_CONT[7].strip()[0].lower()=='d': 
	FracCoord=True
	selec_dynm = False
elif line_CONT[7].strip()[0].lower()=='c':
	FracCoord=False
	selec_dynm = False
else: 
	"Can\'t read cartesian or direct coordinate correctly, line number 8 and 9. \nScript is quiting now!! :( "
	exit() 

if selec_dynm is True:  
	Cord=np.loadtxt(line_CONT[9:9+sum(El_num)], dtype=str)[:,:3].astype(float)
	AtomFrez = (np.loadtxt(line_CONT[9:9+sum(El_num)], dtype=str)[:,3:]=='F').all(axis=1)
	Freez =  AtomFrez*-1
else:
	try: 
		Cord = np.loadtxt(line_CONT[8:8+sum(El_num)], dtype=str).astype(float) 
	except ValueError:  
		print('Trying to read Direct coordinate format, found string in line, \nScript is quiting') 
		exit() 

 
Atom_symbol =[] 
for i in range(len(Elements)): 
 Atom_symbol = Atom_symbol + [Elements[i]]*El_num[i] 

if FracCoord is True:
 Cart_Cord = np.dot(XYZ,Cord.T).T
elif FracCoord is False:
 Cart_Cord = Cord 


if selec_dynm is True: 
 AA = np.c_[np.asarray(Atom_symbol),Freez]  
 Cord_Str = np.c_[AA,np.asarray(Cart_Cord).round(4)] 
elif selec_dynm is False: 
 Cord_Str = np.c_[np.asarray(Atom_symbol),np.asarray(Cart_Cord).round(4)] 
print(name_out)
if Gaus_out == True: 
 Header = '''%%mem=20GB
%%nprocshared=40
%%chk=CONTCAR.chk 

#p opt pbe1pbe/gen nosymm pseudo=read scf=maxcycles=10000

%s

0 1
''' %compname.rstrip('\n')
 
 if os.path.isfile('%s.com'%name_out):
  print ('%s.com is exist'%name_out) 
  name = input('What should be the output file name?:').replace(" ", "").rstrip('.com')  
  if not len(name)==0: 
   name_out= name     
  if os.path.isfile('%s.com'%name_out): 
   os.remove('%s.com'%name_out) 
   print('Overwriting the file %s.com'%name_out)
 with open('%s.com'%name_out, "w") as f:
  f.write(Header)
  if selec_dynm is True: 
   np.savetxt(f, Cord_Str,fmt='%-5.10s %-5.10s %-10.10s  %-10.10s  %-10.10s')
  elif selec_dynm is False: 
   np.savetxt(f, Cord_Str,fmt='%-5.10s  %-10.10s  %-10.10s  %-10.10s')
  PBC = np.c_[np.asarray(['Tv','Tv','Tv']),XYZ.astype('<U10')] 
  np.savetxt(f, PBC,fmt='%-5.10s  %-10.10s  %-10.10s  %-10.10s')
 print('Gaussian input file is written in the file name %s.com' %name_out) 

else :
 if os.path.isfile('%s.xyz'%name_out):
  print ('%s.xyz is exist'%name_out) 
  name = input('What should be the output file name?:').replace(" ", "").rstrip('.xyz')  
  if not len(name)==0: 
   name_out= name     
  if os.path.isfile('%s.xyz'%name_out): 
   os.remove('%s.xyz'%name_out) 
   print('Overwriting the file %s.xyz'%name_out)
 with open('%s.xyz'%name_out.rstrip('.xyz'), "a") as f_out: 
  print(len(Cord_Str))
  f_out.write('%i \n' %len(Cord_Str))
  f_out.write('%s \n' %compname.rstrip('\n')) 
  np.savetxt(f_out,np.c_[Cord_Str[:,0],Cord_Str[:,-3:]],header="%.0f \n"%sum(El_num),comments="",
			fmt='%-5.10s  %-10.10s  %-10.10s  %-10.10s') 
 print('xyz file is written in the file name %s.xyz' %name_out) 
print('Happy Calculation') 
