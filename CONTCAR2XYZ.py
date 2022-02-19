#!//home/mohammed.jabed/bin/anaconda3/bin/python


####################################################
##
## Script to convert VASP-5.0 CONTCAR/POSCAR file to XYZ (default) and 
## gaussian input file (.com) 
## arg 1 should be 'com' to get .com file. 
## 
##$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

import numpy as np
from sys import argv 
import os  
 
#Reading contcar file, if exist, else read the CONTCAR type file name  
if len(argv) < 2: 
 try: 
  name = 'CONTCAR'
  f_CONT=open(name,'r') 
  name_out='CONTCAR'
 except: 
  name=input("No CONTCAR file in the directory, file name:  ") 
  f_CONT=open(name,'r') 
  name_out=name.split('.')[0]

if len(argv) >= 2: 
# try:
 name = argv[1]
 f_CONT = open(name,'r') 
 name_out=name.split('.')[0]
# except: 
#  name = argv[2]
#  f_CONT=open(name,'r')

Gaus_out=False 
if "com" in argv: 
 Gaus_out=True  

# reading file, and cell normalizing factor  
line_CONT = f_CONT.readlines()  
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

if Gaus_out == True: 
 Header = '''%%mem=20GB
%%nprocshared=40
%%chk=CONTCAR.chk 

#p opt pbe1pbe/gen nosymm pseudo=read scf=maxcycles=10000

%s

0 1
''' %name
 
 if os.path.isfile('%s.com'%name_out):
  print ('%s.com is exist'%name_out) 
  name = input('What should be the output file name?:').replace(" ", "")  
  try:
   os.remove(name) 
  except FileNotFoundError: 
   pass
  if len(name)==0: 
   name= '%s.com'%name_out
   print('Overwriting the file %s'%name)   
 else:
  name = '%s.com'%name_out   
 if '.' not in name: 
  name='%s.com'%name
 with open(name, "w") as f:
  f.write(Header)
  if selec_dynm is True: 
   np.savetxt(f, Cord_Str,fmt='%-5.10s %-5.10s %-10.10s  %-10.10s  %-10.10s')
  elif selec_dynm is False: 
   np.savetxt(f, Cord_Str,fmt='%-5.10s  %-10.10s  %-10.10s  %-10.10s')
  PBC = np.c_[np.asarray(['Tv','Tv','Tv']),XYZ.astype('<U10')] 
  np.savetxt(f, PBC,fmt='%-5.10s  %-10.10s  %-10.10s  %-10.10s')
else :  
  np.savetxt('%s.xyz'%name_out,np.c_[Cord_Str[:,0],Cord_Str[:,-3:]],header="%.0f \n"%sum(El_num),comments="",fmt='%-5.10s  %-10.10s  %-10.10s  %-10.10s') 
  
print('Gaussian input file is written in the file name %s' %name_out) 
print('Happy Calculation') 