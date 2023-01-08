# VASPTools
Python Script to read, edit, modify and visualize and file conversion for VASP and  Gaussian input and output! 

This is actively developing version. Current version is very first trial. 

Install: 
    Download the repo and install using pip - 
    pip install -e VASPTools/

Execution: 
* running from Shell command 
eg. python -m vasputils -job job_name --options options 

-job : task to perform  
   vasp2gaus : convert poscar or Contcar file to xyz or gaussian input file format 
   gaus2poscar: Convert gaussian input file to VASP POSCAR file format. Periodic cell (Tv) should be given in the .com file. 


   options: 
     --file : input file name 
     --fout : output file name 
     --ftype: output file format xyz or gaus (gaussian input file) 
     --select : Boolean, select coordinates in POSCAR file
     --direct: direct or cartesian coordinates in POSCAR file format 