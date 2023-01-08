# VaspUtils
Python Script to read, edit, modify and visualize and file conversion for VASP and  Gaussian input and output! 

This is actively developing version. Current version is very first trial. 

### Installation  
* With pip `pip install -e VASPTools/` 

## Running VaspUtils 
```python 
python  python -m vasputils [-job] [vasp2gaus/gaus2poscar] [--file input file] [--fout output file] [--ftype xyz/gaus (gaussian input file)] [--select boolean] [--direct boolean]  
```
* `-job` option to select task to execute, currently available task 
    **  `vasp2gaus` : Convert poscar or Contcar file to xyz or gaussian input file format 
   ** `gaus2poscar` : Convert gaussian input file to VASP POSCAR file format. Periodic cell (Tv) should be given in the .com file. 
* `--file` : The input file name 
* `--fout` : The output file name 
* `--ftype` :  Output file format, `xyz` or `gaus` (gaussian input file) 
* `--select` : `Boolean`, select coordinates in POSCAR file
* `--direct` : `Boolean`, direct or cartesian coordinates in POSCAR file format 




---
#### Disclaimer:
#### It is an incomplete package, still under development.
---
#### License: 

VaspUtils is freely available under an [MIT](https://opensource.org/licenses/MIT) License
