# ReadGV
Python Script to read, edit, modify and visualize and file conversion for VASP and  Gaussian input and output! 

This is actively developing version. Current version is very first trial. 

### Installation  
* pip install -e git+git://github.com/MAJabed/ReadGV.git 
* With pip `pip install readgv`


## Running ReadGV 
```
python -m readgv [-job] [vasp2gaus/gaus2poscar/read_bands**] [--file input file] [--fout output file] [--ftype xyz/gaus (gaussian input file)] [--select boolean] [--direct boolean]  
```
* `-job` option to select task to execute, currently available task 
    * `vasp2gaus` or `vasp2com` : Convert poscar or Contcar file to xyz or gaussian input file format 
    * `gaus2poscar` or `gaus2vasp` : Convert gaussian input file to VASP POSCAR file format. Periodic cell (Tv) should be given in the .com file. 
    * `read_bands` : Read selected or all bands from OUTCAR and write in a file, print or plot (follow `--plot`)
    * `get_geom` : Write or print molecular dynamic trajectory from the OUTCAR files 
    * `abs_gaus` : Read Gaussian TDDFT output file, Dress excited energies using Gaussian fucntion to plot absorption spectra 

* `--file` : The input file name 
* `--fout` : The output file name 
* `--ftype` :  Output file format, `xyz` or `gaus` (gaussian input file) 
* `--select` : `Boolean`, select coordinates in POSCAR file
* `--direct` : `Boolean`, direct or cartesian coordinates in POSCAR file format 
* `--bands` : No of bands to read, default 100 
* `--frames` : MD frames to read,  last frame  `final`, initial frame `initial`, all `all`, last n frames `-n`, first n frame `n`
* `--nthkpoint` : integar, read outcar data of nth Kpoints
* `--coord` : output coordinates in Cartesian `cart` or fractional `frac` 
* `--plot` : Plot output of `-job`, special arguments: `dos_all` to plot DOS of all frames in waterfall graphs (This function may not produce the best figure in come cases) 
* `--shiftfermi` : Band energies shift by Fermi energy 
* `--xlim` : X-axis limit, input type: list - `lower upper` 
* `--sigma` : float, linewidth for Gaussian dressing, default 0.1 eV
* `--theta` : Rotation of X-axis for 3D view of waterfall graph, degree unit. 
* `--gamma` : Rotation of Z axis for 3D view of waterfall graph, degree unit.
* `--alpha` : Adjust the transparency of a matplotlib graph plot`[0-1]`  
* `--unit` : Unit of a calculation, eg. `eV` or `nm` in absorption spectra calculations 
* `--figsize`: Matplotlib graph figure size, default `(7,5)`

---
#### Disclaimer:
#### It is an incomplete package, still under development.
---
#### License: 

ReadGV is freely available under an [MIT](https://opensource.org/licenses/MIT) License
