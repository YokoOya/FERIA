These codes are a computer program for calculating a three-dimentional model of an infalling-rotating envelope and a Keplerian disk. 
Fits files of a cube data and a position-velocity diagram are created.
This program uses [CFITSIO library](https://heasarc.gsfc.nasa.gov/fitsio/).


The formalism adopted in this program and how to use it are summarized in [Oya et al. (2022)](https://ui.adsabs.harvard.edu/abs/2022PASP..134i4301O/abstract). 


# Contents
- Makefile
- feria.h
- feria_PV.cpp
- feria_env.cpp
- feria_fitsio.cpp
- feria_main.cpp
- feria_mesh.cpp
- feria_parameter.cpp
- feria_sky.cpp
- feria_varParameters.py
- readme.txt
- template.in


## Files to be Edited 
You need to edit the following files.
- Makefile
- feria.h
- template.in
- (Optional) feria_varParameters.py


### Makefile
- Please edit the linker for cfitsio according to your environment.

### feria.h
- Header file.
- Please edit the values for 'logNpix' and 'logNvel' to change the number of mesh.
- The number of mesh for the position axes in the cube file is 2^{logNpix}, and that for the velocity axes is 2^{logNvel}.

### template.in
- Template for the input file.
- This file includes the parameters for the model.
- See comments in the file for the details of the parameters.

### (Optional) feria_varParameters.py
- This python script makes fits files for various parameter sets.
- Please edit the lists for the phsical parameters.
- This script overwrite the input file and execute the 'feria' file.


## Modifying the Model
You can edit the following files to modify the model calculation.
- feria_env.cpp
- feria_sky.cpp

### feria_env.cpp
- This file contains the class for the physical conditions of the envelope gas.
- By default, the density profile of the emitter and the gas temperature profile are the power law of the distance from the protostar.
- You can edit the functions "feriaModel::env_dens" and "feriaModel::env_temp" to change these profiles.

### feria_sky.cpp
- This file contains the functions to calculate the observed intensity.
- By default, the emissivity is assumed to be proportional to the column density of the emitter.
- Neither of the excitation effect or the radiation transfer is considered.
You can edit the function "skyPlane::projection" to change how to calculate the emissivity.


# Usage

```
make
./feria < $(your inputfile name)
```


## References 
CFITSIO:  Pence, W. (1999) in ASP Conf. Ser., Vol. 172, Astronomical Data Analysis Software and Systems VIII, ed. D. Mehringer, R. Plante, and D. Roberts (San Francisco: ASP), 487
FERIA: [Oya, Y., Kibukawa, H., Miyake, S., and Yamamoto, S. (2022)](https://ui.adsabs.harvard.edu/abs/2022PASP..134i4301O/abstract)

### Contact 
- GitHub: [https://github.com/YokoOya/feria](https://github.com/YokoOya/feria)
- e-Mail: Yoko Oya (yoko.oya_at_yukawa.kyoto-u.ac.jp) <!-- (oya_at_taurus.phys.s.u-tokyo.ac.jp) -->

