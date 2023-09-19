# SAXS_score
Python script to fit and score theoretical ensemble-averaged SAXS profile against experimental SAXS data


## General description

For a set of conformational models and their corresponding population fractions (weights), this script at first computes theoretical SAXS profiles for an ensemble, or loads previously computed 
SAXS profiles (see repository [SAXS_profile](https://github.com/mpopara/SAXS_profile)), if these are available under pdbPath. Subsequently, the ensemble-averaged SAXS profile is computed as weighted linear 
combination of SAXS profiles of ensemble members, by using population fractions of ensemble members as weights.
 
In the next step, theoretical ensemble-averaged profile is fit to experimentally recorded scattering profile, by optimizing a set of fit parameters (c1, c2, offset) in the following user-defined range:

```
    fit_settings = {
        'min_c1': 0.95, 
        'max_c1': 1.05,
        'min_c2': -2.0, 
        'max_c2': 4.0,
        'use_offset': True,
    }
```

till minimum descrepancy to experimental data, measured by &chi;<sup>2</sup>, is obtained.<sup>1</sup>

Offset accounts for systematic errors in the experimental data due to mismatched buffers, while parameters c1 and c2 are part of the atomic form factor function. Parameter c1 has a role in adjusting the total
excluded volume of the atoms, while c2 has a role in adjusting the difference between the densities of the hydration layer around the protein and the bulk water.<sup>1</sup> 

Computation of SAXS profiles as well as their fitting against experimental data is performed using _saxs_ module<sup>1</sup> of the Integrative Modelling Platform<sup>2</sup> (IMP). 
In principle, saxs module also allows to score theoretical SAXS profile against experimental data without parameter fitting.

Results of the fitting are saved as _saxs.fit_ file, which in the header contains fit parameters as well as the &chi;<sup>2</sup>. This is then followed with four columns containing range of q values, 
experimental scattering intensity, experimental error and fitted theoretical scattering profile.

## Example data and input file requirements

In the folder [example_data](https://github.com/mpopara/SAXS_score/tree/main/example_data), provided are exemplary input files:

* /PDBs/*.pdb- ensemble of structural models, provided as individual .pdb files.
* conformer_weights.dat file containing weights (population fractions) of ensemble members. This space-delimited file is of a size N<sub>conformers</sub> x 2, where the first column contains indices of the ensemble members,
 and the second column contains their corresponding weights. This script assumes that the numbering of esemble members in their file name follows the same order as in the weights file. 
* experimental_SAXS_profile.dat- file containing experimentally recorded scattering profile against which the theoretical SAXS profile computed for an ensemble with be scored.
This is a three column file, containing range of q values, scattering intensity and experimental error.

## Dependencies

saxs_scoring.py is a python script bult on Python 3.8.8. Script was tested with provided examplary input files under the following configuration:

* Windows 10
* Python 3.8.8
* IMP 2.17.0
* numpy 1.23.0


## References

1. Schneidman-Duhovny, D.; Hammel, M.; Tainer, John A.; Sali, A., Accurate SAXS Profile Computation and its Assessment by Contrast Variation Experiments. Biophys. J. 2013, 105 (4), 962-974.

2. Russel, D.; Lasker, K.; Webb, B.; Velázquez-Muriel, J.; Tjioe, E.; Schneidman-Duhovny, D.; Peterson, B.; Sali, A., Putting the Pieces Together: Integrative Modeling Platform Software for Structure Determination of Macromolecular Assemblies. PLoS Biol. 2012, 10 (1), e1001244.



## Authors

* Thomas-Otavio Peulen
* Milana Popara
