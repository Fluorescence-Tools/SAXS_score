# -*- coding: utf-8 -*-
"""
Created on Thu Jul 28 11:35:05 2022

@author: popara
"""

import sys
sys.path.insert(1,'C:\\Program Files\\IMP-2.17.0\\python\\') 
import IMP
import IMP.atom
import IMP.core
import IMP.saxs

import numpy as np
import pathlib
import os
import typing
import tqdm

def get_model_profile(
    pdb_fn: pathlib.Path, 
    model_delta_q = 0.5 / 500,
    model_min_q = 0.0,
    model_max_q = 0.5
) -> IMP.saxs.Profile:
    m = IMP.Model()
    pdb_fn = str(pdb_fn)
    saxs_fn = pdb_fn + '.saxs.dat'
    # calculate SAXS profile
    model_profile = IMP.saxs.Profile(model_min_q, model_max_q, model_delta_q)
    if os.path.exists(saxs_fn):
        model_profile.read_SAXS_file(saxs_fn)
    else:
        mp = IMP.atom.read_pdb(pdb_fn, m, IMP.atom.NonWaterNonHydrogenPDBSelector(), True, True)
        # select particles from the model
        particles = IMP.atom.get_by_type(mp, IMP.atom.ATOM_TYPE)
        # add radius for water layer computation
        ft = IMP.saxs.get_default_form_factor_table()
        for i in range(0, len(particles)):
            radius = ft.get_radius(particles[i])
            IMP.core.XYZR.setup_particle(particles[i], radius)
        # compute surface accessibility
        s = IMP.saxs.SolventAccessibleSurface()
        surface_area = s.get_solvent_accessibility(IMP.core.XYZRs(particles))
        model_profile.calculate_profile_partial(particles, surface_area)
        # Write SAXS curve
        model_profile.write_SAXS_file(saxs_fn)
    return model_profile


def score_weighted_profiles(
    model_weights: typing.List[float],
    fit_file_name: str
) -> float:
    # Compute weighted average over all models
    model_profile = IMP.saxs.Profile(model_min_q, model_max_q, model_delta_q)
    for p, w in zip(model_profiles, model_weights):
        model_profile.add(p, w)
    fit_settings = {
        'min_c1': 0.95, # c1 adjusts the excluded volume
        'max_c1': 1.05,
        'min_c2': -2.0, # c2 adjusts the density of hydration layer
        'max_c2': 4.0,
        'use_offset': True,
        'fit_file_name': fit_file_name
    }
    # calculate chi-square score 
    saxs_score = IMP.saxs.ProfileFitterChi(exp_profile)
    fit_parameters = saxs_score.fit_profile(model_profile, *fit_settings.values())
    chi2 = fit_parameters.get_chi_square()
    return chi2




########## input data path and input parameters


basePath = pathlib.Path('C:/user/folder')
expPath = basePath/'experimental_SAXS_profile.dat'
pdbPath = basePath /'PDBs/'
weightsPath = basePath/'conformer_weights.dat'



model_delta_q = 0.6 / 500 
model_min_q = 0.0
model_max_q = 0.6  

################################################################

# read experimental profile: IMP.saxs.Profile(file_name, fit_file, max_q, units)
#  file_name: profile file name
#  fit_file: if true, intensities are read from column 3
#  max_q: read till maximal q value = max_q, or all if max_q<=0
#  units: gets 1, 2, or 3 for unknown q units, 1/A, or 1/nm


exp_profile = IMP.saxs.Profile(str(expPath), False, 0.55, 3)
print('min_q = ' + str(exp_profile.get_min_q()))
print('max_q = ' + str(exp_profile.get_max_q()))
print('delta_q = ' + str(exp_profile.get_delta_q()))

# Compute model profiles for PDBs or if saxs profiles already exist in pdbPath, it will read those
pdb_fns = sorted(list(pdbPath.glob('*.pdb')))

model_profiles = list()
for fn in tqdm.tqdm(pdb_fns):
    model_profiles.append(
        get_model_profile(fn, model_delta_q, model_min_q, model_max_q)
    )




weights = np.loadtxt(weightsPath)[:,1]
fit_file_name = str(basePath / 'saxs.fit') # file containing fit parameters, as well as the exp and fitted theoretical saxs profile  
chi2 = score_weighted_profiles(weights, fit_file_name)

print('chi2 = ' + str(chi2))
    
       





