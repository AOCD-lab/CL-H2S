#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generates a database of sulfides with estimation of Gibbs free energies
In the end gererates a csv-file ""
@author: tolstop
"""

import pandas as pd
from pymatgen.ext.matproj import MPRester
from pymatgen.entries.computed_entries import GibbsComputedStructureEntry, Composition

def has_excluded_element(formula):
    composition = Composition(formula)
    elements = [str(element) for element in composition.elements]
    return any(element in exclude_elements for element in elements)

def count_atoms(formula):
    composition = Composition(formula)
    total_atoms = sum(composition.values())
    return total_atoms



numbers_elements = [2, 3]
names = ['one_metal_sulfides', 'two_metal_sulfides']
for number_elements, name in zip(numbers_elements, names):
    with MPRester(api_key="your_ip_key") as mpr:
        docs = mpr.summary.search(elements='S',
                                  num_elements=number_elements, 
                                  fields=["material_id", 
                                          "formula_pretty",
                                          'composition_reduced',
                                          "formation_energy_per_atom",
                                          "energy_above_hull",
                                          "structure",
                                          "theoretical",
                                          "symmetry"])
    
       
    '''
    Create a dictionary of this informations excluding 'structure' 
    AND opposite value for 'theoretical' for changing it to "experimentally_observed" column
    create a DataFrame
    '''    
    mpid_dict = {doc.material_id: [doc.formula_pretty, 
                                   doc.symmetry.crystal_system,
                                   doc.symmetry.symbol,
                                   not doc.theoretical, 
                                   doc.energy_above_hull,
                                   doc.formation_energy_per_atom] for doc in docs}
    columns=["pretty_formula",
             "crystal_system",
             "spacegroup",
             "experimentally_observed",
             "energy_above_hull, eV/atom",
             "Ef"]   
    df =  pd.DataFrame.from_dict(mpid_dict,orient='index', columns=columns)
    
    exclude_elements = ['H','O','C','N','Se',
                     'He', 'Ne','Ar','Kr','Xe','Rn','P',
                     'O', 'F', 'Cl', 'I', 'Br', 'I',
                     'Ac','Th','Pa','U','Np','Pu']


    # Filter materials with less than 20 atoms in the formula
    df['total_atoms'] = df['pretty_formula'].apply(count_atoms)
    if name == 'one_metal_sulfides':
        df = df[df['total_atoms'] < 20]
    elif name == 'two_metal_sulfides':
        df = df[df['total_atoms'] < 30]


    #Filter materials which have non-desired elements inside
    df = df[~df['pretty_formula'].apply(has_excluded_element)]
    df = df[(df['energy_above_hull, eV/atom'] < 0.1)]
    #Among polymorphs, leave only those with lowest Ef
    df_sorted = df.sort_values(by='Ef')
    idx_min_ef = df_sorted.groupby('pretty_formula')['Ef'].idxmin()
    df = df_sorted.loc[idx_min_ef]

    print(f'len of dataframe if {len(df)}')
    '''
    Gf calculation for all the materials
    Some of entries give KeyError, we will fill their predicted Gibbs energy with NaN, after delete rows containing NaN 
    Creating a dictionary of {MP-ID: Gf} at certain temperature, adding this column to dataframe
    '''
    temp_list = [i for i in range(400, 2050, 100)]
    
    for temp in temp_list:  
        Gibbs_formation_energy_dict = {}
        for doc in docs:   
            try:
                Gibbs = GibbsComputedStructureEntry(doc.structure, 
                                                    formation_enthalpy_per_atom=doc.formation_energy_per_atom, 
                                                    temp=temp)
            except KeyError:
                Gibbs_value = None
            else:
                comp = Composition(doc.composition_reduced)
                Gibbs_value = Gibbs.gf_sisso() / Gibbs.structure.num_sites * comp.num_atoms
            Gibbs_formation_energy_dict[doc.material_id]= Gibbs_value
        df[f'Gf_{temp}'] = Gibbs_formation_energy_dict
    
    df.dropna(inplace=True)
    
    df.to_csv(f'{name}_MP_database.csv')
    
    
