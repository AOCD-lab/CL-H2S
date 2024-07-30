#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
from pymatgen.entries.computed_entries import Composition
import re
from itertools import permutations
from math import gcd


def one_metal_coefficients(row):
    '''
    Parameters
    ----------
    row : pandas.core.frame.DataFrame
        
    Returns
    -------
    a: int
    b: int
    Returns indexes (a, b) in formula (MaSb)

    '''
    name = row['pretty_formula'].values[0]
    composition_dict = dict(Composition(name).to_reduced_dict)
    for element in composition_dict.keys():
        if element == 'S':
            b = int(composition_dict[element])
        else: 
            a = int(composition_dict[element])
    return a, b
    
def two_metal_coefficients(row):
    '''
    Parameters
    ----------
    row : pandas.core.frame.DataFrame
        
    Returns
    -------
    alpha: int
    beta: int
    gamma: int
    eps: int
    Returns indexes (alpha, beta, gamma, eps) in formula (M1(eps*alpha)M2(eps*beta)Sgamma)

    '''
    composition_reduced = Composition(row['pretty_formula'].values[0])
    Sn = ''.join(re.findall(r'S\d+', str(composition_reduced)))
    composition_without_Sn = ' '.join(list(filter(None, map(lambda x: x.strip(), (str(composition_reduced).split(Sn))))))

    formula, eps = Composition(composition_without_Sn).get_reduced_composition_and_factor()
    composition_dict = dict(formula.to_reduced_dict)
    alpha, beta = composition_dict.values()
    alpha, beta = int(alpha), int(beta)
    

    composition_Sn = dict(Composition(Sn).as_dict())
    gamma = composition_Sn['S']
    gamma = int(gamma)

    return alpha, beta, gamma, eps



def get_coefficients_one_metal_reaction(material_1_row, material_2_row):
    a, b = one_metal_coefficients(material_1_row)
    c, d = one_metal_coefficients(material_2_row)

    material_1_coeff = c / (a * d - c * b)
    material_2_coeff = a / (a * d - c * b)
    
    return [material_1_coeff, material_2_coeff]

def get_coefficients_two_metal_reaction(material_1_row, material_2_row):
    alpha, beta, gamma, i = two_metal_coefficients(material_1_row)
    alpha_2, beta_2, theta, eps = two_metal_coefficients(material_2_row)

    
    material_1_coeff = eps / (i * theta - gamma * eps)
    material_2_coeff = i / (i * theta - gamma * eps)
    
    return [material_1_coeff, material_2_coeff]




def get_material_name(row):
    return row['pretty_formula'].values[0]

def get_name_for_pair(material_1_row, material_2_row):
    return [get_material_name(material_1_row), get_material_name(material_2_row)]

def get_loop(names, coefficients):
    material_1_coeff, material_2_coeff = coefficients
    material_1_name, material_2_name = names
    reaction_I = f'{round(material_1_coeff,2)}{material_1_name} + H2S = {round(material_2_coeff, 2)}{material_2_name} + H2'
    reaction_II = f'{round(material_2_coeff, 2)}{material_2_name} = {round(material_1_coeff,2)}{material_1_name} + (1/2)S(2)'
    return [reaction_I, reaction_II]




def reaction_one_metal(material_1_row, material_2_row):
    '''
    Parameters
    ----------
    material_1_row: pandas.core.frame.DataFrame
    material_2_row: pandas.core.frame.DataFrame


    Returns
    -------
    list
    '''
    coefficients = get_coefficients_one_metal_reaction(material_1_row, material_2_row)
    names = get_name_for_pair(material_1_row, material_2_row)
    reactions = get_loop(names, coefficients)

    return reactions

def reaction_two_metal(material_1_row, material_2_row):
    '''
    Parameters
    ----------
        
    material_1_row: pandas.core.frame.DataFrame
    material_2_row: pandas.core.frame.DataFrame

    Returns
    -------
    list
    '''
    coefficients = get_coefficients_two_metal_reaction(material_1_row, material_2_row)
    names = get_name_for_pair(material_1_row, material_2_row)
    reactions = get_loop(names, coefficients)
    
    return reactions


def get_constants():
    Gf_H2S_list = np.array([-0.387034254,
                   -0.416427424,
                   -0.439436182,
                   -0.45811266,
                   -0.473586568,
                   -0.475410686,
                   -0.42477069,
                   -0.373799036,
                   -0.322630461])

    temp_list = np.array(range(400, 1300, 100))
    
    return Gf_H2S_list, temp_list

def get_Gf_list(row):
    _, temp_list = get_constants()

    return [float(row[f'Gf_{T}']) for T in temp_list]
    
def get_Gf_for_pair(material_1_row, material_2_row):
    return [get_Gf_list(material_1_row), get_Gf_list(material_2_row)]
    
def get_Gr_loop(Gf_pair_lists, coefficients):
    Gf_H2S_list, temp_list = get_constants()
    Gf_1_list, Gf_2_list = Gf_pair_lists
    material_1_coeff, material_2_coeff = coefficients
    
    Gr_I = [material_2_coeff * Gf_2 - Gf_H2S - material_1_coeff * Gf_1 for Gf_1, Gf_2, Gf_H2S in zip(Gf_1_list, Gf_2_list, Gf_H2S_list)]
    Gr_II = [material_1_coeff * Gf_1 - material_2_coeff * Gf_2 for Gf_1, Gf_2 in zip(Gf_1_list, Gf_2_list)]
    return [Gr_I, Gr_II]


def Gr_one_metal(material_1_row, material_2_row):
    '''
    Parameters
    ----------
    Gf_H2S_list: list
    material_1_row: pandas.core.frame.DataFrame
    material_2_row: pandas.core.frame.DataFrame


    Returns
    -------
    list of lists
    Returns Gibbs reaction energies for 2 reaction in sulfur cycle if materials contain only one metal in composition
    '''
    Gf_pair_lists = get_Gf_for_pair(material_1_row, material_2_row)
    coefficients = get_coefficients_one_metal_reaction(material_1_row, material_2_row)
    Gr_reactions = get_Gr_loop(Gf_pair_lists, coefficients)

    return Gr_reactions

def Gr_two_metal(material_1_row, material_2_row):
    '''
    Parameters
    ----------
    Gf_H2S_list: list
    material_1_row: pandas.core.frame.DataFrame
    material_2_row: pandas.core.frame.DataFrame


    Returns
    -------
    list of lists
    Returns Gibbs reaction energies for 2 reaction in sulfur cycle if materials contain two metals in composition
    '''
    
    Gf_pair_lists = get_Gf_for_pair(material_1_row, material_2_row)
    coefficients = get_coefficients_two_metal_reaction(material_1_row, material_2_row)
    Gr_reactions = get_Gr_loop(Gf_pair_lists, coefficients)

    return Gr_reactions
def get_unique_elements(formula):
    composition = Composition(formula)
    return ''.join(sorted(set(str(element) for element in composition.elements)))

def Gf_materials(material_1_row, material_2_row):
    '''
    Parameters
    ----------
    material_1_row: pandas.core.frame.DataFrame
    material_2_row: pandas.core.frame.DataFrame


    Returns
    -------
    list of lists
    Returns Gibbs formation energies for 2 materials in eV/atom
    '''
    _, temp_list = get_constants()
    materials_1_name, material_2_name = get_name_for_pair(material_1_row, material_2_row)
    
    n_atoms_MaSb = Composition(materials_1_name).num_atoms
    n_atoms_McSd = Composition(material_2_name).num_atoms

    G_MaSb_list = [float(material_1_row[f'Gf_{T}'])/n_atoms_MaSb for T in temp_list]
    G_McSd_list = [float(material_2_row[f'Gf_{T}'])/n_atoms_McSd for T in temp_list] 
    
    return [G_MaSb_list, G_McSd_list]

def Hf_0_materials(material_1_row, material_2_row):
    '''
    Parameters
    ----------
    material_1_row: pandas.core.frame.DataFrame
    material_2_row: pandas.core.frame.DataFrame


    Returns
    -------
    list 
    Returns Enthalpy formation energies for 2 materials at 0 K
    '''

    Hf_0_MaSb = float(material_1_row['Ef']) 
    Hf_0_McSd = float(material_2_row['Ef']) 
    
    return [Hf_0_MaSb, Hf_0_McSd]

def MP_id_materials(material_1_row, material_2_row):
    '''
    Parameters
    ----------
    material_1_row: pandas.core.frame.DataFrame
    material_2_row: pandas.core.frame.DataFrame


    Returns
    -------
    list 
    Returns MP-id for 2 materials 
    '''
    MP_id_MaSb = str(material_1_row['MP_id'].iloc[0])
    MP_id_McSd = str(material_2_row['MP_id'].iloc[0])
    
    return [MP_id_MaSb, MP_id_McSd]
      

def get_unique_elements(formula):
    composition = Composition(formula)
    return ''.join(sorted(set(str(element) for element in composition.elements)))

def sulfur_amount(row):
    return row['element_amounts']['S']

def calculate_m1_m2_ratio(row, function_type):    
    row['element_amounts'].pop('S')

    if function_type == 'one_metal':
        m = float(next(iter(row['element_amounts'].values())))
        ratio, gcd_value = 1, m
        return ratio, gcd_value

    elif function_type == 'two_metal':    
        m1, m2 = row['element_amounts'].values()
        ratio = m1 / m2
        gcd_value = gcd(int(m1), int(m2))
        return ratio, gcd_value
    else:
        raise ValueError(f"Invalid function_type: {function_type}")



def process_pair(material_1_row, material_2_row, function_type, exclude_keys):
    reactions = {}
    def calculate_loop(material_1_row, material_2_row, function_type):
        if function_type == 'one_metal':
            return reaction_one_metal(material_1_row, material_2_row)
        elif function_type == 'two_metal':
            return reaction_two_metal(material_1_row, material_2_row)
        else:
            raise ValueError(f"Invalid function_type: {function_type}")

    def calculate_Gr_energy(material_1_row, material_2_row, function_type):
        if function_type == 'one_metal':
            return Gr_one_metal(material_1_row, material_2_row)
        elif function_type == 'two_metal':
            return Gr_two_metal(material_1_row, material_2_row)
        else:
            raise ValueError(f"Invalid function_type: {function_type}")


    def calculate_Gf_materials():
        return Gf_materials(material_1_row, material_2_row)

    def calculate_Hf_0K():
        return Hf_0_materials(material_1_row, material_2_row)

    def calculate_MP_id():
        return MP_id_materials(material_1_row, material_2_row)
    
    if float(material_1_row['Sn_normalized']) < float(material_2_row['Sn_normalized']):
        # Create a dictionary of functions for each key
        function_map = {
            'loop': lambda: calculate_loop(material_1_row, material_2_row, function_type),
            'Gr_energy, eV': lambda: calculate_Gr_energy(material_1_row, material_2_row, function_type),
            'Gf_materials, eV/atom': calculate_Gf_materials,
            'Hf_0K, eV/atom': calculate_Hf_0K,
            'MP_id': calculate_MP_id
        }
        # Exclude specified keys from the reactions dictionary
        if exclude_keys is None:
            include_keys = function_map.keys()
        else:
            include_keys = [key for key in function_map.keys() if key not in exclude_keys]

        # Calculate values for included keys
        for key in include_keys:
            reactions[key] = function_map[key]()

    return reactions


def process_material(material, df_small, reactions_dictionary, function_type, exclude_keys):
    for pair in permutations(df_small['pretty_formula'].tolist(), 2):
        material_1_row = df_small.loc[df_small['pretty_formula'] == pair[0]]
        material_2_row = df_small.loc[df_small['pretty_formula'] == pair[1]]
        reactions = process_pair(material_1_row, material_2_row, function_type, exclude_keys)
        if reactions:
            reactions_dictionary[f'{pair[0]}/{pair[1]}'] = reactions


def process_ratios(material, ratios, df, reactions_dictionary, function_type, exclude_keys):
    for ratio in ratios:
        df_small = df.loc[(df['unique_elements'] == material) & 
                                   (df['M1_M2_ratio'] == ratio)]
        if len(df_small) >= 2:
            process_material(material, df_small, reactions_dictionary, function_type, exclude_keys)
    return reactions_dictionary  # Return the dictionary from the function


def process_compositions(df, function_type='two_metal', exclude_keys=None):
    reactions_dictionary = {}
    
    df['unique_elements'] = df['pretty_formula'].apply(get_unique_elements)
    df['composition'] = df['pretty_formula'].apply(Composition)
    df['element_amounts'] = df['composition'].apply(lambda comp: comp.get_el_amt_dict())
    df['sulfur_amount'] = df.apply(lambda row: sulfur_amount(row), axis=1)
    df[['M1_M2_ratio', 'gcd']] = df.apply(lambda row: calculate_m1_m2_ratio(row, function_type), axis=1, result_type="expand")
    df['Sn_normalized'] = df['sulfur_amount']/df['gcd']
  
    chemspace = set(df['unique_elements'].tolist())

    for subset in chemspace:
        ratios = df.loc[df['unique_elements'] == subset, 'M1_M2_ratio'].unique()
        process_ratios(subset, ratios, df, reactions_dictionary, function_type, exclude_keys)
        
    df.drop(['unique_elements', 'composition', 'element_amounts', 'M1_M2_ratio', 'sulfur_amount', 'gcd'], axis=1, inplace=True)

    return reactions_dictionary



