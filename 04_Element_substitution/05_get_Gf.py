#!/usr/bin/env python
# coding: utf-8

#The script estimates Gibbs formation energies for each new composition

import pandas as pd
from pymatgen.entries.computed_entries import GibbsComputedStructureEntry, Composition
from pymatgen.core.structure import Structure, Lattice
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer




def parse_structure_string(structure_string):
    # Split the structure string into lines
    lines = structure_string.split('\n')

    # Extract lattice parameters
    abc = list(map(float, lines[2].split()[2:5]))
    angles = list(map(float, lines[3].split()[1:4]))

    # Create a lattice using the parameters
    lattice = Lattice.from_parameters(a=abc[0], 
                                      b=abc[1], 
                                      c=abc[2], 
                                      alpha=angles[0], 
                                      beta=angles[1], 
                                      gamma=angles[2])

    # Extract atomic coordinates
    elements = []
    coords = []

    for line in lines[8:]:
        if not line.strip():
            break

        parts = line.split()
        element = parts[1]
        x, y, z = map(float, parts[2:5])
        elements.append(element)
        coords.append((x, y, z))

    # Create a Structure object
    structure = Structure(lattice, elements, coords, to_unit_cell=True)
    return structure


def get_Gf(row, temperature):
    try:
        Gibbs = GibbsComputedStructureEntry(row['relaxed'], 
                                            formation_enthalpy_per_atom=row['Ef_MEGNET_corrected'], 
                                            temp=temperature)
    except KeyError:
        Gibbs_value = None
    else:
        comp = Composition(row['pretty_formula'])
        Gibbs_value = Gibbs.gf_sisso() / Gibbs.structure.num_sites * comp.num_atoms
    return Gibbs_value


def process_csv(file_path):
    df = pd.read_csv(file_path)

    df['relaxed'] = df['relaxed'].apply(parse_structure_string)
    temp_list = [i for i in range(400, 2050, 100)]
    for temperature in temp_list:
        df[f'Gf_{temperature}'] = df.apply(lambda row: get_Gf(row, temperature), axis=1)  

    df.dropna(inplace=True)
    new_file_name = 'new_materials_database.csv'
    df.to_csv(new_file_name, index=False)  # Resave the file without NaN rows



def main():
    csv_files_list_path = 'not_filtered_new_materials_Hd.csv'
    process_csv(csv_files_list_path)

if __name__ == "__main__":
    main()



