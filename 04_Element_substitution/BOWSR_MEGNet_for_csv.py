#!/usr/bin/env python
# coding: utf-8


from pymatgen.core.structure import Structure, Lattice
from pymatgen.core.periodic_table import get_el_sp
from maml.apps.bowsr.model.megnet import MEGNet
from maml.apps.bowsr.optimizer import BayesianOptimizer
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from timeout_decorator import timeout, TimeoutError
import tensorflow as tf
import pandas as pd
import sys


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

@timeout(120)
def get_relaxed_value(row, model):
    try:
        name = row["pretty_formula"]
        compressed_optimizer = BayesianOptimizer(
            model=model,
            structure=row['structure'],
            relax_coords=True,
            relax_lattice=True,
            use_symmetry=True,
            seed=42
        )

        compressed_optimizer.set_bounds()
        compressed_optimizer.optimize(n_init=100, n_iter=100, alpha=0.026 ** 2)
        cutoff_distance = 1
 
        relaxed, _ = compressed_optimizer.get_optimized_structure_and_energy(cutoff_distance=cutoff_distance)
        energy = model.predict_energy(relaxed)
        return relaxed, energy
    except TimeoutError as te:
        print(f"Timeout error for {name}: {te}")
        return None, None

    except Exception as e:
        print(f'Exception occured for {name}: {e}')
        return None, None

def get_crystal_system_value(row):
    try:
        return SpacegroupAnalyzer(row["relaxed"]).get_crystal_system()
    except Exception as e:
        return None

def get_space_group_symbol_value(row):
    try:
        return SpacegroupAnalyzer(row["relaxed"]).get_space_group_symbol()
    except Exception as e:
        return None
    

def main():
    if len(sys.argv) != 2:
        print("Usage: python script_name.py <csv_file>")
        sys.exit(1)

    csv_file = sys.argv[1]
    df=pd.read_csv(csv_file)
    return df, csv_file


if __name__ == "__main__":
    df, csv_filename = main()
    print(f'running {csv_filename}')    
    df['structure'] = df['structure'].apply(parse_structure_string)
    model = MEGNet()

    df[["relaxed", 'Ef_MEGNET']] = df.apply(lambda row: get_relaxed_value(row, model), axis=1, result_type="expand")
    df["crystal_system"] = df.apply(lambda row: get_crystal_system_value(row), axis=1)
    df["space_group_symbol"] = df.apply(lambda row: get_space_group_symbol_value(row), axis=1)
    
    df.to_csv(f'results/{csv_filename}')
    print(f'excecution for {csv_filename} completed')

