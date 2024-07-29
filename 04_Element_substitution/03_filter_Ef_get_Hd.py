#!/usr/bin/env python
# coding: utf-8

#The script goes to every .csv file, leaves only compositions with lowest formation enthalpies among polymorphs, and estimates decomposition energy of each new composition

from pymatgen.entries.computed_entries import ComputedEntry
import pandas as pd
from pymatgen.analysis.phase_diagram import PhaseDiagram, PDPlotter
from mp_api.client import MPRester
from pymatgen.entries.computed_entries import Composition




def filter_polymorphs(df):
    df_sorted = df.sort_values(by='Ef_MEGNET_corrected')
    filtered_df = df_sorted.groupby(['M1_M2_ratio', 'Sn_normalized']).first().reset_index()                
    filtered_df.reset_index(inplace=True, drop=True)
    return filtered_df


def get_reference_energy(comp: Composition, el_refs) -> float:
    """Sum of elemental reference energies over all elements in a composition.

    Args:
        comp (Composition): Input composition.

    Returns:
        float: Reference energy
    """
    return sum(comp[el] * el_refs[el].energy_per_atom for el in comp.elements)

                                             additional_criteria={"thermo_types": ["GGA_GGA+U"]}) 
def get_all_entries(df):
    composition = Composition(df['pretty_formula'][0])
    elements = [str(element) for element in composition.elements]

    with MPRester("your_mp_id") as mpr:
        entries = mpr.get_entries_in_chemsys(elements=elements,
                                             additional_criteria={"thermo_types": ["GGA_GGA+U"]}) 

    return entries

def get_MEGNET_entry(row, el_refs):
    composition = Composition(row['pretty_formula'])
    energy = row['Ef_MEGNET_corrected']
    my_entry = ComputedEntry(composition, 
                             energy * composition.num_atoms, 
                             correction=get_reference_energy(composition, el_refs))
    return my_entry
        

def get_Hd_MEGNET(entries, row, el_refs):
    
    new_entry = get_MEGNET_entry(row, el_refs)
    
    entries.append(new_entry)
    
    pd_new = PhaseDiagram(entries)
    
    Hd_MEGNET = pd_new.get_phase_separation_energy(new_entry, stable_only=False)
    return Hd_MEGNET

def get_Hd(row, entries):
    phase_diagram  = PhaseDiagram(entries)
    el_refs = phase_diagram.el_refs

    Hd_MEGNET = get_Hd_MEGNET(entries, row, el_refs)
    
    return Hd_MEGNET



def process_csv_files(csv_files_list_path):
    with open(csv_files_list_path, 'r') as file:
        csv_files = file.read().splitlines()

    for file_path in csv_files:
        df = pd.read_csv(file_path)
        
        df_filter = filter_polymorphs(df)
        entries = get_all_entries(df_filter)

        df_filter['Hd_MEGNET'] = df_filter.apply(lambda row: get_Hd(row, entries), axis=1, result_type="expand")

        df_filter.to_csv(file_path, index=False)  # Resave the file without NaN rows


def main():
    csv_files_list_path = 'csv_files_list.txt'
    process_csv_files(csv_files_list_path)

if __name__ == "__main__":
    main()

