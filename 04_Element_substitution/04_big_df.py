#!/usr/bin/env python
# coding: utf-8

#makes one database from all .csv files

import pandas as pd

def stack_csv_files(csv_files_list_path, output_file_path):
    with open(csv_files_list_path, 'r') as file:
        csv_files = file.read().splitlines()

    dfs = []
    for file_path in csv_files:
        df = pd.read_csv(file_path)
#        df.drop(['structure', 'spacegroup', 'parent_structure', 'Ef_MEGNET', 'Unnamed: 0.1', 'Unnamed: 0','parent_formula'], axis=1, inplace=True)
        dfs.append(df)

    combined_df = pd.concat(dfs, ignore_index=True)
    combined_df.to_csv(output_file_path, index=False)

def main():
    csv_files_list_path = 'csv_files_list.txt'
    output_file_path = 'not_filtered_new_materials_Hd.csv'
    stack_csv_files(csv_files_list_path, output_file_path)
    print("Combined CSV file saved as", output_file_path)

if __name__ == "__main__":
    main()

