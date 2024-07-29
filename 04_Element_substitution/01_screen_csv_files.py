#!/usr/bin/env python
# coding: utf-8

#get a list of all .csv files in the folders
import os

def get_csv_files(folder_path):
    csv_files = [os.path.join(folder_path, file) for file in os.listdir(folder_path) if file.endswith('.csv')]
    return csv_files

def main():
    folders = ['folder1/results', 'folder2/results', 'folder3/results', 'folder4/results', 'folder5/results']
    #folders = ['test_results']
    all_csv_files = []

    for folder in folders:
        csv_files = get_csv_files(folder)
        all_csv_files.extend(csv_files)

    with open('csv_files_list.txt', 'w') as file:
        file.write('\n'.join(all_csv_files))

if __name__ == "__main__":
    main()

