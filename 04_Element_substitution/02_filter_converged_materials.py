#!/usr/bin/env python
# coding: utf-8

#get a total number of converged materials, remove compositions which did not converge
import pandas as pd

def process_csv_files(csv_files_list_path):
    with open(csv_files_list_path, 'r') as file:
        csv_files = file.read().splitlines()

    total_sum = 0

    for file_path in csv_files:
        df = pd.read_csv(file_path)
        df.dropna(inplace=True)
        df.to_csv(file_path, index=False)
        total_sum += len(df)

    return total_sum

def main():
    csv_files_list_path = 'csv_files_list.txt'
    total_sum = process_csv_files(csv_files_list_path)
    print("Total sum of lengths of all .csv files:", total_sum)

if __name__ == "__main__":
    main()


