#Some .csv files might have not been converged. This script shows these .csv files 

import os

def find_missing_results(folder_path):
    subset_files = []
    result_files = []

    # Get list of subset_*.csv files
    for file in os.listdir(folder_path):
        if file.startswith('subset_') and file.endswith('.csv'):
            subset_files.append(file)

    # Get list of result_subset_*.csv files
    results_folder = os.path.join(folder_path, 'results')
    if os.path.exists(results_folder) and os.path.isdir(results_folder):
        for file in os.listdir(results_folder):
            if file.startswith('result_subset_') and file.endswith('.csv'):
                result_files.append(file)

    # Find missing result_subset_*.csv files
    missing_results = [file for file in subset_files if file.replace('subset_', 'result_subset_') not in result_files]
    return missing_results


def main():
    folders = ['folder1', 'folder2', 'folder3', 'folder4', 'folder5']
    script_dir = os.path.dirname(os.path.realpath(__file__))

    for folder in folders:
        folder_path = os.path.join(script_dir, folder)  
        missing_results = find_missing_results(folder_path)
        print(f'number of missing in {folder} is {len(missing_results)}')
 

if __name__ == "__main__":
    main()
