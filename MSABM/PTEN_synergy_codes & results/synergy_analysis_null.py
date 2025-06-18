import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import csv


def R12_Calculate(group_name, set_list, dose_C_list, dose_A_list, dose_C_index, dose_A_index): 
    output_directory = f"./results/figures/synergy_analysis_null/{group_name}/"
    os.makedirs(output_directory, exist_ok=True)  # Create directory if it doesn't exist
    R12_data_frame = pd.DataFrame(np.zeros((len(dose_C_index), len(dose_A_index))), index=dose_C_index, columns=dose_A_index)
    for n in set_list:
        # Construct the file path
        file_path = f"./results/{group_name}/set_{n}/tc_num.csv"
        
        # Check if the file exists
        if os.path.exists(file_path):
            with open(file_path, "r") as f:
                reader = csv.reader(f)
                data_list = list(reader)
        
            max_cols = 9600
            data_list = [row + [row[-1]] * (max_cols - len(row)) for row in data_list]
            data = pd.DataFrame(data_list, dtype=float)
            row_sums = data.sum(axis=1)
        
        # Calculate the average of the row sums
        average_sum = row_sums.mean()
        R12_data_frame.at[dose_C_list[n], dose_A_list[n]] = np.round(average_sum/2154, 3)

    output_file_path = os.path.join(output_directory, "R12.csv")
    R12_data_frame.to_csv(output_file_path)

    print(f"Data saved to {output_file_path}")
    return R12_data_frame

def R_Calculate(group_name, set_list, dose_list): 
    output_directory = f"./results/figures/synergy_analysis_null/{group_name}/"
    os.makedirs(output_directory, exist_ok=True)  # Create directory if it doesn't exist
    R_data_frame = pd.DataFrame(np.zeros(len(dose_list)), index=dose_list)
    for n in set_list:
        file_path = f"./results/{group_name}/set_{n}/tc_num.csv"
        if os.path.exists(file_path):
            with open(file_path, "r") as f:
                reader = csv.reader(f)
                data_list = list(reader)
            max_cols = 9600
            data_list = [row + [row[-1]] * (max_cols - len(row)) for row in data_list]
            data = pd.DataFrame(data_list, dtype=float)

        row_sums = data.sum(axis=1)
        # Calculate the average of the row sums
        average_sum = row_sums.mean()
        R_data_frame.at[dose_list[n]] = np.round(average_sum/2154, 3)
        
    output_file_path = os.path.join(output_directory, "R.csv")
    R_data_frame.to_csv(output_file_path)

    print(f"Data saved to {output_file_path}")

    return R_data_frame

def Synergy_Index_Calculate(group_name, R12, R1, R2, dose1_name, dose2_name):
    output_directory = f"./results/figures/synergy_analysis_null/{group_name}_index/"
    os.makedirs(output_directory, exist_ok=True)  # Create directory if it doesn't exist    
    
    R = pd.DataFrame(np.minimum(R1, R2), index=R1.index, columns=R1.columns)

    SI = pd.DataFrame(np.zeros(R12.shape), index=R12.index, columns=R12.columns)
    
    for i in SI.index:
        for j in SI.columns:
            ij_rounded = round(i+j, 1)
            if ij_rounded in R.index:
                SI.loc[i, j] = R.loc[ij_rounded].iloc[0] / R12.loc[i, j]
            else:
                print(f"Index {ij_rounded} not found in R")
                
    output_file_path = os.path.join(output_directory, "SI.csv")
    SI.to_csv(output_file_path)
    print(f"Data saved to {output_file_path}")
    
    plt.figure(figsize=(10, 8))  
    sns.heatmap(SI, cmap='viridis', annot=True, fmt='.2f', linewidths=0.5)
    plt.gca().invert_yaxis()
    plt.title(f"Synergy Index {group_name}")
    plt.xlabel(f'{dose2_name} dose')
    plt.ylabel(f'{dose1_name} dose')
    output_figure_path = os.path.join(output_directory, f"SI_{group_name}.pdf")
    plt.savefig(output_figure_path, format='pdf')
    
def main():
    group_name = "p_ten_synergy_null_CA200"
    set_list = list(range(100))
    dose_C_values = np.round(np.arange(0.1, 1.1, 0.1), 1)
    dose_A_values = np.round(np.arange(0.1, 1.1, 0.1), 1) 

    dose_C_index = list(dose_C_values)
    dose_A_index = list(dose_A_values)
    
    dose_C_values = np.round(np.arange(0.1, 1.1, 0.1), 1)
    dose_A_values = np.round(np.arange(0.1, 1.1, 0.1), 1)    
    dose_A_list, dose_C_list = np.meshgrid(dose_A_values, dose_C_values)
    dose_C_list = dose_C_list.flatten()
    dose_A_list = dose_A_list.flatten()
    R_CA = R12_Calculate(group_name, set_list, dose_C_list, dose_A_list, dose_C_index, dose_A_index)

    group_name1 = "p_ten_synergy_null_C200"
    set_list = list(range(20))
    dose_list = list(np.round(np.arange(0.1, 2.1, 0.1), 1))

    R_C = R_Calculate(group_name1, set_list, dose_list)

    group_name2 = "p_ten_synergy_null_A200"
    set_list = list(range(20))
    dose_list = list(np.round(np.arange(0.1, 2.1, 0.1), 1))
    
    R_A = R_Calculate(group_name2, set_list, dose_list)

    
    Synergy_Index_Calculate(group_name, R_CA, R_C, R_A, 'CSF1R_I', 'AKT_I')
    
if __name__ == "__main__":
    main()
