import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt
import csv

def R_Calculate(group_name, set_list, dose_list): 
    output_directory = f"./results/figures/single_dose_impact_ratio_analysis/{group_name}/"
    os.makedirs(output_directory, exist_ok=True) 
    R_list = []
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
        average_sum = row_sums.mean()
        R_list.append(row_sums)
    print(f"Data calculated")

    return R_list

def Dose_Index_Calculate(group_name, R_C_w, R_A_w, R_C_n, R_A_n):
    output_directory = f"./results/figures/single_dose_impact_ratio_analysis/{group_name}_dose_index/"
    os.makedirs(output_directory, exist_ok=True)    

    # calculate R_C_w / R_A_w and R_C_n / R_A_n
    ratio_w = [[R_C_w[i][j] / R_A_w[i][j] for j in range(len(R_C_w[i]))] for i in range(len(R_C_w))]

    ratio_n = [[R_C_n[i][j] / R_A_n[i][j] for j in range(len(R_C_n[i]))] for i in range(len(R_C_n))]

    last_20_w = ratio_w[-20:]
    last_20_n = ratio_n[-20:]
    
    mean_w = [np.mean(sorted(sublist, reverse=True)[:5]) for sublist in last_20_w]
    var_w = [np.var(sorted(sublist, reverse=True)[:5]) for sublist in last_20_w]
    
    mean_n = [np.mean(sublist) for sublist in last_20_n]
    var_n = [np.var(sublist) for sublist in last_20_n]

    dose_values = [round(0.1 + i * 0.1, 1) for i in range(20)]

    plt.figure(figsize=(10, 6))
    plt.errorbar(dose_values, mean_w, yerr=np.sqrt(var_w), label='PTEN-WT', marker='o', capsize=5, linestyle='-', color='#2075B6')
    plt.errorbar(dose_values, mean_n, yerr=np.sqrt(var_n), label='PTEN-Null', marker='s', capsize=5, linestyle='--', color='#FF7F0C')
    plt.title(f'Dose Index for {group_name}')
    plt.xlabel('Dose')
    plt.ylabel('Ratio')
    plt.xticks(dose_values)
    plt.legend()
    output_path = os.path.join(output_directory, f"{group_name}.pdf")
    plt.savefig(output_path, format='pdf')
    plt.close()
    print(f"Plot saved to {output_path}")    

    
def main():
    group_name = "CI"

    group_name1 = "p_ten_synergy_wild_C"
    set_list = list(range(21))
    dose_list = list(np.round(np.arange(0.0, 2.1, 0.1), 1))

    R_C_w = R_Calculate(group_name1, set_list, dose_list)

    group_name2 = "p_ten_synergy_wild_I"
    set_list = list(range(21))
    dose_list = list(np.round(np.arange(0.0, 2.1, 0.1), 1))
    
    R_A_w = R_Calculate(group_name2, set_list, dose_list)

    group_name1 = "p_ten_synergy_null_C200"
    set_list = list(range(20))
    dose_list = list(np.round(np.arange(0.1, 2.1, 0.1), 1))

    R_C_n = R_Calculate(group_name1, set_list, dose_list)

    group_name2 = "p_ten_synergy_null_I200"
    set_list = list(range(20))
    dose_list = list(np.round(np.arange(0.1, 2.1, 0.1), 1))
    
    R_A_n = R_Calculate(group_name2, set_list, dose_list)
    
    Dose_Index_Calculate(group_name, R_C_w, R_A_w, R_C_n, R_A_n)
    
if __name__ == "__main__":
    main()
