import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt

def process_csv_files(group_name, case_list, n_min_list, n_max_list, file_name):
    # Initialize a dictionary to store combined data for each m
    combined_data_dict = {m: [] for m in case_list}
    length = len(case_list)
   
    for m in range(length):
        output_directory = f"./figures/{group_name}/"
        os.makedirs(output_directory, exist_ok=True) 

        for n in range(n_min_list[m], n_max_list[m]):       
            # Construct the file path
            file_path = f"./{group_name}/set_{n}/{file_name}.csv"
            
            # Check if the file exists
            if os.path.exists(file_path):
                # Read the CSV file as a DataFrame
                try:
                    # Attempt to read the CSV file as a DataFrame
                    data = pd.read_csv(file_path, header=None) 
                    
                except pd.errors.ParserError:
                    print("error! skip first line")
                    data = pd.read_csv(file_path, header=None, skiprows=1)
                
                data = data.iloc[-1:, :48 * 30]

                # Check if the number of rows is less than 48 * 30
                if data.shape[1] < 48 * 30:
                    # Get the last element from the last row
                    last_element = data.iloc[-1, -1]  # Get the last element from the last row
                    
                    # Pad to 48 * 30 rows
                    missing_count = 48 * 30 - data.shape[1]
                    
                    # Create new rows based on last_element
                    if last_element > 0:
                        new_rows = pd.DataFrame(last_element, index=[0], columns=np.arange(data.shape[1], data.shape[1] + missing_count))  # Create new rows with last_element (1 row, multiple columns)
                        data = pd.concat([data, new_rows], axis=1, ignore_index=True)

                # Append the DataFrame to the combined data dictionary
                if data.shape[1] == 48 * 30:
                    combined_data_dict[case_list[m]].append(data)  # Append the DataFrame data to the combined list
            else:
                print(f"Warning: File not found - {file_path}")

        # Combine all rows into a DataFrame
        combined_df = pd.DataFrame(np.vstack(combined_data_dict[case_list[m]]))
        print("save_csv.shape: ", combined_df.shape)
        # Determine the normalization denominator
        if 'tc_num' in file_name:
            denominator = 2154
            print(f"{file_name}_max: ", denominator)
        elif 'm1' in file_name or 'm2' in file_name or 'm_num' in file_name:
            denominator = 450
            print(f"{file_name}_max: ", denominator)
        else:
            denominator = 1
            print(f"{file_name}_max: ", denominator)
        # Normalize the DataFrame
        combined_df = combined_df / denominator
        combined_df = combined_df.round(3)
        # Save the combined data to the specified path
        output_file_path = os.path.join(output_directory, f"{case_list[m]}_all_individuals_{file_name}.csv")
        combined_df.to_csv(output_file_path, index=False, header=False)

        print(f"Data of individuals successfully saved to: {output_file_path}")
        
def read_csv_files(group_name, case_list, file_name): 

    color_list = ["black", "#6495ED", "#5F9EA0", "#F4A460", "#3CB371"]
    if 'null' in group_name:
        day_max = 30 * 48
    elif 'wild' in group_name:
        day_max = 30 * 48
        
    k = 0
    plt.figure(figsize=(10, 6))

    for m in case_list:
        read_directory = f"./figures/{group_name}/"
        read_file_path = os.path.join(read_directory, f"{m}_all_individuals_{file_name}.csv")
        combined_df = pd.read_csv(read_file_path, header=None)

        data = combined_df.iloc[:, :-1].to_numpy()
        
        mean_values = np.mean(data, axis=0)
        std_values = np.std(data, axis=0)
        
        plt.plot(mean_values, color=color_list[k], linewidth=1, label=f'{m}')

        indices = np.linspace(0, data.shape[1] - 1, data.shape[1] // 10, dtype=int)
        x_sampled = indices
        mean_sampled = mean_values[indices]
        std_sampled = std_values[indices]
    
        plt.fill_between(
            x_sampled,
            mean_sampled - std_sampled,
            mean_sampled + std_sampled,
            color=color_list[k],
            alpha=0.2
        )
        k = k + 1

    # Set x and y limits
    plt.xlim(0, 48 * 30)  # Set x-axis limits
    if 'Gal9_en' in file_name:
        plt.ylim(0, 1.3)
    else:
        plt.ylim(0, 1)    # Set y-axis limits
    xticks = np.linspace(0, 48 * 30, num=4)  
    xtick_labels = np.linspace(0, 30, num=4).astype(int) 
    plt.xticks(xticks, xtick_labels) 
    # Set labels and title
    plt.xlabel('time (days)')  # Change the label as needed
    if 'tc' in file_name or 'm1' in file_name or 'm2' in file_name or 'm_num' in file_name:
        plt.ylabel('cell density')
    elif '_I' in file_name:
        plt.ylabel('drug concentration')
    else:
        plt.ylabel('cytokine concentration')
    plt.title(f'{group_name}_{file_name}')
    plt.legend(fontsize='small')  # Optionally show a legend  
    plt.xlim([0, day_max])
    # Define the output path
    output_file_path = os.path.join(read_directory, f"{group_name}_{file_name}.pdf")
    
    # Save the figure as a PDF
    plt.savefig(output_file_path, format='pdf', bbox_inches='tight')
    plt.close()  # Close the figure to free up memory
    print(f"Figure of individuals successfully saved to: {output_file_path}")        


def main():
    case_list = ["growth_fit"]
    n_min_list = [0]
    n_max_list = [50]

    group_name_list = ["p_ten_growth_wild", "p_ten_growth_null"]
    
    file_list = ['avgCSF1', 'avgIGF1', 'avgEGF', 'avgGal9_en',
                'avgERK', 'avgAKT', 'avgIGF1R', 'avgEGFR', 
                'avgGSK3b', 'avgIRF1', 'avgGal9','m_num']
         
    for file in file_list:
        process_csv_files(group_name_list[0], case_list, n_min_list, n_max_list, file)
        read_csv_files(group_name_list[0],  case_list, file)
    
        process_csv_files(group_name_list[1], case_list, n_min_list, n_max_list, file)
        read_csv_files(group_name_list[1], case_list, file)
    
if __name__ == "__main__":
    main()
