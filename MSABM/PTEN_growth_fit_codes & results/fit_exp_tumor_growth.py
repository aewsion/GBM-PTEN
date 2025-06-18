import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt

def process_csv_files(group_name, case_list, n_min_list, n_max_list):

    combined_data_dict = {m: [] for m in case_list}
    length = len(case_list)

    for m in range(length):
        output_directory = f"./figures/{group_name}/"
        os.makedirs(output_directory, exist_ok=True) 

        for n in range(n_min_list[m], n_max_list[m]):
            # Construct the file path
            file_path = f"./{group_name}/set_{n}/tc_num.csv"
            
            # Check if the file exists
            if os.path.exists(file_path):
                # Read the CSV file as a DataFrame
                try:
                    # Attempt to read the CSV file as a DataFrame
                    data = pd.read_csv(file_path, header=None) 
                    
                except pd.errors.ParserError:
                    print("error! skip first line")
                    data = pd.read_csv(file_path, header=None, skiprows=1)
            
                data = data.iloc[-1:, :]

                # Check if the number of rows is less than 48 * 100
                if data.shape[1] < 48 * 100:
                    # Get the last element from the last row
                    last_element = data.iloc[-1, -1] 
                    
                    # Pad to 48 * 100 rows
                    missing_count = 48 * 100 - data.shape[1]
                    
                    # Create new rows based on last_element
                    if last_element < 100:
                        new_rows = pd.DataFrame(0, index=[0], columns=np.arange(data.shape[1], data.shape[1] + missing_count))  # Create new rows with zeros (1 row, multiple columns)
                    else:
                        new_rows = pd.DataFrame(last_element, index=[0], columns=np.arange(data.shape[1], data.shape[1] + missing_count))  # Create new rows with last_element (1 row, multiple columns)

                    # Concatenate the new rows to the data
                    data = pd.concat([data, new_rows], axis=1, ignore_index=True)

                # Append the DataFrame to the combined data dictionary
                combined_data_dict[case_list[m]].append(data)  # Append the DataFrame data to the combined list
            else:
                print(f"Warning: File not found - {file_path}")

        # Combine all rows into a DataFrame
        combined_df = pd.DataFrame(np.vstack(combined_data_dict[case_list[m]]))

        # Determine the normalization denominator
        denominator = 2154

        # Normalize the DataFrame
        combined_df = combined_df / denominator
        combined_df = combined_df.round(3)

        # Save the combined data to the specified path
        output_file_path = os.path.join(output_directory, f"{case_list[m]}_all_individuals_TC.csv")
        combined_df.to_csv(output_file_path, index=False, header=False)

        print(f"Data of individuals successfully saved to: {output_file_path}")
        
def read_csv_files(group_name, case_list): 
    combined_data_list = []
    for m in case_list:
        read_directory = f"./figures/{group_name}/"
        read_file_path = os.path.join(read_directory, f"{m}_all_individuals_TC.csv")
        combined_df = pd.read_csv(read_file_path, header=None)
        print("read_csv.shape:", combined_df.shape)
        combined_data_list.append(combined_df)
    
        # plot individuals
        plt.figure(figsize=(10, 6))  # You can adjust the figure size as needed
    
        dead_threshold = 0.85
        alive = combined_df[(combined_df.iloc[:, -4800:-1] <= dead_threshold).all(axis=1)]
        dead = combined_df[~((combined_df.iloc[:, -4800:-1] <= dead_threshold).all(axis=1))]

        for index in range(alive.shape[0]):
            plt.plot(alive.iloc[index, :-1], color='black', linewidth=0.5, label='alive individuals' if index == 0 else "")
        
        for index in range(dead.shape[0]):
            plt.plot(dead.iloc[index, :-1], color='#648CC8', linewidth=0.5, label='dead individuals' if index == 0 else "")


        plt.legend()

        # Set x and y limits
        plt.xlim(0, 48 * 100)  # Set x-axis limits
        plt.ylim(0, 1)    # Set y-axis limits
        xticks = np.linspace(0, 48 * 100, num=11)  
        xtick_labels = np.linspace(0, 100, num=11).astype(int) 
        plt.xticks(xticks, xtick_labels) 

        # Set labels and title
        plt.xlabel('time (days)')  # Change the label as needed
        plt.ylabel('cell density')
        plt.title(f'{m}_all_individuals_TC')
        plt.legend(loc='upper right', fontsize='small') 
        plt.text(500, 0.1, f"Number of alive individuals: {alive.shape[0]}", fontsize=10, color='black')
        plt.text(500, 0.05, f"Number of dead individuals: {dead.shape[0]}", fontsize=10, color='#648CC8')

        # Define the output path
        output_file_path = os.path.join(read_directory, f"{m}_multiple_predictions.pdf")
        
        # Save the figure as a PDF
        plt.savefig(output_file_path, format='pdf', bbox_inches='tight')
        plt.close() 
        print(f"Figure of individuals successfully saved to: {output_file_path}")        
        
    return combined_data_list
def fit_exp_growth(group_name): 

    read_directory = f"./figures/{group_name}/"
    read_file_path = os.path.join(read_directory, f"{group_name}_all_individuals_TC.csv")
    combined_df = pd.read_csv(read_file_path, header=None)
    print("read_csv.shape:", combined_df.shape)

    # Experimental PTEN-WT tumor volume data
    x_WT_exp = [10 * 48, 15 * 48, 23 * 48] 
    y_WT_exp = [0.04, 0.04 * 1.5, 0.04 * 6] 

    # Experimental PTEN-Null tumor volume data
    x_null_exp = [8 * 48, 22 * 48] 
    y_null_exp = [0.04, 0.04 * 19.6] 
    err_null_exp = [0.0269, 0.1704] 

    # plot individuals
    plt.figure(figsize=(10, 6)) 
    if 'null' in group_name:
        for index in range(combined_df.shape[0]):
            plt.plot(combined_df.iloc[index, :48*30], color='#4682B4', linewidth=0.1, label='pten_null individuals' if index == 0 else "")
        plt.errorbar(
            x_null_exp, y_null_exp,
            yerr=err_null_exp,
            fmt='s',              
            color='#4682B4',
            markersize=8,
            capsize=4,             
            elinewidth=1.5,      
            label='pten_null experiment points'
        )
    if 'wild' in group_name:
        for index in range(combined_df.shape[0]):
            plt.plot(combined_df.iloc[index, :48*30], color='#B22222', linewidth=0.1, label='pten_wild individuals' if index == 0 else "")
        plt.scatter(
            x_WT_exp, y_WT_exp,
            color='#B22222',     
            marker='^',          
            s=60,                
            edgecolor='#B22222',   
            label='pten_wild experiment points'   
        )      
    plt.legend()

    # Set x and y limits
    plt.xlim(0, 48*30) 
    plt.ylim(0, 1)   
    xticks = np.linspace(0, 48*30, num=4)  
    xtick_labels = np.linspace(0, 30, num=4).astype(int) 
    plt.xticks(xticks, xtick_labels) 
    # Set labels and title
    plt.xlabel('time (days)') 

    plt.ylabel('cell density')
    plt.title(f'{group_name}_all_individuals_TC')
    plt.legend(loc='upper right', fontsize='small')

    # Define the output path
    output_file_path = os.path.join(read_directory, f"{group_name}_simulation_fit_exp.pdf")
    
    # Save the figure as a PDF
    plt.savefig(output_file_path, format='pdf', bbox_inches='tight')
    plt.close() 
    print(f"Figure of individuals successfully saved to: {output_file_path}")    
    

        
def main():

    case_list = ["p_ten_growth_wild"] # or "p_ten_growth_null"
    n_min_list = [0]
    n_max_list = [50]
    group_name = "p_ten_growth_wild"

    process_csv_files(group_name, case_list, n_min_list, n_max_list)
    fit_exp_growth(group_name)
    
if __name__ == "__main__":
    main()
