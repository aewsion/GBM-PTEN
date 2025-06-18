import pandas as pd
import os
import numpy as np
import lifelines
import matplotlib.pyplot as plt

def process_csv_files(group_name, case_list, n_min_list, n_max_list):

    combined_data_dict = {m: [] for m in case_list}
    length = len(case_list)

    for m in range(length):
        output_directory = f"./results/figures/{group_name}/"
        os.makedirs(output_directory, exist_ok=True)  # Create directory if it doesn't exist

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

                # Check if the number of rows is less than 9600
                if data.shape[1] < 9600:
                    # Get the last element from the last row
                    last_element = data.iloc[-1, -1] 

                    # Pad to 9600 rows
                    missing_count = 9600 - data.shape[1]
                    
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
        
def read_csv_files(group_name, case_list, threshold): 
    combined_data_list = []
    color_list = ["#6495ED", "#4169E1", "#90EE90", "#2E8B57"]
    if 'p_ten_wild_nodrug' in group_name:
        day_max = 100 * 48
    elif 'p_ten_wild_drug' in group_name:
        day_max = 200 * 48
        
    k = 1    
    plt.figure(figsize=(10, 6))

    for m in case_list:
        read_directory = f"./results/figures/{group_name}/"
        read_file_path = os.path.join(read_directory, f"{m}_all_individuals_TC.csv")
        combined_df = pd.read_csv(read_file_path, header=None)
        print("read_csv.shape:", combined_df.shape)
        if len(combined_data_list) == 0:
            combined_data_list = combined_df  
        else:
            combined_data_list = pd.concat([combined_data_list, combined_df], axis=0) 

    dead_threshold = threshold
    alive = combined_data_list[(combined_data_list.iloc[:, -4800:-1] <= dead_threshold).all(axis=1)]
    dead = combined_data_list[~((combined_data_list.iloc[:, -4800:-1] <= dead_threshold).all(axis=1))]

    for index in range(alive.shape[0]):
        plt.plot(alive.iloc[index, :-1], color=color_list[k], linewidth=0.5, label=f'alive individuals' if index == 0 else "")
    
    for index in range(dead.shape[0]):
        plt.plot(dead.iloc[index, :-1], color=color_list[k+2], linewidth=0.5, label=f'dead individuals' if index == 0 else "")


    alive_num = alive.shape[0]
    dead_num = dead.shape[0]
        
    plt.legend()
    plt.xlim(0, 9600)  # Set x-axis limits
    plt.ylim(0, 1)    # Set y-axis limits
    xticks = np.linspace(0, 9600, num=11)  
    xtick_labels = np.linspace(0, 200, num=11).astype(int) 
    plt.xticks(xticks, xtick_labels) 
    # Set labels and title
    plt.xlabel('time (days)')  # Change the label as needed
    plt.ylabel('cell density')
    plt.title(f'{m}_all_individuals_TC')
    plt.legend(loc='upper right', fontsize='small')  # Optionally show a legend  
    plt.xlim([0, day_max])
    
    plt.text(500, 0.8, f"Number of alive individuals: {alive_num}", fontsize=10, color=color_list[k])
    plt.text(500, 0.85, f"Number of dead individuals: {dead_num}", fontsize=10, color=color_list[k+2])

    # Define the output path
    output_file_path = os.path.join(read_directory, f"{group_name}_multiple_predictions.pdf")
    
    # Save the figure as a PDF
    plt.savefig(output_file_path, format='pdf', bbox_inches='tight')
    plt.close()
    print(f"Figure of individuals successfully saved to: {output_file_path}")        
        
    return combined_data_list

    
def logranktest_exp(group_ABM_name, group_ABM_data, threshold_list, day_list):
    exp_data_read_directory = f"./exp_data_p_ten_wild/exp_{group_ABM_name}.csv"
    save_directory = f"./results/figures/{group_ABM_name}/"
    os.makedirs(save_directory, exist_ok=True)
    # Number of individuals
    ABM_data = group_ABM_data
    n_ABM = ABM_data.shape[0]  # Number of rows (individuals)
    threshold = threshold_list
    day = day_list
    if 'drug' in group_ABM_name:
        day_max = 200
    elif 'nodrug' in  group_ABM_name:
        day_max = 100

    # prepare data for group ABM prediction
    ABM_time = []
    ABM_event = []
    
    for i in range(n_ABM):
        # Determine death status
        death_condition = ABM_data.iloc[i, :] > threshold
        death_condition[:day*48] = False
        if ABM_data.iloc[i, -1] > 0:
            if any(death_condition):
                death_time_index = day*48-1 if day*48-1 >= death_condition.idxmax() else death_condition.idxmax()  # Get the first index where death occurs
                ABM_time.append((death_time_index+1) * 0.5 / 24)  # Calculate death time
                ABM_event.append(1)  # Event occurred (death)
            else:
                ABM_time.append(day_max)  # Censoring at the end of follow-up
                ABM_event.append(0)  # Event did not occur (censored)
        else:
            ABM_time.append(day_max)  # Censoring at the end of follow-up
            ABM_event.append(0)  # Event did not occur (censored)  
                
    exp_data = pd.read_csv(f'{exp_data_read_directory}')
    exp_data.iloc[:, 0] = exp_data.iloc[:, 0].astype(int)  
    exp_data.iloc[:, 2] = exp_data.iloc[:, 2].astype(int)  
    n_exp = exp_data.shape[0]
    
    # Prepare data for group basal
    exp_time = exp_data.iloc[:n_exp, 1].tolist()
    exp_event = exp_data.iloc[:n_exp, 2].tolist()

    # Create DataFrame for Log-Rank test
    survival_data = pd.DataFrame({
        'time': ABM_time + exp_time,
        'event': ABM_event + exp_event,
        'group': [f'{group_ABM_name}'] * n_ABM + [f'experiment'] * n_exp
    })
    # Perform Log-Rank test
    group_ABM = survival_data[survival_data['group'] == f'{group_ABM_name}']
    group_exp = survival_data[survival_data['group'] == f'experiment']

    kmf_ABM = lifelines.KaplanMeierFitter()
    kmf_ABM.fit(group_ABM['time'], event_observed=group_ABM['event'], label=f'{group_ABM_name}')
    kmf_exp = lifelines.KaplanMeierFitter()
    kmf_exp.fit(group_exp['time'], event_observed=group_exp['event'], label=f'experiment')
    results = lifelines.statistics.logrank_test(group_ABM['time'], group_exp['time'], event_observed_A=group_ABM['event'], event_observed_B=group_exp['event'])

    fig = plt.figure(figsize=(14, 12))
    kmf_ABM.plot_survival_function(ci_show=False)
    kmf_exp.plot_survival_function(ci_show=False)
    plt.annotate(f'{group_ABM_name} vs experiment log-rank-test p-value: {results.p_value:.4g}', xy=(0.95, 0.05), xycoords='axes fraction', fontsize=12, color='blue', horizontalalignment='right')
    plt.title('survival curves')
    plt.xlabel('time (days)')
    plt.ylabel('percent survival')
    plt.xlim([0, day_max])
    current_ticks = plt.yticks()[0]
    plt.yticks(current_ticks, [int(tick * 100) for tick in current_ticks])
    plt.ylim([0, 1])
    plt.legend()
    save_path = f'{save_directory}/{group_ABM_name}_exp_KM.pdf'
    plt.savefig(save_path, format='pdf')
    plt.close() 

def logranktest(group_ABM_name_nodrug, group_ABM_data_nodrug, group_ABM_name_drug, group_ABM_data_drug, threshold, day_nodrug, day_drug):

    save_directory = f"./results/figures/{group_ABM_name_nodrug}_vs_{group_ABM_name_drug}/"
    os.makedirs(save_directory, exist_ok=True)
    # Number of individuals
    ABM_data_nodrug = group_ABM_data_nodrug
    n_ABM_nodrug = ABM_data_nodrug.shape[0]  # Number of rows (individuals)
    day_max = 100

    # prepare data for group ABM prediction
    ABM_time_nodrug = []
    ABM_event_nodrug = []
    
    for i in range(n_ABM_nodrug):
        # Determine death status
        death_condition = ABM_data_nodrug.iloc[i, :] > threshold
        death_condition[:day_nodrug*48] = False
        if ABM_data_nodrug.iloc[i, -1] > 0:
            if any(death_condition):
                death_time_index = day_nodrug*48-1 if day_nodrug*48-1 >= death_condition.idxmax() else death_condition.idxmax()  # Get the first index where death occurs
                ABM_time_nodrug.append((death_time_index+1) * 0.5 / 24)  # Calculate death time
                ABM_event_nodrug.append(1)  # Event occurred (death)
            else:
                ABM_time_nodrug.append(day_max)  # Censoring at the end of follow-up
                ABM_event_nodrug.append(0)  # Event did not occur (censored)
        else:
            ABM_time_nodrug.append(day_max)  # Censoring at the end of follow-up
            ABM_event_nodrug.append(0)  # Event did not occur (censored)  
                
    # Number of individuals
    ABM_data_drug = group_ABM_data_drug
    n_ABM_drug = ABM_data_drug.shape[0]  # Number of rows (individuals)

    day_max = 200

    # prepare data for group ABM prediction
    ABM_time_drug = []
    ABM_event_drug = []
    
    for i in range(n_ABM_drug):
        # Determine death status
        death_condition = ABM_data_drug.iloc[i, :] > threshold
        death_condition[:day_drug*48] = False
        if ABM_data_drug.iloc[i, -1] > 0:
            if any(death_condition):
                death_time_index = day_drug*48-1 if day_drug*48-1 >= death_condition.idxmax() else death_condition.idxmax()  # Get the first index where death occurs
                ABM_time_drug.append((death_time_index+1) * 0.5 / 24)  # Calculate death time
                ABM_event_drug.append(1)  # Event occurred (death)
            else:
                ABM_time_drug.append(day_max)  # Censoring at the end of follow-up
                ABM_event_drug.append(0)  # Event did not occur (censored)
        else:
            ABM_time_drug.append(day_max)  # Censoring at the end of follow-up
            ABM_event_drug.append(0)  # Event did not occur (censored)  


    # Create DataFrame for Log-Rank test
    survival_data = pd.DataFrame({
        'time': ABM_time_nodrug + ABM_time_drug,
        'event': ABM_event_nodrug + ABM_event_drug,
        'group': [f'{group_ABM_name_nodrug}'] * n_ABM_nodrug + [f'{group_ABM_name_drug}'] * n_ABM_drug
    })
    # Perform Log-Rank test
    group_ABM_nodrug = survival_data[survival_data['group'] == f'{group_ABM_name_nodrug}']
    group_exp_drug = survival_data[survival_data['group'] == f'{group_ABM_name_drug}']

    kmf_ABM_nodrug = lifelines.KaplanMeierFitter()
    kmf_ABM_nodrug.fit(group_ABM_nodrug['time'], event_observed=group_ABM_nodrug['event'], label=f'{group_ABM_name_nodrug}')
    kmf_ABM_drug = lifelines.KaplanMeierFitter()
    kmf_ABM_drug.fit(group_exp_drug['time'], event_observed=group_exp_drug['event'], label=f'{group_ABM_name_drug}')
    results = lifelines.statistics.logrank_test(group_ABM_nodrug['time'], group_exp_drug['time'], event_observed_A=group_ABM_nodrug['event'], event_observed_B=group_exp_drug['event'])
    print("p-value", results.p_value)  # Print Log-Rank test results

    fig = plt.figure(figsize=(14, 12))
    kmf_ABM_nodrug.plot_survival_function(ci_show=False)
    kmf_ABM_drug.plot_survival_function(ci_show=False)
    plt.annotate(f'{group_ABM_name_nodrug} vs {group_ABM_name_drug} log-rank-test p-value: {results.p_value:.4g}', xy=(0.95, 0.05), xycoords='axes fraction', fontsize=12, color='blue', horizontalalignment='right')
    plt.title('survival curves')
    plt.xlabel('time (days)')
    plt.ylabel('percent survival')
    plt.xlim([0, day_max])
    current_ticks = plt.yticks()[0]
    plt.yticks(current_ticks, [int(tick * 100) for tick in current_ticks])
    plt.ylim([0, 1])
    plt.legend()
    save_path = f'{save_directory}/{group_ABM_name_nodrug}_vs_{group_ABM_name_drug}_KM.pdf'
    plt.savefig(save_path, format='pdf')
    plt.close() 

def main():

    case_list = ["p_ten_wild_nonresponders", "p_ten_wild_responders"]
    n_min_list = [0, 100]
    n_max_list = [0, 100]
    threshold = 0.732

    group_name_list = ["p_ten_wild_nodrug", "p_ten_wild_drug"]
    day_list = [0, 28]
        
    process_csv_files(group_name_list[0], case_list, n_min_list, n_max_list)

    data_0 = read_csv_files(group_name_list[0], case_list, threshold)
    logranktest_exp(group_name_list[0], data_0, threshold, day_list[0])

    process_csv_files(group_name_list[1], case_list, n_min_list, n_max_list)
    data_1 = read_csv_files(group_name_list[1], case_list, threshold)
    logranktest_exp(group_name_list[1], data_1, threshold, day_list[1])

    logranktest(group_name_list[0], data_0, group_name_list[1], data_1, threshold, day_list[0], day_list[1])
if __name__ == "__main__":
    main()
