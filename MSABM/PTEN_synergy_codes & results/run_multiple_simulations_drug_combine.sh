MAX_CPU_USAGE=85 # maximum CPU usage
CHECK_INTERVAL=60 # check interval = 60s
SIMU_N=1

# Precompute the index mapping based on dose_1 and dose_2 sequences
DOSE_1_ARRAY=($(seq 0.1 0.1 1.0)) # Array of dose_1 values
DOSE_2_ARRAY=($(seq 0.1 0.1 1.0)) # Array of dose_2 values

# Calculate total number of combinations of dose_1 and dose_2
total_combinations=$((${#DOSE_1_ARRAY[@]} * ${#DOSE_2_ARRAY[@]}))

for set in $(seq 0 $((SIMU_N - 1))) # for each individual
do  
  for i in $(seq 0 $((${#DOSE_1_ARRAY[@]} - 1))) # loop through index of dose_1
  do
    for j in $(seq 0 $((${#DOSE_2_ARRAY[@]} - 1))) # loop through index of dose_2
    do
      dose_1=${DOSE_1_ARRAY[$i]} # Get the actual dose_1 value
      dose_2=${DOSE_2_ARRAY[$j]} # Get the actual dose_2 value
      
      # Calculate the unique index based on dose_1 and dose_2
      index=$((set * total_combinations + i * ${#DOSE_2_ARRAY[@]} + j))
      
      while true
      do
          cpu_usage=$(top -bn1 | grep "Cpu(s)" | awk '{print 100 - $8}') # current CPU usage
          current_jobs=$(ps aux | grep -E "p_ten_synergy" | grep -v grep | wc -l) # current number of 'p_ten_synergy' jobs
      
          echo "current CPU usage: ${cpu_usage}% current number of jobs: $current_jobs"
          if (( $(echo "$cpu_usage < $MAX_CPU_USAGE" | bc -l) )); then
              break  
          fi

          sleep $CHECK_INTERVAL
      done     
      echo "current CPU usage: ${cpu_usage}% current number of jobs: $current_jobs"
      
      # start jobs: use two types of drug for combination
      nohup ./bin/p_ten_synergy p_ten_synergy_null_IA200 $index 5 200 0 0 0 $dose_1 $dose_2 &
    done
  done
done

  
  # ./bin/p_ten_synergy \       # path to the executable file for continuous CSF1R_I treatment case
  # p_ten_synergy_wild_A \      # output folder name to store the population results (argv[1])
  # $index \                    # output set name to store each individual results (argv[2])
  # 5 \                         # number of simulations for an individual (argv[3])
  # 200 \                       # simulation time 200 days (argv[4])
  # 0 \                         # 0 for non-responders and 1 for responders (argv[5])
  # 1 \                         # pten = 1 for p_ten_wild; pten = 0 for p_ten_null
  # 0.0                         # CSF1R_I dose
  # $dose_1                     # IGF1R_I dose
  # $dose_2                     # AKT_I dose


