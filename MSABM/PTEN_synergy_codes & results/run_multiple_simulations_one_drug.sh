MAX_CPU_USAGE=90 # maximum CPU usage
CHECK_INTERVAL=60 # check interval = 60s
SIMU_N=1

DOSE_ARRAY=($(seq 0.0 0.1 2.0)) # Array of dose values

for set in $(seq 0 $((SIMU_N - 1))) # for each individual
do  
  for i in $(seq 0 $((${#DOSE_ARRAY[@]} - 1))) # loop through index of dose_I
  do
      dose=${DOSE_ARRAY[$i]} # Get the actual dose_I value

      
      # Calculate the unique index based on dose_I and dose_A
      index=$i 

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
      
      # start jobs: use only one type of drug
      nohup ./bin/p_ten_synergy p_ten_synergy_wild_A $index 5 200 0 1 0 0 $dose &

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
  # 0.0                         # IGF1R_I dose
  # $dose                       # AKT_I dose


