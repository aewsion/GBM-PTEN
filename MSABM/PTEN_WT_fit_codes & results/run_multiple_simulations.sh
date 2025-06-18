MAX_CPU_USAGE=95 # maximum CPU usage
CHECK_INTERVAL=60 # check interval = 60s
SIMU_N=100
for set in $(seq 0 $((SIMU_N - 1))) # for each individual
do  
  while true
  do
      cpu_usage=$(top -bn1 | grep "Cpu(s)" | awk '{print 100 - $8}') # current CPU usage
      current_jobs=$(ps aux | grep -E "p_ten" | grep -v grep | wc -l) # current number of 'p_ten' jobs
  
      echo "current CPU usage: ${cpu_usage}% current number of jobs: $current_jobs"
      if (( $(echo "$cpu_usage < $MAX_CPU_USAGE" | bc -l) )); then
          break  
      fi

      sleep $CHECK_INTERVAL
  done     
  echo "current CPU usage: ${cpu_usage}% current number of jobs: $current_jobs"
  
  # start p_ten jobs 
  nohup ./bin/p_ten p_ten_wild_drug $((set + SIMU_N*0)) 1 200 1 1 1.0 &
  nohup ./bin/p_ten p_ten_wild_nodrug $((set + SIMU_N*0)) 1 100 1 1 0.0 &

  # ./bin/p_ten \          # path to the executable file for p_ten jobs under continuous CSF1R_I treatment
  # p_ten_wild_drug \      # output folder name to store the population results (argv[1])
  # $set \                 # output set name to store each individual results (argv[2])
  # 1 \                    # number of simulations for an individual (argv[3])
  # 200 \                  # simulation time, 200 days for CSF1R_I-treatment-case and 100 days for no-treatment-case (argv[4])
  # 1 \                    # whether to enable randomness for the key parameters, 1 for enable and 0 for disable (argv[5])
  # 1 \                    # 0 for PTEN-Null, 1 for PTEN-WT (argv[6])
  # 1.0                    # drug dose of CSF1R_I (argv[7])
done

