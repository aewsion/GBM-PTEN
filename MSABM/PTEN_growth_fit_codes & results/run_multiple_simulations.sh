MAX_CPU_USAGE=95 # maximum CPU usage
CHECK_INTERVAL=60 # check interval = 60s

for set in {0..49} # for each individual
do  
  while true
  do
      cpu_usage=$(top -bn1 | grep "Cpu(s)" | awk '{print 100 - $8}') # current CPU usage
      current_jobs=$(ps aux | grep -E "p_ten_growth" | grep -v grep | wc -l) # current number of 'p_ten_growth' jobs
  
      echo "current CPU usage: ${cpu_usage}% current number of jobs: $current_jobs"
      if (( $(echo "$cpu_usage < $MAX_CPU_USAGE" | bc -l) )); then
          break  
      fi

      sleep $CHECK_INTERVAL
  done     
  echo "current CPU usage: ${cpu_usage}% current number of jobs: $current_jobs"
  
  # start p_ten_growth jobs 
  nohup ./bin/p_ten_growth p_ten_growth_wild $set 1 100 0 1 &
  nohup ./bin/p_ten_growth p_ten_growth_null $set 1 50 0 0 &
  # ./bin/p_ten_growth \       # path to the executable file for p_ten_growth case
  # p_ten_growth_wild \        # output folder name to store the population results (argv[1])
  # $set \                     # output set name to store each individual results (argv[2])
  # 1 \                        # number of simulations for an individual (argv[3])
  # 100 \                      # simulation time 100 or 50 days (argv[4])
  # 1                          # whether to enable randomness for the key parameters, 1 for enable and 0 for disable (argv[5])
done

