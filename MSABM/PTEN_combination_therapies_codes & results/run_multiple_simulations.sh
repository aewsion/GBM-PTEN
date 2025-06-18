MAX_CPU_USAGE=95 # maximum CPU usage
CHECK_INTERVAL=60 # check interval = 60s
SIMU_N=5
for set in $(seq 0 $((SIMU_N - 1))) # for each individual
do  
  while true
  do
      cpu_usage=$(top -bn1 | grep "Cpu(s)" | awk '{print 100 - $8}') # current CPU usage
      current_jobs=$(ps aux | grep -E "p_ten_combine" | grep -v grep | wc -l) # current number of 'p_ten_combine' jobs
  
      echo "current CPU usage: ${cpu_usage}% current number of jobs: $current_jobs"
      if (( $(echo "$cpu_usage < $MAX_CPU_USAGE" | bc -l) )); then
          break  
      fi

      sleep $CHECK_INTERVAL
  done     
  echo "current CPU usage: ${cpu_usage}% current number of jobs: $current_jobs"
  
  # start p_ten_combine jobs 
  nohup ./bin/p_ten_combine p_ten_combine_null $((set + SIMU_N*0)) 1 200 0 0 2 0 0 &
  nohup ./bin/p_ten_combine p_ten_combine_null $((set + SIMU_N*1)) 1 200 0 0 1 1 0 &
  nohup ./bin/p_ten_combine p_ten_combine_null $((set + SIMU_N*2)) 1 200 0 0 0 1 1 &
  nohup ./bin/p_ten_combine p_ten_combine_null $((set + SIMU_N*3)) 1 200 0 0 1 0 1 &

  nohup ./bin/p_ten_combine p_ten_combine_wild $((set + SIMU_N*0)) 1 200 0 1 2 0 0 &
  nohup ./bin/p_ten_combine p_ten_combine_wild $((set + SIMU_N*1)) 1 200 0 1 1 1 0 &
  nohup ./bin/p_ten_combine p_ten_combine_wild $((set + SIMU_N*2)) 1 200 0 1 0 1 1 &
  nohup ./bin/p_ten_combine p_ten_combine_wild $((set + SIMU_N*3)) 1 200 0 1 1 0 1 &  
  # ./bin/p_ten_combine \     # path to the executable file for continuous combination treatment case
  # p_ten_combine_null \      # output folder name to store the population results (argv[1])
  # $set \                    # output set name to store each individual results (argv[2])
  # 1 \                       # number of simulations for an individual (argv[3])
  # 200 \                     # simulation time 200 days (argv[4])
  # 0 \                       # 0 for non-responders and 1 for responders (argv[5])
  # 1 \                       # pten = 1 for p_ten_wild; pten = 0 for p_ten_null (argv[6])
  # 1.0                       # CSF1R_I dose (argv[7])
  # 1.0                       # IGF1R_I dose (argv[8])
  # 1.0                       # AKT_I dose (argv[9])
done 

