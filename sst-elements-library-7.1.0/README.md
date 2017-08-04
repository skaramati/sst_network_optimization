--------------------------------------------
File list
--------------------------------------------
path: sst-elements-library-7.1.0/src/sst/elements/ember/run/configurations/network_optimization

topology_optimizer.py		Topology optimization
main_opt_LB.py			Latency and bandwidth optimizations


--------------------------------------------
Command Line Arguments-topology_optimizer.py
--------------------------------------------

Example:  ./topology_optimizer.py --motif "FFT3D iteration=1 nx=64 ny=64 nz=64" --output_file=out.csv

--motif (required) - Input motif and its arguments 
--output_file (required) - output file name to save the simulation result 


These were the steps to the process:
1. Implement the motif
2. Add the following line to your motif:
	printf("###execution_time %f\n", execution time);
3. Run topology optimizer




--------------------------------------------
Command Line Arguments-main_opt_LB.py
--------------------------------------------
Example: ./main_opt_LB.py  --output_file=out.csv --bandwidth=200 --latency=1e-9 --iterations=5 --tolerance=0.1


--output_file (required) - output file name to save the simulation result
--bandwidth (required) - upper limit on bandwidth value
--latency (required) - lower limit on latency value
--iterations (required) - number of iterations to run the algorithm with increased tolerance value in each iteration
--tolerance (required) - set the tolerance limit for execution time

These were the steps to the process:
1. Implement the motif
2. Add the following line to your motif:
	printf("###execution_time %f\n", execution time)
3. In main_opt_LB.py, update def run_sst(w) based on the motif and system parameters
4. Run main_opt_LB.py



#### Copyright (c) 2009-2017, Sandia Corporation
Portions are copyright of other developers:
See the file CONTRIBUTITORS.TXT in the top level directory
of this repository for more information.

##### [LICENSE](https://github.com/sstsimulator/sst-elements/blob/devel/LICENSE)
