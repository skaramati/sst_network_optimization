#! /usr/bin/Python
import os, sys
import subprocess
import math
import csv
import numpy as np 
from numpy import linalg as LA
from time import sleep 



N = np.array([[0.0, 0.0], [0.0, 0.0]],dtype=np.float64) #stores 
w = np.array([0.0, 0.0],dtype=np.float64)
w_try = np.array([0.0, 0.0],dtype=np.float64)
w_init = np.array([0.0, 0.0],dtype=np.float64)

#initialization	
npop = 2      # optimization directions
sigma = 1    # noise standard deviation  
topo = "fattree"
shape="16,16:16"





#####################################################################################

def run_sst(w,topo,shape):
	'''the function runs the sst simulation for current set-up
	   @param w: the vector w specifies the latency and bandwidth of network for current simulation [latency, bandwidth]
	   @return exec_time: the fuction return the execution time for current simulation
	'''
	
	option = (
	'\'--cmdLine=Init --cmdLine=\"%s\" --cmdLine=Fini ' 
	'--numNodes=256 '
	'--numRanksPerNode=1 '
	'--topo=%s '
	'--shape=%s '
	'--width=1x1 '
	'--netLAT=%.16fs '
	'--netBW=%.16fGB/s '
	'--rtrArb=merlin.xbar_arb_lru ' 
	'--emberRankmapper=fff\'') % ("FFT3D iterations=1 npRow=16 nx=512 ny=512 nz=512", topo, shape, w[0], w[1])
    	command = 'mpirun -np 32 sst --model-options={0} ../../script/emberLoadCmd.py'.format(option)
	#command = " --model-options=%s ../../script/emberLoadCmd.py" % (option)	
    	print command
    	os.system(command)
	
	with open('output.csv', 'r') as f: #the execution time of sst simulation is stored in output.csv
		exec_time = float(f.readline())
	return exec_time

	
def find_optimal_point(data,w_init,topo,shape,runs,exec_time_init, sim_count, tolerance):
    '''
    golden section search

    to find the minimum of f on [a,b]

    f: a strictly unimodal function on [a,b]



    example:

    >>> f = lambda x: (x-2)**2

    >>> x = gss(f, 1, 5)

    >>> x

    2.000009644875678

    '''
    a = data[0][:2]
    b = data[2][:2]
    m = data[1][:2]
    cost_a = data[0][4]
    cost_b = data[2][4]
    cost_m = data[1][4]
    exec_time_a = data[0][2]
    exec_time_b = data[2][2]
    exec_time_m = data[1][2]
	
   
    
    while True:
	c = (a + m) / 2
    	d = (b + m) / 2
	sigma = 1
        N = [c[0]*2,c[1]]
	w_diff = (N-c) 
	w = c
        curr_error=0
	while (curr_error < tolerance) and (sigma > 0.0001):
		w_try = w + sigma * w_diff
		exec_time = run_sst(w_try,topo,shape)
		sim_count +=1
		curr_error = (exec_time - exec_time_init) / exec_time_init
		while (sigma > 0.0001) and (curr_error >= tolerance):
			sigma = sigma / 3.0
			w_try = w + sigma * w_diff
			exec_time = run_sst(w_try,topo,shape)
			sim_count +=1
			curr_error = (exec_time - exec_time_init) / exec_time_init
		w = w_try
	cost_c = cost_model(w,w_init)		
  	exec_time_c = exec_time
	
	
	sigma = 1
        N = [d[0]*2,d[1]]
	w_diff = (N-d) 
	w = d
        curr_error=0
	while (curr_error < tolerance) and (sigma > 0.0001):
		w_try = w + sigma * w_diff
		exec_time = run_sst(w_try,topo,shape)
		sim_count +=1
		curr_error = (exec_time - exec_time_init) / exec_time_init
		while (sigma > 0.0001) and (curr_error >= tolerance):
			sigma = sigma / 3.0
			w_try = w + sigma * w_diff
			exec_time = run_sst(w_try,topo,shape)
			sim_count +=1
			curr_error = (exec_time - exec_time_init) / exec_time_init
		w = w_try
	cost_d = cost_model(w,w_init)		
	exec_time_d = exec_time
	
	
	
        if (cost_c > cost_m) and (cost_d > cost_m):
            a = c
	    b = d
	    exec_time_a = exec_time_c 
            exec_time_b = exec_time_d 
	    cost_a = cost_c 
            cost_b = cost_d 
        else:
            if (cost_c < cost_m):
            	b = m
		m = c
                exec_time_m = exec_time_c 
                exec_time_b = exec_time_m
		cost_b = cost_m 
            	cost_m = cost_c 
	    elif (cost_d < cost_m):
            	a = m
		m = d
                exec_time_a = exec_time_m 
                exec_time_m = exec_time_d 
		cost_a = cost_m 
            	cost_m = cost_d 
	    else:
		print "errorrrrrrrrrrrrrrr"
	    	return
	if LA.norm(a - b) < 0.001: 
	 	#runs.append(np.append(np.append(np.append(a,exec_time_a),sim_count),cost_a))
		#runs.append(np.append(np.append(np.append(b,exec_time_b),sim_count),cost_b))
		runs.append(np.append(np.append(np.append(m,exec_time_m),sim_count),cost_m))
		#runs.append(np.append(np.append(np.append(m,exec_time_m),sim_count),cost_m))
		print "---------------------------------------------------------\n"
		print "Greedy Search Method\n"
		print "Best configuration\n"
		print "[Latency, Bandwidth, Execution time, Simulation runs, Cost saving]\n"
		print runs[-1]
		return runs
        
        
	
	
    
    
def find_initial_search_range(runs, first_point_on_boundry, w_init, exec_time_init, sim_count,topo,shape, tolerance):
	trace = []
	trace.append(first_point_on_boundry[-1])
	
	
	
	w = first_point_on_boundry[-1][:2]
	interval = (w_init - w) / 100.0
	N[1] = [w[0], w[1] / 2.0]
	w_diff = (N[1] - w) 
	old_cost = first_point_on_boundry[-1][4]
	cost_decrease = 1
	while cost_decrease > 0:
		sigma =1
		curr_error=0
		w = [w[0] + interval[0], w[1]]
		while (curr_error < tolerance) and (sigma > 0.0001):
			w_try = w + sigma * w_diff
			exec_time = run_sst(w_try,topo,shape)
			sim_count += 1
			curr_error = (exec_time - exec_time_init) / exec_time_init
			while (sigma > 0.0001) and (curr_error >= tolerance):
				sigma = sigma / 3.0
				w_try = w + sigma * w_diff
				exec_time = run_sst(w_try,topo,shape)
				sim_count += 1				
				curr_error = (exec_time - exec_time_init) / exec_time_init
			w = w_try
		cost = cost_model(w,w_init)
		cost_decrease =  old_cost - cost;
		old_cost = cost 	
		trace.append(np.append(np.append(np.append(w,exec_time),sim_count),cost))
 
	
	print "trace length"
	print len(trace)	
	
	if len(trace)==2:
		trace.append(first_point_on_boundry[-1])
		w = first_point_on_boundry[-1][:2]
		N[1] = [w[0]*2,w[1]]
		w_diff = (N[1]-w) 
		old_cost = first_point_on_boundry[-1][4]
		cost_decrease = 1
		while cost_decrease > 0:
			sigma =1
			curr_error=0
			w = [w[0],w[1]+interval[1]]
			while (curr_error < tolerance) and (sigma > 0.0001):
				w_try = w + sigma * w_diff
				exec_time = run_sst(w_try,topo,shape)
				sim_count += 1
				curr_error = (exec_time - exec_time_init) / exec_time_init
				while (sigma > 0.0001) and (curr_error >= tolerance):
					sigma = sigma / 3.0
					w_try = w + sigma * w_diff
					exec_time = run_sst(w_try,topo,shape)
					sim_count += 1
					curr_error = (exec_time - exec_time_init) / exec_time_init
				w = w_try		
			cost = cost_model(w,w_init)
			cost_decrease =  old_cost - cost;
			old_cost = cost 
			trace.append(np.append(np.append(np.append(w,exec_time),sim_count),cost))
        
	
	data = trace[-3:]
	print data
	
	#data1 = sorted(data[:2], key=lambda x: x[0])
	runs1 = find_optimal_point(data,w_init,topo,shape,runs,exec_time_init, sim_count, tolerance)
	return runs1



			
def greedy_search(runs, w_init, sim_count, topo, shape, tolerance):
	
	exec_time_init= run_sst(w_init, topo, shape)
	exec_time = exec_time_init
	cost = cost_model(w_init,w_init)
	sigma = 1.0000
	curr_error = 0
	w = w_init
	while (curr_error < tolerance) and (sigma>0):
		N[0] = [w[0]*2,w[1]]
		N[1] = [w[0],w[1]/2]
		cost_reduction_per_time_reduction = np.zeros(npop) #stores reward for each optimization direction
		reward = np.zeros(npop)
  		
		# perform the parameter update. The multiply below
  		# is to sum up all the rows of the noise matrix N,
		# where each row N[j] is weighted by standardize reward in direction j
		for j in range(npop):
    	      		time_reduction = run_sst(N[j],topo,shape) - exec_time + 1e-16#change w in direction j and run simulation for current change
			cost_reduction = cost - cost_model(N[j],w_init) + 1e-16
			cost_reduction_per_time_reduction[j] = cost_reduction/time_reduction
			print cost_reduction_per_time_reduction[j]
			sim_count += 1
		# standardize the rewards
		
		reward_accum = cost_reduction_per_time_reduction[0]+cost_reduction_per_time_reduction[1]
		if reward_accum == 0:
			break;
		reward[0] = cost_reduction_per_time_reduction[0]/reward_accum
		reward[1] = cost_reduction_per_time_reduction[1]/reward_accum
		
		
		w_diff = reward[1] * (N[1]-w) + reward[0] * (N[0]-w)
		w_try = w + sigma * w_diff
		exec_time = run_sst(w_try,topo,shape)
		sim_count += 1
			
			
		curr_error = (exec_time - exec_time_init) / exec_time_init
			
				
				
		while (sigma > 0.00001) and (curr_error >= tolerance):
			sigma = sigma/3.0
			w_try = w + sigma * w_diff
			exec_time = run_sst(w_try,topo,shape)
			sim_count += 1
			#runs.append(np.append(np.append(np.append(np.append(w,exec_time),count),0),color))
			curr_error = (exec_time - exec_time_init) / exec_time_init
		w = w_try
		cost = cost_model(w,w_init)
		#runs.append(np.append(np.append(np.append(w,exec_time),sim_count),cost))
		
				
			
	#cost_saving = cost_saving_model(w,w_init)
	
	#runs.append(np.append(np.append(np.append(w,exec_time),sim_count),cost))
	#runs.append(np.append(np.append(np.append(w,exec_time),sim_count),cost))
	
	
	first_point_on_boundry =[]
	first_point_on_boundry.append(np.append(np.append(np.append(w,exec_time),sim_count),cost))
	runs1 = find_initial_search_range(runs, first_point_on_boundry, w_init, exec_time_init, sim_count,topo,shape, tolerance)
	return runs1

def cost_model(w,w_init):
	return 1e-6/w[0] + w[1]
def exhaustive_search(w_init, sim_count, topo, shape):
	'''
    	exhaustive search
    	to find the boundry points on a grid
    	w_init: initial point to start simulation (best possible preformance)
	sim_count: count the number of run_sst calls
        '''
	runs=[] # stores simulation input and result
	exec_time_init = run_sst(w_init,topo,shape); # execution time for best possible configuration
	
	directions = 100 # number of direction we want to try in exhaustive search
	cost = 0.00
	for d in range(directions+1):
		cost = float(d) / float(directions)
		w = w_init
		sigma = 1.0000
		curr_error = 0
		while (curr_error < float(error.value)) and (sigma>0):
			N[0] = [w[0]*2,w[1]]
			N[1] = [w[0],w[1]/2]
			
		
			w_diff = cost * (N[1]-w) + (1.00 - cost) * (N[0]-w)
			w_try = w + sigma * w_diff
			exec_time = run_sst(w_try,topo,shape)
			sim_count += 1
			
			curr_error = (exec_time - exec_time_init) / exec_time_init
			
			while (sigma > 0.0001) and (curr_error >= float(error.value)):
				sigma = sigma/3.0
				w_try = w + sigma * w_diff
				exec_time = run_sst(w_try,topo,shape)
				sim_count += 1
				curr_error = (exec_time - exec_time_init) / exec_time_init
			w = w_try	
				
			
		cost = cost_model(w,w_init)
		runs.append(np.append(np.append(np.append(w,exec_time),sim_count),cost))
	min_index = min(enumerate(runs), key=(lambda x: float(x[1][4])))[0]
	print "---------------------------------------------------------\n"
	print "Exhaustive Search Method (found %d boundry point)\n" %(directions)
	print "Best configuration\n"
	print "[Latency, Bandwidth, Execution time, Simulation runs, Cost saving]\n"
	print runs[min_index]
	return runs
	
	
	
	
	
def update():
	
	w_init = np.array([1e-9, 200],dtype=np.float64) #initial point to start simulation
	
	sim_count = 0 # count number of simulations 
        #runs = exhaustive_search(w_init, 0,topo,shape)
	#find_optimal_point(runs, w_init, npop, sigma,0,cost,sim_count,topo,shape,"black")
	tolerance = 0.1
	runs=[]	
	for i in range(5):
		runs = greedy_search(runs, w_init, 0, topo, shape, 0.1)
		w_init = runs[-1][:2]
		tolerance += 0.1
	runs.append(runs[-1])
	with open('results.csv', 'wb') as myfile:
    		wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
    		for run in runs:
    			wr.writerow(run)


    

update()

