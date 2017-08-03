#! /usr/bin/Python
import os, sys, re
import subprocess 
import math
import csv
import numpy as np 
from numpy import linalg as LA
from time import sleep 
import argparse
from subprocess import Popen, PIPE
class buildingBlocks(object):
	"""Helper class for global settings."""
	costPerNode = 3000
	nic_switch_cable_cost=50
	s1_s2_cable_cost=200
	switches = [12, 24, 36, 48, 60]
	switch_cost = {12 : 160, 24 : 160*2, 36 : 160*3, 48:160*4, 60:160*5 }
	totalBudget = 5000000

def run_sst(topo,shape,numNodes,motif):
	'''the function runs the sst simulation for current set-up
	   @return exec_time: the fuction return the execution time for current simulation
	'''
	
	option = (
	'\'--cmdLine=Init --cmdLine=\"%s\" --cmdLine=Fini ' 
	'--numNodes=%d '
	'--numRanksPerNode=1 '
	'--topo=%s '
	'--shape=%s '
	'--width=1x1 '
	'--netBW=25GB/s '
        '--input_latency=30ns '
	'--rtrArb=merlin.xbar_arb_lru ' 
	'--emberRankmapper=fff\'') % (motif, numNodes, topo, shape)
    	command = 'mpirun -np 60 sst --model-options={0} ../../script/emberLoadCmd.py'.format(option)
	#command = " --model-options=%s ../../script/emberLoadCmd.py" % (option)	
    	print command
    	#os.system(command)
	p = subprocess.Popen(command, stdout=PIPE, shell=True)
	exec_time = -1
	for line in iter(p.stdout.readline,""):
     
		m = re.search('(?<=###execution_time )\d+\.\d+', line)
		if m:
			exec_time = float(m.group(0))
	#with open('output.csv', 'r') as f: #the execution time of sst simulation is stored in output.csv
	#	exec_time = float(f.readline())
	if exec_time == -1:
		print ("sst call was not successful");
		exit(2)
	
	return exec_time
class dragonfly_optimization(object):
	def build_optimal_dragonfly(self,motif, networkBudget, numNodes, npRow, optimal_configurations):
		"""Given number of compute nodes(i.e. numNodes) and availbe budget for cables and routers (i.e. network budget), 
                   search for a dragonfly configuration with highest bisection bandwidth which supports 
                   minimum number of compute nodes and cost of system is less than availble budget.
                """

		confs = []

		#brute-force search for a,h,p value in dragonlfy
		#p: number of compute node connected to each router
		#a: number of routers in each group
		#h: number of channels within each router used to connect to other groups
		for h in range(1,61):
			for p in range(1,61-h):
				a = int(math.ceil(float(math.sqrt(p*p+4*h*p*numNodes)-p)/float(2*h*p)))
				numGrps = (a * h + 1) #number of groups
				numRtrs = (a * h + 1) * a #number of routers
				numPorts = a + h + p - 1 #number of required ports 
				numMaxNodes = numRtrs * p
				if (numMaxNodes >= numNodes) and (numPorts <= 60):
			
					r=0
					networkCost = 0
					while r < 5:
						if buildingBlocks.switches[r] >= numPorts:
							break;
						r = r + 1
					radix = buildingBlocks.switches[r]
					routerCost = numRtrs * buildingBlocks.switch_cost[radix]
					cableCost = (numNodes * buildingBlocks.nic_switch_cable_cost)+((a*(a-1))/2)*numGrps*buildingBlocks.s1_s2_cable_cost+((numGrps*(numGrps-1))/2)*buildingBlocks.s1_s2_cable_cost
                                	networkCost = routerCost + cableCost

					if (networkCost <= networkBudget):
				
						bisection_d1 = float((a*h+1)*(a*h+1))/float(4*p) #bisection bandwidth in dimension 1
						bisection_d2 = float((a*a)*(a*h+1))/float(4*p)	#bisection bandwidth in dimension 2
						bisection = min(bisection_d1,bisection_d2) #minimum bisection bandwidth
						totalCost = networkCost + numNodes *  buildingBlocks.costPerNode
						conf = (radix, p,h,a,totalCost,bisection)
						confs.append(conf)

		#sort configurations based on bisection bandwidth
		sconfs = sorted(confs, key=lambda tup: tup[5])
			
		for c in sconfs[-1:]:	
			exe_t=0	
			shape = str(c[0]) +":"+str(c[1]) +":"+str(c[2]) +":"+str(c[3])
			new_motif = motif + " npRow=%d" %(npRow)	
			exe_t = run_sst("dragonfly" , shape, numNodes,new_motif)
			configuration = {'topo': "dragonfly", 'shape': shape, 'cost':c[4], 'numNodes':numNodes, 'execution_time':exe_t}
			optimal_configurations.append(configuration)
	def build_dragonfly(self, motif, numNodes, npRow, optimal_configurations):
		"""Increase the number of compute nodes in each iteration and
                   Search to build the dragonfly configuration with highest bisection bandwidth for each switch type on availble budget
                """

		scale = 1.1
                nodesCost = buildingBlocks.costPerNode * numNodes
                networkBudget = buildingBlocks.totalBudget - nodesCost
                while networkBudget > 0 and (int(scale*npRow) > int(npRow)):
                	self.build_optimal_dragonfly(motif, networkBudget, numNodes, npRow, optimal_configurations)
	
                     	old_npRow = npRow
                     
			npRow = int(math.ceil(npRow * scale))
                        numNodes = npRow * npRow
                        nodesCost = buildingBlocks.costPerNode * numNodes
                        networkBudget = buildingBlocks.totalBudget - nodesCost
                        while (networkBudget < 0) and (int(scale*npRow) > int(npRow)):
                                npRow = old_npRow
                                scale = float(1 + scale)/float(2)
                                npRow = int(math.ceil(npRow * scale))
                                numNodes = npRow * npRow
                                nodesCost = buildingBlocks.costPerNode * numNodes
                                networkBudget = buildingBlocks.totalBudget - nodesCost

class fattree_optimization(object):
	def build_shape(self, down, up, switch, numNodes, level):
		"""
		Given the shape of fatree for levels lower than level l, bulid the smallest fattree that supports numNodes 
		"""
		# build fatree with highest oversubscription for levels > l
		radix = switch / 2	
		hosts = numNodes
		for l in range(0,level+1):
        		hosts = (hosts + down[l] - 1) / down[l]
		if hosts== 1:
			up = up[:-1]
		else:
        		while hosts > radix:
        			down.append(radix)
                		up.append(radix)
                		level += 1
                		hosts = (hosts + radix - 1) / radix
        		down.append(int(hosts))
			level += 1
		maxHost = 1
		for d in down:
			maxHost *= d
		
		l = level
		# for each level check if reducing the number of downlinks is possible or not
		while l >= 0:
			i = down[l] - 1
			maxHost = 1
	        	for d in down:
        	        	maxHost *= d
			while (maxHost/down[l])*i > numNodes:
				i -= 1
			down[l] = i + 1
			if l < level:
				if up[l] > down[l]:
					up[l] = down[l]
			l -= 1
     
		shape =""	
		#based on downs and ups build shape string (shape is an input to sst)
		for l in range(0,level):
			shape += "%d,%d:" % (down[l], up[l])
		shape += "%d" %(down[level])

        	return shape, level
	def fattree_switches_cost(self, downs, ups, switch, numNodes):
		"""
		given the shape of a fattree (ups,downs), calculate the cost of routers and cables
		"""
		routers_per_level = []  # keeps the number of routers in each level of fattree
		maxNodes = 1
		levels = 0
		for d in downs:
			maxNodes *= d 
			levels += 1
    	
		routers_per_level.append(maxNodes / downs[0])
		cableCost = numNodes * buildingBlocks.nic_switch_cable_cost #cable cost of leaf level
		#calculate the number of routers and cable cost in each level
    		for i in range(1, levels):
			numRtrs = math.ceil((float(routers_per_level[i-1]) * float(ups[i-1])) / float(downs[i]))
			cableCost += numRtrs*downs[i] * buildingBlocks.s1_s2_cable_cost
			routers_per_level.append(numRtrs)
		total_switch_numbers = 0
		for r in routers_per_level:
			total_switch_numbers += r
		routerCost = total_switch_numbers * buildingBlocks.switch_cost[switch]
		networkCost = routerCost + cableCost
		return networkCost

 	
	def build_optimal_fattree(self, motif, networkBudget, numNodes, npRow, optimal_configurations):
		"""Given number of compute nodes(i.e. numNodes) and availbe budget for cables and routers (i.e. network budget), 
		   for each switch type, search for a fattree configuration with highest bisection bandwidth which supports 
		   minimum number of compute nodes and cost of system is less than availble budget.
                """
		confs=[]

		for switch in buildingBlocks.switches:
			down = [] #keeps number of downlinks in each level of fattree
			up = [] #keeps number of uplinks in each level of fattree
			radix = switch/2
			down.append(radix)
			up.append(switch - radix)
			shape, level = self.build_shape(down, up, switch, numNodes, 0)
                	networkCost = self.fattree_switches_cost(down, up, switch, numNodes)
                
			sub = np.empty(level); 
			sub.fill(radix) #
			l = 0
			while (networkCost > networkBudget):
				sub[l] = sub[l] - 1 # increase oversubscription ratio (downlinks/uplinks) for routers in level l
				if (sub[l] == 0):   # when the oversubscription ratio for level l reach its maximum value, set level to l+1
					sub[l] = 1
					l += 1;
					if l == level:
						break
				up = [1] * l
				down = [switch-1] * l 				
				down.append(switch - sub[l])
        			up.append(sub[l])


				shape, level2 = self.build_shape(down, up, switch, numNodes, l)
				networkCost = self.fattree_switches_cost(down, up, switch, numNodes)
		
			if networkCost <= networkBudget: # add optimal configuration to list
				totalCost = networkCost + numNodes * buildingBlocks.costPerNode
				configuration = {'topo': "fattree", 'shape': shape, 'cost':totalCost, 'numNodes':numNodes}
				confs.append(configuration)
	
		#run sst for all optimal configurations
		for conf in confs:
			exe_t = 0.0

			new_motif = motif + " npRow=%d" %(npRow) 

			exe_t = run_sst("fattree" , conf["shape"], numNodes,new_motif)
			conf['execution_time'] = exe_t	    
			optimal_configurations.append(conf)    	
		
	




	def build_fattree(self, motif, numNodes, npRow, optimal_configurations):
		"""Increase the number of compute nodes in each iteration and
		   Search to build the fattree configurations with highest bandsection bandwidth for each switch type on availble budget
		"""
		scale = 1.1 #scale of increaseing the number of compute nodes in each iteration
		nodesCost = buildingBlocks.costPerNode * numNodes 
		networkBudget = buildingBlocks.totalBudget - nodesCost
		while networkBudget > 0 and (int(scale*npRow) > int(npRow)):
        		self.build_optimal_fattree(motif, networkBudget, numNodes, npRow, optimal_configurations) 
               		old_npRow = npRow
                	npRow = int(math.ceil(npRow * scale))
			numNodes = npRow * npRow
                	nodesCost = buildingBlocks.costPerNode * numNodes
                	networkBudget = buildingBlocks.totalBudget - nodesCost      
                	while (networkBudget < 0) and (int(scale*npRow) > int(npRow)):
				npRow = old_npRow
				scale = float(1 + scale)/float(2)
				npRow = int(math.ceil(npRow * scale))
                        	numNodes = npRow * npRow
				nodesCost = buildingBlocks.costPerNode * numNodes
				networkBudget = buildingBlocks.totalBudget - nodesCost


def main():


	"""command-line interface to optimization 
	   Example: 
	   topology_optimizer.py --motif "FFT3D iteration=1 nx=64 ny=64 nz=64" --output_file=result.csv

	"""
        parser = argparse.ArgumentParser()
        parser.add_argument('--motif', type=str, required=True)
        parser.add_argument('--output_file', type=str, required=True)
        args = parser.parse_args()
        motif = args.motif
       	output_file = args.output_file

	#availble topologies
	topos = ["fattree", "dragonfly"]

	optimal_configurations = [] #keeps optimal configurations for all topology
	for topo in topos:
		scale = 1.1
	        npRow = 32 #fft3d motif argument
        	numNodes = 1024  #initial number of hosts
  
		
		if topo == 'fattree':
			fattree_optimization().build_fattree(motif, numNodes, npRow, optimal_configurations)
		if topo == 'dragonfly': 
			dragonfly_optimization().build_dragonfly(motif, numNodes, npRow, optimal_configurations)
 
	#write all optimal configurations to output file
	with open(output_file, 'w') as csvfile:
    		fieldnames = ['topo', 'shape', 'numNodes', 'cost', 'execution_time']
		writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
		writer.writeheader()
    		for conf in optimal_configurations:
			writer.writerow(conf)


				
if __name__ == "__main__":
        main()

