#!/usr/bin/env python

import numpy as np
import argparse
from collections import defaultdict
import pandas as pd
from tqdm import tqdm
import concurrent.futures
from concurrent import futures

def simulate(N, p=0.001):
	'''
	N is the coverage to simulate

	p is the probability of an illumina error

	Returns three probabilities, which are the probabilities of error for the three bases
	'''

	#Assume "A" is the correct nucleotide.
	result = np.random.choice([0,1,2,3], size=N, p=[1-p,p/3,p/3,p/3])
	num_G = np.sum(result == 1)
	num_C = np.sum(result == 2)
	num_T = np.sum(result == 3)

	nucleotides = [num_G,num_C,num_T]
	return nucleotides

def full_simulate(N, bootstraps, p=0.001):
	err_obs = defaultdict(int)

	for bt in range(1, bootstraps):
		# Do simulation
		res = simulate(N, p=p)
		err_obs[res[0]] += 1
		err_obs[res[1]] += 1
		err_obs[res[2]] += 1


	table = defaultdict(list)
	running_total=0
	for i in reversed(range(1, max(list(err_obs.keys())) + 1)):
		running_total += err_obs[i]

		table['coverage'].append(N)
		table['error_base_coverage'].append(i)
		table['probability'].append(running_total/bootstraps)

	return pd.DataFrame(table)

def main(args):
	'''
	Run the simulation
	'''
	bootstraps = args.bootstraps
	coverage = args.max_coverage
	probability = args.probability_of_error
	threads = args.threads

	ex = concurrent.futures.ProcessPoolExecutor(max_workers=threads)
	wait_for = [ex.submit(prob_wrapper, cmd) for cmd in iterate_commands(coverage, bootstraps, probability)]
	total_cmds = coverage

	dbs = []
	for f in tqdm(futures.as_completed(wait_for), total=total_cmds, desc='Calculating null'):
		db = f.result()
		dbs.append(db)


	db = pd.concat(dbs).sort_values(['coverage', 'error_base_coverage'])
	db = db.pivot(index='coverage', columns='error_base_coverage', values='probability').fillna(0).reset_index()
	db.to_csv('NullModel.txt', sep='\t', index=False)

def prob_wrapper(cmd):
	N, bootstraps, probability = cmd
	return full_simulate(N, bootstraps, p=probability)

def iterate_commands(coverage, bootstraps, probability):
	for N in range(1,coverage, 1):
		cmd = (N, bootstraps, probability)
		yield cmd

if __name__ == '__main__':
	 parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	 parser.add_argument("-c", "--max_coverage", help="maximum coverage to profile", type=int, default=10000)
	 parser.add_argument("-b", "--bootstraps", help="number of bootstraps", type=int, default=10000)
	 parser.add_argument("-p", "--probability_of_error", help="probability of a sequencing error", type=float, default=0.001)
	 parser.add_argument("-t", "--threads", help="threads", type=int, default=6)

	 args = parser.parse_args()
	 main(args)
