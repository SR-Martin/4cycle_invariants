#!/usr/bin/env python

import sys, getopt, errno, time
import numpy as np
import itertools
import random
import skbio
from skbio import TabularMSA, DNA

start = time.time()

msa_length = 1000
outputFilename = "msa.phylip"

edges = ["a", "b", "c", "d", "e", "f", "g", "h"]
nucleotides = ["A", "C", "G", "T"]
substitutionRates = dict()
generageEdgeRates = True
generateGamma = True
rate = "random"
networkString = "(0,1,2,3)"

try:
	opts, args = getopt.getopt(sys.argv[1:],"hs:p:g:l:o:r:n:")
except getopt.GetoptError:
	print("Option not recognised.")
	print("python simulate.py -l <MSA length> -s <seed> -g <gamma parameter> -o <output filename> -n <network string>")
	print("python simulate.py -h for further usage instructions.")
	sys.exit(2)
for opt, arg in opts:
	if opt == "-h":
		print("python simulate.py -l <MSA length> -s <seed> -p <edge parameters> -g <gamma parameter> -o <output filename>")
		print("-l <MSA length>\t\t Length of MSA to output")
		print("-s <seed>\t\t Integer value for seeding random numbers")
		print("-p <edge parameters>\t\t Comma-separated list of values giving the probability of mutation along each edge.")
		print("-g <gamma parameter>\t\t Floating point value giving the probability of a site evolving along edge e in the network.")
		print("-o <output filename>\t Filename for output MSA in phylip format")
		sys.exit()
	elif opt in ("-o"):
		outputFilename = arg 
	elif opt in ("-l"):
		msa_length = int(arg)
	elif opt in ("-s"):
		random.seed(int(arg))
	elif opt in ("-g"):
		val = float(arg)
		if val < 0 or val > 1:
			print ("Error: User supplied gamma parameter is not in interval [0,1]: " + param)
			sys.exit(2)
		treeRate = val
		generateGamma = False
	elif opt in ("-p"):
		parameters = arg.split(",")
		if len(parameters) != 8*3:
			print ("Error: User must supply all edge paramters or none.")
			sys.exit(2)
		for i in range(8):
			edge = edges[i]
			substitutionRates[edge] = [0,0,0]
			for j in range(3):
				val = float(parameters[3 * i + j])
				if val < 0 or val > 1:
					print ("Error: User supplied edge parameter is not in interval [0,1]: " + param)
					sys.exit(2)
				substitutionRates[edge][j]  = val
			if substitutionRates[edge][0] + substitutionRates[edge][2] + 2 * substitutionRates[edge][1] != 1.0:
				print ("Error: User supplied edge parameters for edge " + str(i) + " do not follow K2P format.")
				sys.exit(2)
		generageEdgeRates = False	
	elif opt in ("-r"):
		rate = arg
	elif opt in ("-n"):
		networkString = arg

networkPermutationString = networkString.strip("()").split(",")
networkPermutation = [int(leaf) for leaf in networkPermutationString]
inversePermutation = np.argsort(networkPermutation)

if generateGamma:
	treeRate = random.random()

if generageEdgeRates:
	if rate == "mixed":
		for edge in ["a","b","d","e","g"]:
			alpha = random.uniform(0.95, 0.99)
			beta = random.uniform(0,0.01)
			gamma = random.uniform(0,0.01)
			total = alpha + 2*beta + gamma
			substitutionRates[edge] = [alpha/total, beta/total, gamma/total]
		for edge in ["c","f","h"]:
			alpha = random.uniform(0.85, 0.95)
			beta = random.uniform(0,0.05)
			gamma = random.uniform(0,0.05)
			total = alpha + 2*beta + gamma
			substitutionRates[edge] = [alpha/total, beta/total, gamma/total]
	else:
		for edge in edges:
			if rate == "very_low":
				alpha = random.uniform(0.99, 1)
				beta = random.uniform(0,0.0025)
				gamma = random.uniform(0,0.0025)				
			elif rate == "low":
				alpha = random.uniform(0.95, 0.99)
				beta = random.uniform(0,0.01)
				gamma = random.uniform(0,0.01)
			elif rate == "medium":
				alpha = random.uniform(0.85, 0.95)
				beta = random.uniform(0,0.05)
				gamma = beta = random.uniform(0,0.05)		
			else:
				alpha = random.random()
				beta = random.random()
				gamma = random.random()
			total = alpha + 2*beta + gamma
			substitutionRates[edge] = [alpha/total, beta/total, gamma/total]

print("Simulating alignments with parameters: ")
print("Tree Ratio = " + str(treeRate))
print("\nedge\talpha\tbeta\tgamma")
for edge in edges:
	edge_params = str(edge)
	for param in substitutionRates[edge]:
		edge_params += "\t" + str(param)
	print(edge_params)

sequence_0 = ""
sequence_1 = ""
sequence_2 = ""
sequence_3 = ""

def mutate_K2P(start_nucl, edge_param):
	rand = random.random()
	if start_nucl == "A":
		if rand < edge_param[0]:
			return "A"
		elif rand < edge_param[0] + edge_param[2]:
			return "G"
		elif rand < edge_param[0] + edge_param[2] + edge_param[1]:
			return "C"
		else:
			return "T"
	elif start_nucl == "G":
		if rand < edge_param[2]:
			return "A"
		elif rand < edge_param[2] + edge_param[0]:
			return "G"
		elif rand < edge_param[2] + edge_param[0] + edge_param[1]:
			return "C"
		else:
			return "T"
	elif start_nucl == "C":
		if rand < edge_param[1]:
			return "A"
		elif rand < edge_param[1] + edge_param[1]:
			return "G"
		elif rand < edge_param[1] + edge_param[1] + edge_param[0]:
			return "C"
		else:
			return "T"
	else: # start_nucl == "T":
		if rand < edge_param[1]:
			return "A"
		elif rand < edge_param[1] + edge_param[1]:
			return "G"
		elif rand < edge_param[1] + edge_param[1] + edge_param[2]:
			return "C"
		else:
			return "T"

for i in range(msa_length):

	# start with a uniform distribution at the root
	rv_root = nucleotides[random.randrange(4)]

	rv_4 = mutate_K2P(rv_root, substitutionRates["h"])
	rv_6 = mutate_K2P(rv_4, substitutionRates["g"])
	rv_5 = 0
	if random.random() < treeRate:
		rv_5 = mutate_K2P(rv_root, substitutionRates["e"])
	else:
		rv_5 = mutate_K2P(rv_6, substitutionRates["f"])

	sequence_0 += mutate_K2P(rv_5, substitutionRates["a"])
	sequence_1 += mutate_K2P(rv_6, substitutionRates["b"])
	sequence_2 += mutate_K2P(rv_4, substitutionRates["c"])
	sequence_3 += mutate_K2P(rv_root, substitutionRates["d"])

seqs = [DNA(sequence_0, metadata={"id":"Taxon" + str(networkPermutation[0])}), 
		DNA(sequence_1, metadata={"id":"Taxon" + str(networkPermutation[1])}), 
		DNA(sequence_2, metadata={"id":"Taxon" + str(networkPermutation[2])}), 
		DNA(sequence_3, metadata={"id":"Taxon" + str(networkPermutation[3])})]

permutedSeqs = [seqs[inversePermutation[0]],seqs[inversePermutation[1]], seqs[inversePermutation[2]], seqs[inversePermutation[3]]]
aln = TabularMSA(permutedSeqs, minter="id")

f = open(outputFilename, "w")
aln.write(f, format='phylip')
f.close()

end = time.time()
print("simulate.py time: " + str(end - start))


