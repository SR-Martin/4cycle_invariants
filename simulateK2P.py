#!/usr/bin/env python

import sys, getopt, errno, time
import numpy as np
import itertools
import random
import skbio
from skbio import TabularMSA, DNA

start = time.time()

msa_length = 10
model = "K2P"
outputFilename = "msa.phylip"

edges = ["a", "b", "c", "d", "e", "f", "g", "h"]
nucleotides = ["A", "C", "G", "T"]
substitutionRates = dict()
generageEdgeRates = True
generateGamma = True
rate = "random"

try:
	opts, args = getopt.getopt(sys.argv[1:],"hm:s:p:g:l:o:r:")
except getopt.GetoptError:
	print("Option not recognised.")
	print("python simulate.py -m <model> -l <MSA length> -s <seed> -g <gamma parameter> -o <output filename>")
	print("python simulate.py -h for further usage instructions.")
	sys.exit(2)
for opt, arg in opts:
	if opt == "-h":
		print("python simulate.py -m <model> -l <MSA length> -s <seed> -p <edge parameters> -g <gamma parameter> -o <output filename>")
		print("-m <model>\t\t One of JC, K2P, or K3P")
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
#	 elif opt in("-m"):
#	 	model = arg.upper()
#	 	if model == "JC":
#
#		elif model == "K2P":
#
#		elif model == "K3P":
#
#		else:
#			print("Could not understand model(-m) argument. Use either JC, K2P, or K3P.")
#			sys.exit(2)

if generateGamma:
	treeRate = random.random()

if generageEdgeRates:
	for edge in edges:
		if rate == "low":
			alpha = random.uniform(0.99, 1)
			gamma = random.uniform(0,0.001)
			beta = random.uniform(0,0.001)
		elif rate == "medium":
			alpha = random.uniform(0.85, 0.95)
			gamma = random.uniform(0,0.03)
			beta = random.uniform(0,0.01)
		else:
			alpha = random.random()
			gamma = random.random()
			beta = random.random()
		randomSum = alpha + gamma + 2 * beta
		substitutionRates[edge] = [alpha/randomSum, beta/randomSum, gamma/randomSum]

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
			return "C"
		elif rand < edge_param[0] + edge_param[2] + edge_param[1]:
			return "G"
		else:
			return "T"
	elif start_nucl == "C":
		if rand < edge_param[2]:
			return "A"
		elif rand < edge_param[2] + edge_param[0]:
			return "C"
		elif rand < edge_param[2] + edge_param[0] + edge_param[1]:
			return "G"
		else:
			return "T"
	elif start_nucl == "G":
		if rand < edge_param[1]:
			return "A"
		elif rand < edge_param[1] + edge_param[1]:
			return "C"
		elif rand < edge_param[1] + edge_param[1] + edge_param[0]:
			return "G"
		else:
			return "T"
	else: # start_nucl == "T":
		if rand < edge_param[1]:
			return "A"
		elif rand < edge_param[1] + edge_param[1]:
			return "C"
		elif rand < edge_param[1] + edge_param[1] + edge_param[2]:
			return "G"
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

seqs = [DNA(sequence_0, metadata={"id":"Taxon0"}), 
				DNA(sequence_1, metadata={"id":"Taxon1"}), 
				DNA(sequence_2, metadata={"id":"Taxon2"}), 
				DNA(sequence_3, metadata={"id":"Taxon3"})]
aln = TabularMSA(seqs, minter="id")

f = open(outputFilename, "w")
aln.write(f, format='phylip')
f.close()

end = time.time()
print("simulate.py time: " + str(end - start))


