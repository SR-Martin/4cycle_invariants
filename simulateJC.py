#!/usr/bin/env python

import sys, getopt, errno, time
import numpy as np
import itertools
import random
import skbio
from skbio import TabularMSA, DNA

start = time.time()

msa_length = 1000
model = "JC"
outputFilename = "msa.phylip"

edges = ["a", "b", "c", "d", "e", "f", "g", "h"]
nucleotides = ["A", "C", "G", "T"]
substitutionRates = dict()
generageEdgeRates = True
generateGamma = True
biologyModel = False
rate = ""

try:
	opts, args = getopt.getopt(sys.argv[1:],"hm:s:p:g:l:o:br:")
except getopt.GetoptError:
	print("Option not recognised.")
	print("python simulate.py -m <model> -l <MSA length> -s <seed> -p <edge parameters> -g <gamma parameter> -o <output filename>")
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
	elif opt in ("-p"):
		parameters = arg.split(",")
		if len(parameters) != 8:
			print ("Error: User must supply all edge paramters or none.")
			sys.exit(2)
		i = 0
		for param in parameters:
			val = float(param)
			if val < 0 or val > 1:
				print ("Error: User supplied edge parameter is not in interval [0,1]: " + param)
				sys.exit(2)
			substitutionRates[edges[i]] = val
			i += 1
		generageEdgeRates = False
	elif opt in ("-g"):
		val = float(arg)
		if val < 0 or val > 1:
			print ("Error: User supplied gamma parameter is not in interval [0,1]: " + param)
			sys.exit(2)
		gamma = val
		generateGamma = False
	elif opt in ("-b"):
		biologyModel = True	
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
	gamma = random.random()
if generageEdgeRates:
	for edge in edges:
		if rate == "very_low":
			substitutionRates[edge] = random.uniform(0,0.001)
		elif rate == "low":
			substitutionRates[edge] = random.uniform(0,0.01)
		else:
			substitutionRates[edge] = random.random()

if biologyModel:
	for edge in edges:
		length = substitutionRates[edge]
		transitionProb = (3 - 3 * np.exp(-4 * length /3.))/4.

simType = "from transition matrices"
if biologyModel:
	simType = "from rate matrices"
print("Simulating alignments " + simType + " with parameters: ")
print("Tree ratio: " + str(gamma))
print("\nedge\talpha")
print("-------------------------------")
for edge in edges:
	edge_params = str(edge) + "\t" + str(substitutionRates[edge])
	print(edge_params)


sequence_0 = ""
sequence_1 = ""
sequence_2 = ""
sequence_3 = ""

def mutate_JC(start_nucl, edge_param):
	rand = random.random()
	if rand > edge_param:
		return start_nucl
	else:
		nucl = ["A", "C", "G", "T"]
		nucl.remove(start_nucl)
		return nucl[random.randrange(3)]

for i in range(msa_length):

	# start with a uniform distribution at the root
	rv_root = nucleotides[random.randrange(4)]

	rv_4 = mutate_JC(rv_root, substitutionRates["h"])
	rv_6 = mutate_JC(rv_root, substitutionRates["g"])
	rv_5 = 0
	if random.random() < gamma:
		rv_5 = mutate_JC(rv_4, substitutionRates["e"])
	else:
		rv_5 = mutate_JC(rv_6, substitutionRates["f"])

	sequence_0 += mutate_JC(rv_5, substitutionRates["a"])
	sequence_1 += mutate_JC(rv_6, substitutionRates["b"])
	sequence_2 += mutate_JC(rv_root, substitutionRates["c"])
	sequence_3 += mutate_JC(rv_4, substitutionRates["d"])

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


