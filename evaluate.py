#!/usr/bin/env python

import sys, getopt, errno, time
import numpy as np
import itertools
import mpmath
import skbio
from skbio import DNA, TabularMSA
import utils.polynomials
from utils.polynomials import Polynomial
import operator
import os

start = time.time()

#Metric Functions:
def get1Norm(values):
	total = mpmath.mpf(0)
	for val in values:
		total = mpmath.fadd(total, mpmath.fabs(val))
	return total

def get1NormNormalised(values):
	total = get1Norm(values)
	return total/len(values)

def getMaxNorm(values):
	maxVal = 0
	for val in values:
		if abs(val) > maxVal:
			maxVal = abs(val)
	return maxVal

def getEuclideanNorm(values):
	total = 0
	for val in values:
		total += val*val
	return np.sqrt(total)

def getProdScore(values):
	return mpmath.fprod(values)

def getClimbingScore(networkString):
	score = 1
	for poly in networkPolys:
		polyValues = invariant_values[poly.getPolyString()]
		sortedPolyValues = sorted(polyValues.items(), key=operator.itemgetter(1))
		for i in range(0,len(sortedPolyValues)):
			if sortedPolyValues[i][0] == networkString:
				score *= i+1
				break
	return score

Nucl = {
  "A": 0, # = (0,0)
  "C": 1, # = (0,1)
  "G": 2, # = (1,0)
  "T": 3  # = (1,1)
}

# Function to return character values when calculating Fourier transform
def Chi(char, nucl):
	if char == 0: # Chi_A
		return 1
	elif char == 1: # Chi_C
		if nucl == 1 or nucl == 3:
			return -1
		else:
			return 1
	elif char == 2: # Chi_G
		if nucl == 2 or nucl == 3:
			return -1
		else:
			return 1
	else: # Chi_T
		if nucl == 1 or nucl == 2:
			return -1
		else:
			return 1

MSAFilename = ""
scoringFunction = get1Norm
climbingScore = False
useDir = False
directoryStr = ""

multiplierForTrees1 = 3
multiplierForTrees2 = 2

try:
	opts, args = getopt.getopt(sys.argv[1:],"ha:i:s:d:")
except getopt.GetoptError:
	print("Option not recognised.")
	print("python evaluate.py -a <MSA file> -d <directory> -i <invariants file> -s <scoring method>")
	print("python evaluate.py -h for further usage instructions.")
	sys.exit(2)
for opt, arg in opts:
	if opt == "-h":
		print("python evaluate.py -a <MSA file> -d <directory> -i <invariants file> -s <scoring method>")
		print("-a <MSA file>\t\t Multiple sequence alignment file.")
		print("-d <directory>\t\t Directory of MSAs to use. All MSAs must have taxa in the same order.")
		print("-i <invariants file>\t\t File containing list of polynomial invariants to use in Fourier coordinates.")
		print("-s <scoring method>\t\t Method of scoring networks from invariants.")
		sys.exit()
	elif opt in ("-a"):
		MSAFilename = arg
	elif opt in ("-i"):
		NetworkGBFilename = arg
	elif opt in ("-s"):
		if arg.lower() == "1norm":
			scoringFunction = get1Norm
		elif arg.lower() == "1normnormalised":
			scoringFunction = get1NormNormalised
		elif arg.lower() == "euclideannorm":
			scoringFunction = getEuclideanNorm
		elif arg.lower() == "maxnorm":
			scoringFunction = getMaxNorm
		elif arg.lower() == "product":
			scoringFunction = getProdScore
		elif arg.lower() == "climbing":
			climbingScore = True
			scoringFunction = getClimbingScore
		else:
			print("Error: Could not understand scoring method " + arg)
			sys.exit(2)
	elif opt in ("-d"):
		useDir = True
		directoryStr = arg

if (not useDir and len(MSAFilename) == 0) or (useDir and len(MSAFilename) > 0):
	print("Error: You must provide either a directory (-d) or a single file (-a).")
	sys.exit(2)

if useDir and len(directoryStr) == 0:
	print("Error: You must provide a directory with -d <directory>.")
	sys.exit(2)	

if useDir:
	# Open every file in the directory
	first = True
	directory = os.fsencode(directoryStr)
	for file in os.listdir(directory):
		if os.fsdecode(file).startswith('.'):
			continue
		filename = directoryStr + "/" + os.fsdecode(file)
		seqs = []
		with open(filename, 'r') as infile:
			for line in infile:
				if(line[0] == '>'):
					fields = line.split()
					seqs.append(fields[0][1:])
		if first:
			try:
				fullMSA = TabularMSA.read(filename, constructor=DNA)
				fullMSA.index = seqs
				fullMSA.sort()
			except (ValueError, TypeError, skbio.io.UnrecognizedFormatError) as e:
				print(e)
				sys.exit(2)
			first = False
		else:
			try:
				nextMSA = TabularMSA.read(filename, constructor=DNA)
				nextMSA.index = seqs
				nextMSA.sort()
				fullMSA = fullMSA.join(nextMSA)
			except (ValueError, TypeError, skbio.io.UnrecognizedFormatError) as e:
				print(e)
				sys.exit(2)
	fullMSA.write("all.fasta")
	print(fullMSA)
	if not fullMSA or len(fullMSA) < 4:
		print("Error: MSA files must be a multiple sequence alignment of at least 4 sequences.")
		sys.exit(2)
else:
	# Open the user-specified file
	seqs = []
	with open(MSAFilename, 'r') as infile:
		for line in infile:
			if(line[0] == '>'):
				fields = line.split()
				seqs.append(fields[0][1:])
	try:
		fullMSA = TabularMSA.read(MSAFilename, constructor=DNA)
		fullMSA.index = seqs
		fullMSA.sort()
		print(fullMSA)
	except (ValueError, TypeError, skbio.io.UnrecognizedFormatError) as e:
		print(e)
		sys.exit(2)

if not fullMSA or len(fullMSA) < 4:
	print("Error: MSA files must be a multiple sequence alignment of at least 4 sequences.")
	sys.exit(2)

# Read in polynomials
#print("Reading invariants...")
try:
	with open(NetworkGBFilename, 'r') as infile:
		line = infile.readline()
		line = line.strip("| ")
		fields = line.split()
		networkPolys = list()
		for polynomial in fields:
			if polynomial.split() != "0":
				polyObject = Polynomial(polynomial)
				if polyObject not in networkPolys:
					networkPolys.append(polyObject)
					#print("Added poly: " + polyObject.getPolyString())
except (OSError, IOError) as e: 
	if getattr(e, 'errno', 0) == errno.ENOENT:
		print("Could not find file " + NetworkGBFilename)
	print(e)
	sys.exit(2)

subsets = list(itertools.combinations(range(fullMSA.shape[0]), 4))
for subset in subsets:
	print("-----------------------------------------------------------------------")
	print("Evaluating subset {" + str(fullMSA.index[subset[0]]) + ", " + str(fullMSA.index[subset[1]]) + ", " + str(fullMSA.index[subset[2]]) + ", " + str(fullMSA.index[subset[3]]) + "}")
	originalMSA = fullMSA.iloc[[subset[0],subset[1],subset[2],subset[3]]]

	# Permute the sequences in MSA to test all possible leaf labellings 
	# (instead of calculating invariants for all possible leaf labellings)
	msas = list()
	#permutations = itertools.permutations([0,1,2,3])
	#permutations = [(0,1,2,3),(3,2,1,0)]
	permutations = [(0,1,2,3),(0,2,1,3),(0,1,3,2),(1,2,0,3),(1,0,2,3),(1,0,3,2),(2,1,0,3),(2,0,1,3),(2,0,3,1),(3,1,0,2),(3,0,1,2),(3,0,2,1)]
	for perm in permutations:
		seqs = [originalMSA.iloc[perm[0]], originalMSA.iloc[perm[1]], originalMSA.iloc[perm[2]], originalMSA.iloc[perm[3]]]
		ind = [originalMSA.index[perm[0]], originalMSA.index[perm[1]], originalMSA.index[perm[2]], originalMSA.index[perm[3]]]
		msas.append(TabularMSA(seqs, index=ind))

	networkNorms = dict()

	#print("Calculating scores...")
	invariant_values = dict()
	for poly in networkPolys:
		invariant_values[poly.getPolyString()] = dict()

	for MSA in msas:
		labelString = "(Leaf 0: " + str(MSA.index[0]) + ", Leaf 1: " + str(MSA.index[1]) + ", Leaf 2: " + str(MSA.index[2]) + ", Leaf 3: " + str(MSA.index[3]) + ")"
		dictionaryString = "("+ str(MSA.index[0]) + "," + str(MSA.index[1]) + "," + str(MSA.index[2]) + "," + str(MSA.index[3]) + ")"
		outputString = "("+ str(MSA.index[0])[-1] + str(MSA.index[1])[-1] + str(MSA.index[2])[-1] + str(MSA.index[3])[-1] + ")"
		
		for poly in networkPolys:
				invariant_values[poly.getPolyString()][outputString] = 0

		# The 2-cycle (13) gives the same network, so we would like the result to be the same in both cases
		# We achieve this by calculating the frequency scores for both networks.
		symmetricSeqs = [MSA.iloc[0], MSA.iloc[3], MSA.iloc[2], MSA.iloc[1]]
		symmetricInd = [MSA.index[0], MSA.index[3], MSA.index[2], MSA.index[1]]
		symmetricMSA = TabularMSA(symmetricSeqs, index=symmetricInd)

		for msa in [MSA, symmetricMSA]:

			# Create frequencies array
			frequencies = np.zeros(shape=(4,4,4,4), dtype=np.longdouble)
			count = 0
			for col in msa.iter_positions(ignore_metadata=True):
				if "-" not in col:
					frequencies[Nucl[str(col[0])], Nucl[str(col[1])], Nucl[str(col[2])], Nucl[str(col[3])]] += 1
					count += 1

			for i in range(4):
				for j in range(4):
					for k in range(4):
						for l in range(4):
							freq = frequencies[i,j,k,l]
							if freq != 0:
								frequencies[i,j,k,l] = mpmath.fdiv(freq, count)

			# Perform Fourier transformation
			transformed = np.zeros(shape=(4,4,4,4), dtype=np.longdouble)
			for i in range(4):
				for j in range(4):
					for k in range(4):
						for l in range(4):
							transformed_value = mpmath.mpf(0)
							for w in range(4):
								for x in range(4):
									for y in range(4):
										for z in range(4):
											if frequencies[w,x,y,z] > 0:
												transformed_value = mpmath.fadd(transformed_value, mpmath.fprod([Chi(i,w), Chi(j,x), Chi(k,y), Chi(l,z), frequencies[w,x,y,z]]))
							transformed[i,j,k,l] = mpmath.fdiv(transformed_value, mpmath.power(4,4))

			print(transformed)
			sys.exit(0)
			for poly in networkPolys:
				invariantSum = mpmath.fadd(invariant_values[poly.getPolyString()][outputString], mpmath.fabs(poly.evaluate_mpm(transformed)))
				invariant_values[poly.getPolyString()][outputString] = invariantSum

		# divide by the number of symmetric MSAs for each invariant, to get the mean abs.
		for poly in networkPolys:
			invariant_values[poly.getPolyString()][outputString] = mpmath.fdiv(invariant_values[poly.getPolyString()][outputString], 2)

	for MSA in msas:
		labelString = "(Leaf 0: " + str(MSA.index[0]) + ", Leaf 1: " + str(MSA.index[1]) + ", Leaf 2: " + str(MSA.index[2]) + ", Leaf 3: " + str(MSA.index[3]) + ")"
		dictionaryString = "("+ str(MSA.index[0]) + "," + str(MSA.index[1]) + "," + str(MSA.index[2]) + "," + str(MSA.index[3]) + ")"
		outputString = "("+ str(MSA.index[0])[-1] + str(MSA.index[1])[-1] + str(MSA.index[2])[-1] + str(MSA.index[3])[-1] + ")"		
		networkValues = list()

		for poly in networkPolys:
			val = invariant_values[poly.getPolyString()][outputString]
			networkValues.append(val)

		if not climbingScore:
			score = scoringFunction(networkValues)
		else:
			score = getClimbingScore(outputString)
		networkNorms[dictionaryString] = score

	sortedScores = sorted(networkNorms.items(), key=operator.itemgetter(1))

	for i in range(0,12):
		print(str(i+1) + ": " + sortedScores[i][0] + "\t" + mpmath.nstr(sortedScores[i][1]))

	#firstDiff = sortedScores[7][1] - sortedScores[0][1]
	#secondDiff = sortedScores[-1][1] - sortedScores[7][1]

	if sortedScores[7][1] < multiplierForTrees1 * sortedScores[0][1] and sortedScores[8][1] > multiplierForTrees2 * sortedScores[7][1]:
		network1 = sortedScores[8][0]
		taxa = network1.strip("()").split(",")
		network2 = "(" + taxa[2] + "," + taxa[1] + "," + taxa[0] + "," + taxa[3] + ")"
		network3_1 = "(" + taxa[1] + "," + taxa[0] + "," + taxa[3] + "," + taxa[2] + ")"
		network3_2 = "(" + taxa[1] + "," + taxa[2] + "," + taxa[3] + "," + taxa[0] + ")"
		network4_1 = "(" + taxa[3] + "," + taxa[0] + "," + taxa[1] + "," + taxa[2] + ")"
		network4_2 = "(" + taxa[3] + "," + taxa[2] + "," + taxa[1] + "," + taxa[0] + ")"
		networks = [network1, network2, network3_1, network3_2, network4_1, network4_2]
		if sortedScores[9][0] in networks and sortedScores[10][0] in networks and sortedScores[11][0] in networks:
			print("Tree-like evolution detected with tree ((" + taxa[0] + "," + taxa[2] + "),(" + taxa[1] + "," + taxa[3] + ")).")

end = time.time()
print("evaluate_mp.py time: " + str(end - start))


