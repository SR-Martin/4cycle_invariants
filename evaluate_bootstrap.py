#!/usr/bin/env python

import sys, getopt, errno, time
import numpy as np
import mpmath
from skbio import DNA, TabularMSA
from multiprocessing import Pool
import utils.polynomials
from utils.polynomials import Polynomial
import operator
import random

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

def getFrequenciesFromPermutation(frequencyArray, permutation):
	returnArray = np.zeros(shape=(4,4,4,4), dtype=np.longdouble)
	for i in range(4):
		for j in range(4):
			for k in range(4):
				for l in range(4):
					index=(i,j,k,l)
					permutedIndex=(index[permutation[0]], index[permutation[1]], index[permutation[2]], index[permutation[3]])
					returnArray[permutedIndex] = frequencyArray[index]
	return returnArray

def evaluateBootstrap(MSA, invariants, numSamples, seed):
	#print("Sampling " + str(x+1) + " of " + str(numberOfBootstraps))
	invariant_values = dict()
	permutations = [(0,1,2,3),(0,2,1,3),(0,1,3,2),(1,2,0,3),(1,0,2,3),(1,0,3,2),(2,1,0,3),(2,0,1,3),(2,0,3,1),(3,1,0,2),(3,0,1,2),(3,0,2,1)]
	random.seed(seed)
	allScores = dict()
	multiplierForTrees1 = 3
	multiplierForTrees2 = 2

	for poly in invariants:
		invariant_values[poly.getPolyString()] = dict()

	# Create frequencies array from MSA
	originalFrequencies = np.zeros(shape=(4,4,4,4), dtype=np.longdouble)
	count = 0
	while count < numSamples:
		i = random.randrange(MSA.shape[1])
		col = []
		for j in range(4):
			col.append(str(MSA[j,i]))
		if "-" not in col:
			originalFrequencies[Nucl[col[0]], Nucl[col[1]], Nucl[col[2]], Nucl[col[3]]] += 1
			count += 1

	for i in range(4):
		for j in range(4):
			for k in range(4):
				for l in range(4):
					freq = originalFrequencies[i,j,k,l]
					if freq != 0:
						originalFrequencies[i,j,k,l] = mpmath.fdiv(freq, count)

	originalTransformed = np.zeros(shape=(4,4,4,4), dtype=np.longdouble)
	for i in range(4):
		for j in range(4):
			for k in range(4):
				for l in range(4):
					transformed_value = mpmath.mpf(0)
					for w in range(4):
						for x in range(4):
							for y in range(4):
								for z in range(4):
									if originalFrequencies[w,x,y,z] > 0:
										transformed_value = mpmath.fadd(transformed_value, mpmath.fprod([Chi(i,w), Chi(j,x), Chi(k,y), Chi(l,z), originalFrequencies[w,x,y,z]]))
					originalTransformed[i,j,k,l] = mpmath.fdiv(transformed_value, mpmath.power(4,4))

	for perm in permutations:
		dictionaryString = "("+ str(MSA.index[perm[0]]) + "," + str(MSA.index[perm[1]]) + "," + str(MSA.index[perm[2]]) + "," + str(MSA.index[perm[3]]) + ")"
		
		for poly in invariants:
			invariant_values[poly.getPolyString()][dictionaryString] = 0

		# The 2-cycle (24) gives the same network, so we would like the result to be the same in both cases
		# We achieve this by calculating the frequency scores for both networks.
		symmetricPerm = (perm[0], perm[3], perm[2], perm[1])
		for permutation in [perm, symmetricPerm]:

			# Create frequencies array for permuted network
			transformed = getFrequenciesFromPermutation(originalTransformed, permutation)
			for poly in invariants:
				invariantSum = mpmath.fadd(invariant_values[poly.getPolyString()][dictionaryString], mpmath.fabs(poly.evaluate_mpm(transformed)))
				invariant_values[poly.getPolyString()][dictionaryString] = invariantSum

		# divide by the number of symmetric MSAs for each invariant, to get the mean abs.
		for poly in invariants:
			invariant_values[poly.getPolyString()][dictionaryString] = mpmath.fdiv(invariant_values[poly.getPolyString()][dictionaryString], 2)

	minScore = 1
	minPermString = ""
	for perm in permutations:
		dictionaryString = "("+ str(MSA.index[perm[0]]) + "," + str(MSA.index[perm[1]]) + "," + str(MSA.index[perm[2]]) + "," + str(MSA.index[perm[3]]) + ")"
		networkValues = list()

		for poly in invariants:
			val = invariant_values[poly.getPolyString()][dictionaryString]
			networkValues.append(val)

		score = get1Norm(networkValues)
		allScores[dictionaryString] = score
		if score < minScore:
			minScore = score
			minPermString = dictionaryString

		#Check for trees
	sortedScores = sorted(allScores.items(), key=operator.itemgetter(1))
	if sortedScores[7][1] < mpmath.fmul(multiplierForTrees1, sortedScores[0][1]) and sortedScores[8][1] > mpmath.fmul(multiplierForTrees2, sortedScores[7][1]):
		network1 = sortedScores[8][0]
		taxa = network1.strip("()").split(",")
		network2 = "(" + taxa[2] + "," + taxa[1] + "," + taxa[0] + "," + taxa[3] + ")"
		network3_1 = "(" + taxa[1] + "," + taxa[0] + "," + taxa[3] + "," + taxa[2] + ")"
		network3_2 = "(" + taxa[1] + "," + taxa[2] + "," + taxa[3] + "," + taxa[0] + ")"
		network4_1 = "(" + taxa[3] + "," + taxa[0] + "," + taxa[1] + "," + taxa[2] + ")"
		network4_2 = "(" + taxa[3] + "," + taxa[2] + "," + taxa[1] + "," + taxa[0] + ")"
		networks = [network1, network2, network3_1, network3_2, network4_1, network4_2]
		if sortedScores[9][0] in networks and sortedScores[10][0] in networks and sortedScores[11][0] in networks:
			split1 = ""
			split2 = ""
			if taxa[0] < taxa[2]:
				split1 = "(" + taxa[0] + "," + taxa[2] + ")"
			else:
				split1 = "(" + taxa[2] + "," + taxa[0] + ")"
			if taxa[1] < taxa[3]:
				split2 = "(" + taxa[1] + "," + taxa[3] + ")"
			else:
				split2 = "(" + taxa[3] + "," + taxa[1] + ")"
			if split1 < split2:
				minPermString = "(" + split1 + "," + split2 + ")"
			else:
				minPermString = "(" + split2 + "," + split1 + ")"
	
	return minPermString


if __name__ == '__main__':
	start = time.time()
	MSAFilename = ""
	scoringFunction = get1Norm
	climbingScore = False

	multiplierForTrees1 = 3
	multiplierForTrees2 = 2

	numberOfBootstraps = 100
	doPlots = False
	numProcesses = 4

	try:
		opts, args = getopt.getopt(sys.argv[1:],"ha:i:s:t:b:")
	except getopt.GetoptError:
		print("Option not recognised.")
		print("python evaluate.py -a <MSA file> -i <invariants file> -s <scoring method> -b <bootstrap replicates> -t <threads>")
		print("python evaluate.py -h for further usage instructions.")
		sys.exit(2)
	for opt, arg in opts:
		if opt == "-h":
			print("python evaluate.py -a <MSA file> -i <invariants file> -s <scoring method>")
			print("-a <MSA file>\t\t Multiple sequence alignment file.")
			print("-i <invariants file>\t\t File containing list of polynomial invariants to use in Fourier coordinates.")
			print("-s <scoring method>\t\t Method of scoring networks from invariants.")
			print("-b <sbootstrap replicates>\t Number of bootstrap replicates to perform. Default 100.")
			print("-t <threads>\t\t Number of threads (i.e. simultaneous bootstrap replicates). Default 1.")
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
		elif opt in ("-t"):
			numProcesses = int(arg)
		elif opt in ("-b"):
			numberOfBootstraps = int(arg)

	if len(MSAFilename) == 0:
		print("Error: You must provide an MSA file with -a.")
		sys.exit(2)

	try:
		originalMSA = TabularMSA.read(MSAFilename, constructor=DNA)
	except (ValueError, TypeError, skbio.io.UnrecognizedFormatError) as e:
		print(e)
		sys.exit(2)
	if not originalMSA or len(originalMSA) != 4:
		print("Error: MSA file must be a multiple sequence alignment of exactly 4 sequences.")
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

	winners = dict()
	sortedDictStrings = sorted([str(originalMSA.index[0]), str(originalMSA.index[1]), str(originalMSA.index[2]), str(originalMSA.index[3])])
	winners["((" + sortedDictStrings[0] + "," + sortedDictStrings[1] + "),(" + sortedDictStrings[2] +  "," + sortedDictStrings[3] + "))"] = 0
	winners["((" + sortedDictStrings[0] + "," + sortedDictStrings[2] + "),(" + sortedDictStrings[1] +  "," + sortedDictStrings[3] + "))"] = 0
	winners["((" + sortedDictStrings[0] + "," + sortedDictStrings[3] + "),(" + sortedDictStrings[1] +  "," + sortedDictStrings[2] + "))"] = 0
	#Permutations of the 4-cycle network: = S_4 / <(24)>
	permutations = [(0,1,2,3),(0,2,1,3),(0,1,3,2),(1,2,0,3),(1,0,2,3),(1,0,3,2),(2,1,0,3),(2,0,1,3),(2,0,3,1),(3,1,0,2),(3,0,1,2),(3,0,2,1)]
	for perm in permutations:
		dictionaryString = "("+ str(originalMSA.index[perm[0]]) + "," + str(originalMSA.index[perm[1]]) + "," + str(originalMSA.index[perm[2]]) + "," + str(originalMSA.index[perm[3]]) + ")"
		winners[dictionaryString] = 0

	numSamples =  0 
	for col in originalMSA.iter_positions(ignore_metadata=True):
		if "-" not in col:
			numSamples += 1

	print("Performing analysis on " + str(numberOfBootstraps) + " replicates...")
	results = dict()
	pool = Pool(processes=numProcesses)
	for x in range(numberOfBootstraps):
		results[x] = pool.apply_async(evaluateBootstrap, args=(originalMSA, networkPolys, numSamples, x))

	pool.close()
	pool.join()

	for x in range(numberOfBootstraps):
		winners[results[x].get()] += 1

	for network in winners.keys():
		print(network + ": " + str(winners[network] * 100.0 / numberOfBootstraps) + "%")

	end = time.time()
	print("evaluate_v2_with_bootstrap.py time: " + str(end - start))


