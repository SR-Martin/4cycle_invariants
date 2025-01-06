#!/usr/bin/env python

import sys, getopt, errno, time
import numpy as np
import mpmath
from skbio import DNA, TabularMSA
import utils.polynomials
from utils.polynomials import Polynomial
import operator

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

MSAFilename = ""
scoringFunction = get1Norm
climbingScore = False

multiplierForTrees1 = 3
multiplierForTrees2 = 2

fourierValues = ""

try:
	opts, args = getopt.getopt(sys.argv[1:],"ha:i:s:f:")
except getopt.GetoptError:
	print("Option not recognised.")
	print("python evaluate.py -a <MSA file> -i <invariants file> -s <scoring method>")
	print("python evaluate.py -h for further usage instructions.")
	sys.exit(2)
for opt, arg in opts:
	if opt == "-h":
		print("python evaluate.py -a <MSA file> -i <invariants file> -s <scoring method>")
		print("-a <MSA file>\t\t Multiple sequence alignment file.")
		print("-i <invariants file>\t\t File containing list of polynomial invariants to use in Fourier coordinates.")
		print("-s <scoring method>\t\t Method of scoring networks from invariants.")
		sys.exit()
	elif opt in ("-a"):
		MSAFilename = arg
	elif opt in ("-i"):
		NetworkGBFilename = arg
	elif opt in ("-f"):
		fourierValues = arg
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

if len(MSAFilename) == 0 and len(fourierValues) == 0:
	print("Error: You must provide an MSA file with -a.")
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


networkNorms = dict()

#print("Calculating scores...")
invariant_values = dict()
for poly in networkPolys:
	invariant_values[poly.getPolyString()] = dict()

if len(MSAFilename) > 0:
	try:
		originalMSA = TabularMSA.read(MSAFilename, constructor=DNA)
	except (ValueError, TypeError, skbio.io.UnrecognizedFormatError) as e:
		print(e)
		sys.exit(2)
	if not originalMSA or len(originalMSA) != 4:
		print("Error: MSA file must be a multiple sequence alignment of exactly 4 sequences.")
		sys.exit(2)

	# Create frequencies array from MSA
	originalFrequencies = np.zeros(shape=(4,4,4,4), dtype=np.longdouble)
	count = 0
	for col in originalMSA.iter_positions(ignore_metadata=True):
		if "-" not in col:
			originalFrequencies[Nucl[str(col[0])], Nucl[str(col[1])], Nucl[str(col[2])], Nucl[str(col[3])]] += 1
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
else:
	# Read in Fourier values from string
	stringValues = fourierValues.split(",")
	values = [0] * len(stringValues)
	for i in range(len(stringValues)):
		values[i] = mpmath.mpmathify(stringValues[i])
	originalTransformed = np.full((4,4,4,4), 0, dtype=np.longdouble)

	originalTransformed[0,0,0,0] = values[0]

	originalTransformed[0,0,1,1] = values[1]
	originalTransformed[0,0,2,2] = values[1]
	originalTransformed[0,0,3,3] = values[1]

	originalTransformed[0,1,0,1] = values[2]
	originalTransformed[0,2,0,2] = values[2]
	originalTransformed[0,3,0,3] = values[2]

	originalTransformed[0,1,1,0] = values[3]
	originalTransformed[0,2,2,0] = values[3]
	originalTransformed[0,3,3,0] = values[3]

	originalTransformed[0,1,2,3] = values[4]
	originalTransformed[0,1,3,2] = values[4]
	originalTransformed[0,2,1,3] = values[4]
	originalTransformed[0,2,3,1] = values[4]
	originalTransformed[0,3,1,2] = values[4]
	originalTransformed[0,3,2,1] = values[4]

	originalTransformed[1,0,0,1] = values[5]
	originalTransformed[2,0,0,2] = values[5]
	originalTransformed[3,0,0,3] = values[5]

	originalTransformed[1,0,1,0] = values[6]
	originalTransformed[2,0,2,0] = values[6]
	originalTransformed[3,0,3,0] = values[6]

	originalTransformed[1,0,2,3] = values[7]
	originalTransformed[1,0,3,2] = values[7]
	originalTransformed[2,0,1,3] = values[7]
	originalTransformed[2,0,3,1] = values[7]
	originalTransformed[3,0,1,2] = values[7]
	originalTransformed[3,0,2,1] = values[7]

	originalTransformed[1,1,0,0] = values[8]
	originalTransformed[2,2,0,0] = values[8]
	originalTransformed[3,3,0,0] = values[8]

	originalTransformed[1,1,1,1] = values[9]
	originalTransformed[2,2,2,2] = values[9]
	originalTransformed[3,3,3,3] = values[9]

	originalTransformed[1,2,0,3] = values[10]
	originalTransformed[1,3,0,2] = values[10]
	originalTransformed[2,1,0,3] = values[10]
	originalTransformed[2,3,0,1] = values[10]
	originalTransformed[3,1,0,2] = values[10]
	originalTransformed[3,2,0,1] = values[10]

	originalTransformed[1,2,1,2] = values[11]
	originalTransformed[1,3,1,3] = values[11]
	originalTransformed[2,1,2,1] = values[11]
	originalTransformed[2,3,2,3] = values[11]
	originalTransformed[3,1,3,1] = values[11]
	originalTransformed[3,2,3,2] = values[11]

	originalTransformed[1,2,3,0] = values[12]
	originalTransformed[1,3,2,0] = values[12]
	originalTransformed[2,1,3,0] = values[12]
	originalTransformed[2,3,1,0] = values[12]
	originalTransformed[3,1,2,0] = values[12]
	originalTransformed[3,2,1,0] = values[12]

	originalTransformed[1,1,2,2] = values[13]
	originalTransformed[1,1,3,3] = values[13]
	originalTransformed[2,2,1,1] = values[13]
	originalTransformed[2,2,3,3] = values[13]
	originalTransformed[3,3,1,1] = values[13]
	originalTransformed[3,3,2,2] = values[13]

	originalTransformed[1,2,2,1] = values[14]
	originalTransformed[1,3,3,1] = values[14]
	originalTransformed[2,1,1,2] = values[14]
	originalTransformed[2,3,3,2] = values[14]
	originalTransformed[3,1,1,3] = values[14]
	originalTransformed[3,2,2,3] = values[14]

	count = 0


#Permutations of the 4-cycle network: = S_4 / <(24)>
permutations = [(0,1,2,3),(0,2,1,3),(0,1,3,2),(1,2,0,3),(1,0,2,3),(1,0,3,2),(2,1,0,3),(2,0,1,3),(2,0,3,1),(3,1,0,2),(3,0,1,2),(3,0,2,1)]
for perm in permutations:
	dictionaryString = "("+ str(perm[0]) + "," + str(perm[1]) + "," + str(perm[2]) + "," + str(perm[3]) + ")"
	
	for poly in networkPolys:
			invariant_values[poly.getPolyString()][dictionaryString] = 0

	# The 2-cycle (24) gives the same network, so we would like the result to be the same in both cases
	# We achieve this by calculating the frequency scores for both networks.
	symmetricPerm = (perm[0], perm[3], perm[2], perm[1])
	for permutation in [perm, symmetricPerm]:

		# Create frequencies array for permuted network
		transformed = getFrequenciesFromPermutation(originalTransformed, permutation)
		for poly in networkPolys:
			invariantSum = mpmath.fadd(invariant_values[poly.getPolyString()][dictionaryString], mpmath.fabs(poly.evaluate_mpm(transformed)))
			invariant_values[poly.getPolyString()][dictionaryString] = invariantSum

	# divide by the number of symmetric MSAs for each invariant, to get the mean abs.
	for poly in networkPolys:
		invariant_values[poly.getPolyString()][dictionaryString] = mpmath.fdiv(invariant_values[poly.getPolyString()][dictionaryString], 2)

print("Evaluating alignment with total length (without gaps): " + str(count))	
for perm in permutations:
	dictionaryString = "("+ str(perm[0]) + "," + str(perm[1]) + "," + str(perm[2]) + "," + str(perm[3]) + ")"
	if len(MSAFilename) > 0:
		labelString = "("+ str(originalMSA.index[perm[0]]) + "," + str(originalMSA.index[perm[1]]) + "," + str(originalMSA.index[perm[2]]) + "," + str(originalMSA.index[perm[3]]) + ")"
	else:
		labelString = "(Taxon"+ str(perm[0]+1) + ",Taxon" + str(perm[1]+1) + ",Taxon" + str(perm[2]+1) + ",Taxon" + str(perm[3]+1) + ")"
	networkValues = list()

	for poly in networkPolys:
		val = invariant_values[poly.getPolyString()][dictionaryString]
		networkValues.append(val)

	if not climbingScore:
		score = scoringFunction(networkValues)
	else:
		score = getClimbingScore(dictionaryString)
	networkNorms[labelString] = score

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
				treeString = "(" + split1 + "," + split2 + ")"
			else:
				treeString = "(" + split2 + "," + split1 + ")"
			print("Tree-like evolution detected with tree " + treeString + ".")

end = time.time()
print("evaluate_v3.py time: " + str(end - start))


