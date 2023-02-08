#!/usr/bin/env python

import sys, getopt, errno, time
import numpy as np
import itertools
import mpmath
from fractions import Fraction

#---------------------------------------------------------------
class Monomial:
	def __init__(self, coefficient):
		self.powers = np.zeros(shape=(4,4,4,4), dtype=int)
		self.coefficient = coefficient
	
	def __eq__(self, other):
		if isinstance(other, Monomial):
			return self.__key() == other.__key()
		return NotImplemented

	def __ne__(self, other):
		return not __eq__(other)
	

	def __key(self):
		pows = list()
		pows.append(self.coefficient)
		for i in range(4):
			for j in range(4):
				for k in range(4):
					for l in range(4):
						pows.append(self.powers[i,j,k,l])
		return tuple(pows)

	def __hash__(self):
		return hash(self.__key())


	def getPrintString(self):
		print_string = str(self.coefficient) + "*"
		for i in range(4):
			for j in range(4):
				for k in range(4):
					for l in range(4):
						variable_string = "q(" + str(i) + str(j) + str(k) + str(l) + ")"
						if self.powers[i,j,k,l] == 1:
							print_string += variable_string
						elif self.powers[i,j,k,l] > 1:
							print_string += (variable_string + "^" + str(self.powers[i,j,k,l]))
		return print_string

	def getTotalDegree(self):
		degree = 0
		for i in range(4):
			for j in range(4):
				for k in range(4):
					for l in range(4):
						degree += self.powers[i,j,k,l]
		return degree

	def evaluate(self, values):
		value = np.longdouble(self.coefficient)

		# loop over all elements since most powers are zero and will not contribute
		# instead of np.power(values, power)
		for i in range(4):
			for j in range(4):
				for k in range(4):
					for l in range(4):
						power = self.powers[i,j,k,l]
						if power > 0:
							value *= np.power(values[i,j,k,l], power)
		return value

	def evaluate_mpm(self, values):
		value = mpmath.mpf(self.coefficient)

		# loop over all elements since most powers are zero and will not contribute
		# instead of np.power(values, power)
		for i in range(4):
			for j in range(4):
				for k in range(4):
					for l in range(4):
						power = self.powers[i,j,k,l]
						if power > 0:
							value = mpmath.fmul(value, mpmath.power(values[i,j,k,l], power))
		return value
#---------------------------------------------------------------
class Polynomial:
	def __init__(self, polynomialString=None, monomialSet=None):
		self.monomials = set()
		if monomialSet:
			if isinstance(monomialSet, set):
				for monomial in monomialSet:
					self.monomials.add(monomial)
			else:
				print("Error, must pass a set of Monomial objects")
				sys.exit(2)
		else:
			try:
				index = 0
				while index < len(polynomialString):

					# get the coefficient
					c = polynomialString[index]
					coefficient_string = ""
					while (c != 'q' and c != 'p') and index < len(polynomialString):
						coefficient_string += c
						index += 1
						if index < len(polynomialString):
							c = polynomialString[index]
						
					coefficient = 1
					if len(coefficient_string) > 0:
						if coefficient_string == '+':
							coefficient = 1
						elif coefficient_string == '-':
							coefficient = -1
						else:
							try:
								coefficient = float(Fraction(coefficient_string))
							except(ValueError) as e:
								print(e)
								print(coefficient_string)
					monomial = Monomial(coefficient)

					while c != '+' and c != '-' and index < len(polynomialString):
						if c == 'q' or c == 'p':
							var = polynomialString[index:index+11]
							i = int(str(var[3:4]))
							j = int(str(var[5:6]))
							k = int(str(var[7:8]))
							l = int(str(var[9:10]))
							index = index + 11
							power = 1
							if index < len(polynomialString) and polynomialString[index] == '^':
								index += 1
								power_string = ""
								while index < len(polynomialString) and polynomialString[index].isdigit():
									power_string += polynomialString[index]
									index += 1
								power = int(power_string)

							try:
								monomial.powers[i,j,k,l] = power
							except (IndexError) as e: 
								print(e)

							if index >= len(polynomialString):
								break
							else:
								c = polynomialString[index]
					self.monomials.add(monomial)
			except IndexError as e:
				print(e)
				print("Could not parse polynomial " + polynomialString)
				sys.exit(2)

	def __eq__(self, other):
		return self.monomials == other.monomials

	def getPolyString(self):
		print_string = ""
		for monomial in self.monomials:
				print_string += monomial.getPrintString() + "+"
		return(print_string[:-1])

	def isHomogeneous(self):
		element = next(iter(self.monomials))
		firstDegree = element.getTotalDegree()
		for monomial in self.monomials:
			if monomial.getTotalDegree() != firstDegree:
				return False
		return True

	def getTotalDegree(self):
		element = next(iter(self.monomials))
		totalDegree = element.getTotalDegree()
		for monomial in self.monomials:
			monomialTotalDegree = monomial.getTotalDegree()
			if monomialTotalDegree > totalDegree:
				totalDegree = monomialTotalDegree
		return totalDegree

	def evaluate(self, values):
		totalSum = np.longdouble(0)
		for monomial in self.monomials:
			totalSum += monomial.evaluate(values)
		return totalSum

	def evaluate_mpm(self, values):
		totalSum = mpmath.mpf(0)
		for monomial in self.monomials:
			totalSum = mpmath.fadd(totalSum, monomial.evaluate_mpm(values))
		return totalSum
#---------------------------------------------------------------

def Permute(poly, perm):
	permutedMonomials = set()
	for monomial in poly.monomials:
		permuted = Monomial(monomial.coefficient)
		for i in range(4):
			for j in range(4):
				for k in range(4):
					for l in range(4):
						permuted.powers[perm[i],perm[j],perm[k],perm[l]] = monomial.powers[i,j,k,l]
		permutedMonomials.add(permuted)	
	return Polynomial(monomialSet=permutedMonomials)

def Negate(poly):
	negatedMonomials = set()
	for monomial in poly.monomials:
		negated = Monomial(-monomial.coefficient)
		for i in range(4):
			for j in range(4):
				for k in range(4):
					for l in range(4):
						negated.powers[i,j,k,l] = monomial.powers[i,j,k,l]
		negatedMonomials.add(negated)	
	return Polynomial(monomialSet=negatedMonomials)

def add(monomial1, monomial2):
	newMonomial = Monomial(coefficient=monomial1.coefficient + monomial2.coefficient)
	for i in range(4):
			for j in range(4):
				for k in range(4):
					for l in range(4):
						newMonomial.powers[i,j,k,l] = monomial1.powers[i,j,k,l]
						if monomial1.powers[i,j,k,l] != monomial2.powers[i,j,k,l]:
							return Polynomial(monomialSet = {monomial1, monomial2})

	return Polynomial(monomialSet={newMonomial})




