import random

# Uses genetic algorithms to return an arithmetic expression containing numbers, 0 through 9, and operators, +, -, *, and /, which evaluate to n.
def operator_sequence(n):
	assert pop_size % 2 == 0, "Population size must be even."
	assert max_chromosome_length % 4 == 0, "Maximum chromosome length must be a multiple of 4, because chromosomes are made up of 4-bit chunks."
	# This is the initial population of chromosomes, and will be updated every "generation", ie iteration of the algorithm.
	population = initialize_population(pop_size, max_chromosome_length)
	generation = 0
	while type(population) != str:
		# print("Generation " + str(generation) + ", Average fitness in population: " + str(sum([fitness(p, n) for p in population]) / len(population)))
		# print("Generation " + str(generation) + ", Maximum fitness in population: " + str(max([fitness(p, n) for p in population])))
		population = next_generation(population, n)
		generation += 1
	print("Solution found on generation " + str(generation) + ": " + population)
	return population

# HELPER FUNCTIONS AND VARIABLES:

# This mapping allows for solutions to be converted (bijectively) to bitstrings.
encoding = {"0000": 0, 
			"0001": 1, 
			"0010": 2, 
			"0011": 3, 
			"0100": 4, 
			"0101": 5, 
			"0110": 6, 
			"0111": 7, 
			"1000": 8, 
			"1001": 9, 
			"1010": "+", 
			"1011": "-", 
			"1100": "*", 
			"1101": "/"}
# This is the size of the population of chromosomes. Must be even.
pop_size = 500
# The maximum length of a chromosome. Must be divisible by 4.
max_chromosome_length = 300
# This is the probability that two chromosomes will crossover.
crossover_rate = 0.7
# This is the probability that a given "gene", ie a bit, in a chromosome will be inverted.
mutation_rate = 0.1

# Initializes the population by returning a list of size random bitstrings, each with maximum length of chromo_length.
def initialize_population(size, chromo_length):
	ret = []
	for i in range(size):
		chromosome = ""
		for j in range(int((random.random() * chromo_length / 4) + 0.5) + 1):
			chunk = ""
			for k in range(4):
				bit = random.random()
				if bit < 0.5:
					chunk += "0"
				else:
					chunk += "1"
			chromosome += chunk
		ret.append(chromosome)
	return ret

# Converts a bitstring into an arithmetic expression.
def decode(bitstring):
	string = ""
	last = ""
	i = 0
	while i < len(bitstring):
		if bitstring[i : i + 4] == "1110" or bitstring[i : i + 4] == "1111":
			pass
		elif type(last) == int and type(encoding[bitstring[i : i + 4]]) == int:
			pass
		elif type(last) == str and type(encoding[bitstring[i : i + 4]]) == str:
			pass
		else:
			last = encoding[bitstring[i : i + 4]]
			string += str(last) + " "
		i += 4
	if type(last) == str:
		string = string[0 : len(string) - 2]
	return string

# Evaluates an arithmetic expression sequentially from left to right, as opposed to the normal PEMDAS order of operations.
def evaluate(operator_string):
	try:
		arr = operator_string.split()
		if len(arr) == 0:
			return 0
		elif len(arr) == 1:
			return eval(operator_string)
		ret = eval(arr[0])
		ret = eval(arr[0] + arr[1] + arr[2])
		i = 3
		while i < len(arr):
			ret = eval(str(ret) + arr[i] + arr[i + 1])
			i += 2
		return ret
	except ZeroDivisionError:
		return float("inf")

# Fitness function that maps a chromosome (a string with only 1's and 0's, whose length must be divisible by 4) to a fitness based on how close the string 
# is to the correct solution. Returns 0 if the chromosome is a perfect match.
def fitness(chromosome, n):
	value = n - evaluate(decode(chromosome))
	if value == 0:
		return float('inf')
	return 1 / abs(value)

# This generates mapping from a chromosome in the inputted population to the fitness scores, or returns a single chromosome (ie a bitstring) if it's a perfect 
# solution. 
def fitness_scores(pop, n):
	ret = {}
	for chromosome in pop:
		value = fitness(chromosome, n)
		if value == float('inf'):
			return decode(chromosome)
		ret[chromosome] = value
	return ret

# Performs a crossover on two chromosomes, and returns the resulting two chromosomes in a tuple. A crossover consists of randomly choosing a position on the two
# chromosomes and swapping the bits after that point.
def positional_crossover(chromosome1, chromosome2):
	new_chromosome1, new_chromosome2 = chromosome1, chromosome2
	if (random.random() <= crossover_rate):
		position = int(random.random() * min(len(chromosome1), len(chromosome2)))
		new_chromosome1 = chromosome1[0 : position] + chromosome2[position : len(chromosome2)]
		new_chromosome2 = chromosome2[0 : position] + chromosome1[position : len(chromosome1)]
	return (new_chromosome1, new_chromosome2)

# Alternate mechanism for crossover that consists of swapping corresponding bits between two chromosomes at a 50% chance of swap per bit. This is known as uniform
# crossover and is more similar to the biological analog.
def uniform_crossover(chromosome1, chromosome2):
	new_chromosome1, new_chromosome2 = chromosome1, chromosome2
	if (random.random() <= crossover_rate):
		for i in range(min(len(chromosome1), len(chromosome2))):
			if random.random() <= 0.5:
				new_chromosome1 = new_chromosome1[0 : i] + chromosome2[i] + new_chromosome1[i + 1 : len(new_chromosome1)]
				new_chromosome2 = new_chromosome2[0 : i] + chromosome1[i] + new_chromosome2[i + 1 : len(new_chromosome2)]
	return (new_chromosome1, new_chromosome2)

# Mutates a chromosome.
def mutate(chromosome):
	new_chromosome = chromosome
	for i in range(len(chromosome)):
		if random.random() <= mutation_rate:
			new_chromosome = new_chromosome[0 : i] + str((eval(chromosome[i]) + 1) % 2) + new_chromosome[i + 1 : len(new_chromosome)]
	return new_chromosome

# This is what crossover mechanism to use.
crossover = uniform_crossover

# A probabilistic way of creating a new generation from the current population, using the roulette wheel method in which the chance of two chromosomes being
# selected for reproduction into the new population is proportional to their fitness. If the current population contains a bitstring that is a full
# solution, then that bitstring's representation as an arithmetic sequence (stored in a String) is immediately returned instead. 
def next_generation(pop, n):
	scores = fitness_scores(pop, n)
	if type(scores) == str:
		return scores
	probabilities = {}
	total = sum([scores[chromosome] for chromosome in scores])
	last = 0
	for chromosome in scores:
		probabilities[(last, last + scores[chromosome] / total)] = chromosome
		last += scores[chromosome] / total
	new_population = []
	while len(new_population) < len(pop):
		prob = random.random()
		chromosome1 = chromosome2 = ""
		for prob_range in probabilities:
			if prob >= prob_range[0] and prob <= prob_range[1]:
				chromosome1 = probabilities[prob_range]
				break
		while (chromosome2 != chromosome1):
			prob = random.random()
			for prob_range in probabilities:
				if prob >= prob_range[0] and prob <= prob_range[1]:
					chromosome2 = probabilities[prob_range]
					break
		new_chromosomes = crossover(chromosome1, chromosome2)
		chromosome1, chromosome2 = mutate(new_chromosomes[0]), mutate(new_chromosomes[1])
		new_population += [chromosome1, chromosome2]
	return new_population[0 : len(pop)]
