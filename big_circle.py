import random
import math
from graphics import *

# Given a list of non-overlapping circles (defined by the position of their center, on a graph, and their radius), finds the circle with maximum area that does
# not intersect any of the existing circles. If no list of circles is inputted, a random one is generated.
def big_circle(num_generations = 150, circle_list = None):
	global circles
	if circle_list == None:
		initialize_graph(num_circles)
	else:
		circles = circle_list
	population = initialize_population(pop_size, chromosome_length)
	generation = 0
	for i in range(num_generations):
		# print("Generation " + str(generation) + ", Average fitness: " + str(sum([fitness(p) for p in population]) / len(population)))
		print("Generation " + str(generation))
		population = next_generation(population)
		generation += 1
	solution = population[0]
	max_fitness = fitness(solution)
	for organism in population:
		if fitness(organism) > max_fitness:
			solution = organism
			max_fitness = fitness(organism)
	solution = Circle(Point((int(solution[0 : 32], 2) % x_size), int(solution[32 : 64], 2) % y_size), int(solution[64 : 96], 2) % max_radius)
	print(fitness(solution))
	draw_graph(solution, generation)
	return solution

# The size of the population. Must be even.
pop_size = 300

# The probability of a crossover occurring for any two given chromosomes.
crossover_rate = 0.7

# For a given chromosome, the probability of any particular bit being inverted (ie mutated).
mutation_rate = 0.08

# The maximum length a chromosome can be.
chromosome_length = 96

# The horizontal size of the window.
x_size = 500

# The vertical size of the window.
y_size = 500

# The list of initial circles on the graph.
circles = None

# The number of circles to initialize the graph with.
num_circles = 20

# The maximum radius a circle can be.
max_radius = min(x_size, y_size) / 5

# Initializes the list of circles (ie the global variable circles). If the list already contains circles, it is erased and refilled.
def initialize_graph(n):
	global circles
	circles = []
	if n >= 5:
		num_circles = int(random.random() * (n - 5)) + 5
	else:
		num_circles = int(random.random() * n)
	for i in range(num_circles):
		radius = random.random() * max_radius
		x = radius + (random.random() * (x_size - radius))
		y = radius + (random.random() * (y_size - radius))
		c = Circle(Point(x, y), radius)
		while intersects(c):
			radius = random.random() * max_radius
			x = radius + (random.random() * (x_size - radius))
			y = radius + (random.random() * (y_size - radius))
			c = Circle(Point(x, y), radius)
		circles.append(c)

# Returns whether or not the inputted circle intersects any other circle. Assumes that the inputted circle has not yet been added to the graph.
def intersects(circle):
	if circles == None:
		return False
	x = center(circle).getX()
	y = center(circle).getY()
	for c in circles:
		if (center(c).getX() - x)**2 + (center(c).getY() - y)**2 <= (c.getRadius() + circle.getRadius())**2:
			return True
	return False

# Checks if a circle is within the screen, and returns True if it is.
def in_bounds(circle):
	cen = center(circle)
	radius = circle.getRadius()
	return cen.getX() - radius >= 0 and cen.getX() + radius <= x_size and cen.getY() - radius >= 0 and cen.getY() + radius <= y_size

# Returns the center (as a Point object) of an inputted circle.
def center(circle):
	return Point(circle.getP1().getX() + circle.getRadius(), circle.getP1().getY() + circle.getRadius())

# Displays the graph, along with the inputted solution, filled in with red. Titles the graph using the given generation (on which the solution was found).
def draw_graph(solution, generation):
	graph = GraphWin("Big Circle - Solution found on generation " + str(generation), x_size, y_size)
	graph.setBackground("white")
	graph.setCoords(0, 0, x_size, y_size)
	for circle in circles:
		circle.draw(graph)
	solution.setOutline('red')
	solution.setFill('red')
	solution.draw(graph)

# Initializes the population of bitstrings.
# TEST
def initialize_population(pop_size, chromo_length):
	population = []
	for i in range(pop_size):
		chromosome = ""
		for j in range(chromo_length):
			bit = random.random()
			if bit > 0.5:
				bit = "1"
			else:
				bit = "0"
			chromosome += bit
		population.append(chromosome)
	return population

# Returns the fitness of a given chromosome. Circles are encoded in 96-bit bitstrings as follows: the bitstring is partitioned into three substrings, each with
# length 32 bits. Each substring represents an integer in binary, and is converted to decimal. The first subtring's integer representation is used for the 
# circle's x-coordinate, which is the integer modulo x_size; similarly, the second subtring's integer representation is the y-coordinate and is used modulo
# y_size, and the last substring of 32 bits represents the radius, modulo max_radius.
# TEST
def fitness(chromosome):
	radius = int(chromosome[64 : 96], 2) % max_radius
	c = Circle(Point(int(chromosome[0 : 32], 2) % x_size, int(chromosome[32 : 64], 2) % y_size), radius)
	if intersects(c) or not in_bounds(c):
		return 0
	return radius**2

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
	if random.random() <= crossover_rate:
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

# Assigns each chromosome in the population a probability of selection proportional to its fitness.
# TEST
def roulette_wheel(population):
	fitness_scores = {}
	for chromosome in population:
		fitness_scores[chromosome] = fitness(chromosome)
	total = sum(fitness_scores[chromosome] for chromosome in fitness_scores)
	probabilities = {}
	last = 0
	for chromosome in fitness_scores:
		probabilities[(last, last + fitness_scores[chromosome] / total)] = chromosome
		last += fitness_scores[chromosome] / total
	return probabilities

# Randomly selects two chromosomes from the inputted population and performs crossover and mutation on both of them to produce two children. Uses the inputted
# probability mapping from chromosomes to probabilities to randomly select chromosomes.
# TEST
def mate(population, probabilities):
	chromosome1 = chromosome2 = ""
	prob = random.random()
	for prob_range in probabilities:
		if prob >= prob_range[0] and prob <= prob_range[1]:
			chromosome1 = probabilities[prob_range]
			break
	while chromosome2 != chromosome1:
		prob = random.random()
		for prob_range in probabilities:
			if prob >= prob_range[0] and prob <= prob_range[1]:
				chromosome2 = probabilities[prob_range]
				break
	children = positional_crossover(chromosome1, chromosome2)
	return [mutate(children[0]), mutate(children[1])]

# Performs the same roulette wheel selection process as the original mate function, but implements another methodology known as stochastic acceptance. Whereas 
# the "brute force" roulette wheel method, ie a linear search through a partitioned interval of the real number line, has linear time complexity with respect to
# the population size, that of the stochastic acceptance method is roughly constant.
# Citation: A. Lipowski and D. Lipowska, Roulette-wheel selection via stochastic acceptance, e-print: arXiv:1109.3627
def stochastic_acceptance_mate(population, max_fitness):
	chromosome1 = population[int(random.random() * len(population))]
	while random.random() > fitness(chromosome1) / max_fitness:
		chromosome1 = population[int(random.random() * len(population))]
	chromosome2 = population[int(random.random() * len(population))]
	while random.random() > fitness(chromosome2) / max_fitness:
		chromosome2 = population[int(random.random() * len(population))]
	children = uniform_crossover(chromosome1, chromosome2)
	return [mutate(children[0]), mutate(children[1])]

# Uses the inputted population to generate a new population representing the next generation. The new population is descended from the old, but takes into 
# account each chromosome's fitness, as well as introduces crossover for mating and mutations.
# TEST
def next_generation(population):
	# probabilities = roulette_wheel(population)
	max_fitness = max([fitness(p) for p in population])
	new_population = []
	while len(new_population) < len(population):
		new_population += stochastic_acceptance_mate(population, max_fitness)
	return new_population



big_circle()