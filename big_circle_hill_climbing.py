import random
import math
from graphics import *

def big_circle():
	initialize_graph(12)
	initial = initialize()
	while fitness(initial) == 0:
		initial = initialize()
	next = step(initial)
	iteration = 0
	while next != None:
		iteration += 1
		print(iteration)
		initial = next
		next = step(next)
		if iteration > 100:
			break
	solution = decode(initial)
	draw_graph(solution, iteration)
	

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

#

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

# Checks if a circle is within the screen, and returns true if it is.
def in_bounds(circle):
	cen = center(circle)
	radius = circle.getRadius()
	return cen.getX() - radius >= 0 and cen.getX() + radius <= x_size and cen.getY() - radius >= 0 and cen.getY() + radius <= y_size

# Returns the center (as a Point object) of an inputted circle.
def center(circle):
	return Point(circle.getP1().getX() + circle.getRadius(), circle.getP1().getY() + circle.getRadius())

# Initializes the first, randomly selected solution. Solutions are encoding in arrays of length 3, containing 3 32-bit bitstrings each, which represent the binary
# representations of decimal values which in turn represent the x-coordinate, y-coordinate, and radius, respectively; the strings are taken modulo x_size, y_size,
# and max_radius, respectively.
def initialize():
	arr = []
	for i in range(3):
		string = ""
		for j in range(32):
			if random.random() <= 0.5:
				string += "0"
			else:
				string += "1"
		arr.append(string)
	return arr

# Decodes an array of bitstrings into the corresponding circle.
def decode(solution):
	return Circle(Point(int(solution[0], 2) % x_size, int(solution[1], 2) % y_size), int(solution[2], 2) % max_radius)

# Computes the fitness of an inputted solution, which is set to 0 if the circle represented if out of bounds or intersects initial circles, and is proportional
# to the area otherwise.
def fitness(solution):
	c = decode(solution)
	if intersects(c) or not in_bounds(c):
		return 0
	return c.getRadius()**2


# Given an input solution, computes the most optimal step to take. A step is a one-bit mutation that yields the highest, if any, increase in fitness. Returns
# None if the inputted solution is already a (local) optimum.
def step(solution):
	max_fitness = fitness(solution)
	best_step = None
	for i in range(len(solution[0])):
		new_x = solution[0][0 : i] + str((eval(solution[0][i]) + 1) % 2) + solution[0][i + 1 : len(solution[0])]
		for j in range(len(solution[1])):
			new_y = solution[1][0 : j] + str((eval(solution[1][j]) + 1) % 2) + solution[1][j + 1 : len(solution[1])]
			for k in range(len(solution[2])):
				new_radius = solution[2][0 : k] + str((eval(solution[2][k]) + 1) % 2) + solution[2][k + 1 : len(solution[0])]
				test = [new_x, new_y, new_radius]
				test_fitness = fitness(test)
				if test_fitness > max_fitness:
					max_fitness = test_fitness
					best_step = test
	return best_step

# Displays the graph, along with the inputted solution, filled in with red. Titles the graph using the given optional argument of iteration (on which the solution
# was found).
def draw_graph(solution = None, iteration = -1):
	if iteration == -1:
		name = "Big Circle"
	else:
		name = "Big Circle - Solution found on iteration " + str(iteration)
	graph = GraphWin(name, x_size, y_size)
	graph.setBackground("white")
	graph.setCoords(0, 0, x_size, y_size)
	for circle in circles:
		circle.draw(graph)
	if solution != None:
		solution.setOutline('red')
		solution.setFill('red')
		solution.draw(graph)

