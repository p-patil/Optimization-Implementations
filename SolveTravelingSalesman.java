import java.util.List;
import java.util.Arrays;
import java.util.ArrayList;
import java.util.HashMap;
import java.lang.Math;

/**
 * Class implementing genetic algorithms for use in approximating the Traveling Salesman Problem.
 */
public class SolveTravelingSalesman {

	final int population_size = 800; // How big the population is; larger populations take longer to evolve but have higher chance of producing an optimal solution per generation
	final double crossover_rate = 0.7; // Probability of a crossover occurring
	final double mutation_rate = 0.07; // Probability of mutation occurring
	final int tournament_size = 30; // Size of the tournament to use in tournament selection; bigger size means more bias against weak individuals
	final int num_generations = 1000; // Number of generations to run for
	
	int chromosome_length; // How many alleles a chromosome contains
	Organism[] population; // The population

	/**
	 * Class that represents an organism or individual in this population.
	 */
	public class Organism {
		City[] gene; // This organism's genetic code, represented by a list of cities whose order is the order in which to visit the cities

		/**
		 * Constructor.
		 * Initializes an organism with a random genome.
		 */
		public Organism() {
			if (SolveTravelingSalesman.this.cities == null) {
				throw new RuntimeException("Instance has not been initialized");
			}
			this.gene = new City[SolveTravelingSalesman.this.chromosome_length];
			List<City> temp = new ArrayList<City>(Arrays.asList(SolveTravelingSalesman.this.cities));
			int i = 0;
			while (temp.size() > 0) {
				int index = (int) (Math.random() * temp.size());
				this.gene[i] = temp.get(index);
				temp.remove(index);
				i++;
			}
		}

		/**
		 * Constructor.
		 * Initializes an organism with the specified genome.
		 */
		public Organism(City[] gene) {
			this.gene = gene;
		}
	}

	City[] cities; // List of all cities considered in this problem.

	/**
	 * Class representing a city.
	 */
	public static class City {
		double x; // x coordinate of this city
		double y; // y coordinate of this city
		int id; // The ID of this city, to keep track of cities

		/**
		 * Constructor.
		 * Initializes a city at the inputted location and with the specified ID.
		 * @param x The x-coordinate of this city.
		 * @param y The y-coordinate of this city.
		 * @param id The ID of this city.
		 */
		public City(double x, double y, int id) {
			this.x = x;
			this.y = y;
			this.id = id;
		}

		/**
		 * Whether or not the inputted list of cities contains the given city. Assumes that the list has no null cities.
		 * @param list The list of cities.
		 * @param c The city to check.
		 * @return Whether or not c is in list.
		 */
		public static boolean containsCity(List<City> list, City c) {
			if (c == null) {
				return false;
			}
			for (City d : list) {
				if (d.equals(c)) {
					return true;
				}
			}
			return false;
		}

		@Override
		/**
		 * Whether or not two cities are in fact the same.
		 * @param city The city to check this city against.
		 * @return If this city and city have the same id.
		 */
		public boolean equals(Object city) {
			return this.id == ((City) city).id;
		}

		@Override
		/**
		 * String representation of a city.
		 * @return The city's ID.
		 */
		public String toString() {
			// return "City " + this.id + " located at (" + this.x + ", " + this.y + ")";
			return Integer.toString(this.id);
		}
	}

	/**
	 * Constructor.
	 * Creates an instance of this class with the given list of cities.
	 * @param cities The list of cities to use.
	 */
	public SolveTravelingSalesman(List<double[]> cities) {
		this.chromosome_length = cities.size();
		this.cities = new City[cities.size()];
		for (int i = 0; i < cities.size(); i++) {
			this.cities[i] = new City(cities.get(i)[0], cities.get(i)[1], i);
		}
	}

	/**
	 * Basic constructor.
	 */
	public SolveTravelingSalesman() {

	}

	/**
	 * Initializes a population of organisms with random genomes.
	 * @return A list of organisms representing the new population.
	 */
	public Organism[] initializePopulation() {
		Organism[] pop = new Organism[this.population_size];
		for (int i = 0; i < pop.length; i++) {
			pop[i] = new Organism();
		}
		return pop;
	}

	/**
	 * Computes the fitness of the inputted organism.
	 * @param o The organism.
	 * @return The total distance covered by the cities specified by o's genome.
	 */
	public double fitness(Organism o) {
		double total_distance = 0.0;
		for (int i = 0; i < o.gene.length - 1; i++) {
			total_distance += Math.sqrt(Math.pow(o.gene[i].x - o.gene[i + 1].x, 2) + Math.pow(o.gene[i].y - o.gene[i + 1].y, 2));
		}
		return 1.0 / (total_distance + Math.sqrt(Math.pow(o.gene[o.gene.length - 1].x - o.gene[0].x, 2) + Math.pow(o.gene[o.gene.length - 1].y - o.gene[0].y, 2)));
	}

	/**
	 * Computes the fitness of a candidate solution.
	 * @param cities The ordered list of cities to check.
	 * @return The fitness function used is the reciprocal of the sum of Euclidean distances between consecutive cities in the inputted path.
	 */
	public static double fitness(City[] cities) {
		double total_distance = 0.0;
		for (int i = 0; i < cities.length - 1; i++) {
			total_distance += Math.sqrt(Math.pow(cities[i].x - cities[i + 1].x, 2) + Math.pow(cities[i].y - cities[i + 1].y, 2));
		}
		total_distance += Math.sqrt(Math.pow(cities[cities.length - 1].x - cities[0].x, 2) + Math.pow(cities[cities.length - 1].y - cities[0].y, 2));
		return 1.0 / total_distance;
	}

	/**
	 * Given two parents, produces a child whose gene is the result of crossing over of the two parents. Crossing over happens with probability crossover_rate; 
	 * otherwise, the more fit of the two parents is returned. Crossing over entails choosing a subset (of random size) of the first parent's genome using that 
	 * exact selection of the genome in the new genome; the remainder is filled in with the genome of the second parent in the original order that it was in, 
	 * skipping any alleles that are already contained in the chosen subset of the first parent's genome.
	 * @param o1 The first parent.
	 * @param o2 The second parent.
	 * @return The child.
	 */
	public Organism crossover(Organism o1, Organism o2) {
		if (Math.random() < this.crossover_rate) {
			City[] gene = new City[this.chromosome_length];
			int index1 = (int) (Math.random() * o1.gene.length);
			int index2 = index1;
			while (index2 == index1) {
				index2 = (int) (Math.random() * o2.gene.length);
			}
			ArrayList<City> subset = new ArrayList<>();
			for (int i = Math.min(index1, index2); i <= Math.max(index1, index2); i++) {
				subset.add(o1.gene[i]);
				gene[i] = o1.gene[i];
			}
			for (int i = 0, j = 0; i < o2.gene.length; i++) {
				if (!City.containsCity(subset, o2.gene[i])) {
					while (City.containsCity(subset, gene[j])) {
						j++;
					}
					gene[j] = o2.gene[i];
					j++;
				}
			}
			return new Organism(gene);
		} else {
			double fitness1 = fitness(o1);
			double fitness2 = fitness(o2);
			if (fitness1 >= fitness2) {
				return o1;
			} else {
				return o2;
			}
		}
	}

	/**
	 * With probability mutation_rate, mutates this organism. A mutation entails randomly choosing two locations on the organism's genome and swapping them.
	 * @param o The organism to mutate.
	 */
	public void mutate(Organism o) {
		if (Math.random() <= this.mutation_rate) {
			int index1 = (int) (Math.random() * o.gene.length);
			int index2 = index1;
			while (index2 == index1) {
				index2 = (int) (Math.random() * o.gene.length);
			}
			City temp = o.gene[index1];
			o.gene[index1] = o.gene[index2];
			o.gene[index2] = temp;
		}
	}

	/**
	 * Implements the tournament selection algorithm to semi-randomly choose an individual from the population.
	 */
	public Organism tournamentSelection() {
		Organism[] tournament = new Organism[this.tournament_size];
		for (int i = 0; i < tournament.length; i++) {
			tournament[i] = this.population[(int) (Math.random() * this.population.length)];
		}
		Organism bestOrgranism = tournament[0];
		double max_fitness = fitness(bestOrgranism);
		for (Organism o : tournament) {
			double curr_fitness = fitness(o);
			if (curr_fitness > max_fitness) {
				bestOrgranism = o;
				max_fitness = curr_fitness;
			}
		}
		return bestOrgranism;
	}

	/**
	 * Uses Roulette wheel selection (AKA fitness proportionate selection) by implementing the stochastic acceptance algorithm in order to produce a new population
	 * of individuals based on this instance's current population. Accepts a parameter that specifies whether or not to always keep the fittest organism in the
	 * current population.
	 * @param keepBestOrganism Whether or not to guarantee that the fittest organism in the population is chosen for the next population.
	 * @return The new population.
	 */
	public Organism[] nextGeneration(boolean keepBestOrganism) {
		HashMap<Organism, Double> fitnesses = new HashMap<>();
		Organism[] generation = new Organism[this.population_size];
		int offset = 0;
		if (keepBestOrganism) {
			offset = 1;
		}
		Organism bestOrgranism = this.population[0];
		double currFitness = fitness(bestOrgranism);
		for (Organism o : this.population) {
			double fitness = fitness(o);
			if (fitness > currFitness) {
				bestOrgranism = o;
				currFitness = fitness;
			}
			fitnesses.put(o, fitness);
		}
		double maxFitness = fitnesses.get(bestOrgranism);
		if (keepBestOrganism) {
			generation[0] = bestOrgranism;
		}
		Organism[] selection_population = this.population;
		// Stochastic acceptance selection:
		for (int i = offset; i < generation.length; i++) {
			double p = Math.random();
			Organism chosen = selection_population[(int) (Math.random() * selection_population.length)];
			while (p >= fitnesses.get(chosen) / maxFitness) {
				p = Math.random();
				chosen = selection_population[(int) (Math.random() * selection_population.length)];
			}
			generation[i] = chosen;
		}
		return generation;
	}

	/**
	 * Evolves an initially random population of organisms for the specified number of generations; implements the crossover and mutation genetic operators, 
	 * described above.
	 * @param displayGenerations Whether or not to display the maximum and average fitness of each generation in real-time.
	 * @return The fittest organism in the population after evolving the population.
	 */
	public Organism evolve(boolean displayGenerations) {
		this.population = initializePopulation();
		for (int i = 0; i < this.num_generations; i++) {
			if (displayGenerations) {
				System.out.println("Generation " + i);
				double max_fitness = fitness(this.population[0]);
				double total = 0.0;
				for (Organism o : this.population) {
					double curr_fitness = fitness(o);
					total += curr_fitness;
					if (curr_fitness > max_fitness) {
						max_fitness = curr_fitness;
					}
				}
				System.out.println("\tMax fitness: " + max_fitness + ". Avg fitness: " + (total / this.population.length));
			}
			this.population = nextGeneration(true);
			Organism[] mated_population = new Organism[this.population.length];
			for (int j = 0; j < mated_population.length; j++) {
				Organism father = tournamentSelection();
				Organism mother = tournamentSelection();
				Organism child = crossover(father, mother);
				mutate(child);
				mated_population[j] = child;
			}
			this.population = mated_population;
		}
		Organism bestOrgranism = this.population[0];
		double max_fitness = fitness(bestOrgranism);
		for (Organism o : this.population) {
			double curr_fitness = fitness(o);
			if (curr_fitness > max_fitness) {
				bestOrgranism = o;
				max_fitness = curr_fitness;
			}
		}
		return bestOrgranism;
	}

	/**
	 * Solves the traveling salesman problem within max_iterations generations.
	 * @param max_iterations The maximum number of generations.
	 * @return The final list of cities, in order, as the solution.
	 */
	public City[] solve(boolean displayGenerations) {
		return evolve(displayGenerations).gene;
	}
}