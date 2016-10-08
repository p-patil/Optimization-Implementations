public interface GeneticAlgorithm {
	int population_size = 500;
	int chromosome_length = 300;
	Organism[] population = null;
	double crossover_rate = 0.7;
	double mutation_rate = 0.01;
	public static class Organism {

	}


	public Organism[] initializePopulation();

	public double fitness(Organism o);

	public Organism[] crossover(Organism o1, Organism o2);

	public void mutate(Organism o);

	public Organism[] nextGeneration();

	public Organism evolve(int max_iterations);

}