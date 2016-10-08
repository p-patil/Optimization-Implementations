import java.util.List;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.HashMap;
import java.lang.Math;

/**
 * Tester class for running the SolveTravelingSalesman class.
 */
public class SolveTravelingSalesmanTest {
	public static void main(String[] args) {
		test(Integer.parseInt(args[0]));
	}

	/**
	 * Randomly generates a list of cities and uses the SolveTravelingSalesman class to approximate a solution; compares the approximation to
	 * the correct solution, as well as time taken.
	 */
	public static void test(int n) {
		ArrayList<double[]> cities = new ArrayList<>();
		for (int i = 0; i < n; i++) {
			double[] coordinates = {Math.random() * 100, Math.random() * 100};
			cities.add(coordinates);
		}

		SolveTravelingSalesman solver = new SolveTravelingSalesman(cities);

		long start = System.nanoTime();
		SolveTravelingSalesman.City[] actual = solver.solve(false);
		long finish = System.nanoTime();
		String time = getTimeElapsed(start, finish);
		System.out.print("Genetic algorithm took " + time);

		start = System.nanoTime();
		SolveTravelingSalesman.City[] expected = bruteForceSolve(solver.cities);
		finish = System.nanoTime();
		time = getTimeElapsed(start, finish);
		System.out.println(", brute force took " + time);

		double expectedDistance =  1.0 / SolveTravelingSalesman.fitness(expected);
		double gotDistance = 1.0 / SolveTravelingSalesman.fitness(actual);
		System.out.println("Expected total distance was " + expectedDistance + ", got " + gotDistance);
		System.out.println("Error: " + (100 * Math.abs(expectedDistance - gotDistance) / expectedDistance) + "%");
	}

	/**
	 * Implementation of the Held-Karp dynamic programming algorithm to solve TSP in O(n^2 * 2^n) time, for n = number of cities.
	 */
	// public static SolveTravelingSalesman.City[] heldKarp(SolveTravelingSalesman.City[] a) {
	// 	// Class used to index the HashMap used in the dynamic programming.
	// 	public static class Index {
	// 		HashSet<Integer> subset;
	// 		int j;

	// 		public Index(HashSet<Index> subset, int j) {
	// 			this.subset = subset;
	// 			this.j = j;
	// 		}

	// 		@Override
	// 		public boolean equals(Object obj) {
	// 			return this.subset.equals(((Index) obj).subset) && (this.subset.j == ((Index) obj).j);
	// 		}

	// 		@Override
	// 		// The hash code is simply the hash code of subset spliced with j at the end.
	// 		public int hashCode() {
	// 			int count = 1;
	// 			for (int n = this.j; n > 0; n /= 10, count *= 10); // count now equals 10^(number of digits in n)
	// 			return this.subset.hashCode() * count + j;
	// 		}
	// 	}

	// 	// Initialization for dynamic programming.
	// 	HashMap<Index, Integer> c = new HashMap<>();
	// 	HashSet<Integer> temp = new HashSet<>();
	// 	temp.add(1);
	// 	c.put(new Index(temp, 1), 0);

	// 	for (int s = 2; s <= a.length; s++) {
			
	// 	}
	// }

	/**
	 * Solves TSP by exhaustively searching through all possible permutations. Runs in O(n!) time, for n = number of cities. Used for testing 
	 * purposes. Uses the quickselect algorithm to iteratively generate permutations of a, and checks them on each iteration.
	 * @param a An instance of the TSP problem (a list of city locations).
	 * @return The optimal solution, as an ordered list of cities to visit.
	 */
	public static SolveTravelingSalesman.City[] bruteForceSolve(SolveTravelingSalesman.City[] a) {
		SolveTravelingSalesman.City[] best = a;
		double max_fitness = SolveTravelingSalesman.fitness(a);
		int[] p = new int[a.length];
		int i = 1;
		int j;
		while (i < a.length) {
			if (p[i] < i) {
				if (i % 2 == 1) {
					j = p[i];
				} else {
					j = 0;
				}
				SolveTravelingSalesman.City temp = a[i];
				a[i] = a[j];
				a[j] = temp;
				double curr_fitness = SolveTravelingSalesman.fitness(a);
				if (curr_fitness > max_fitness) {
					max_fitness = curr_fitness;
					SolveTravelingSalesman.City[] b = new SolveTravelingSalesman.City[a.length];
					for (int k = 0; k < b.length; k++) {
						b[k] = a[k];
					}
					best = b;
				}
				p[i]++;
				i = 1;
			} else {
				p[i] = 0;
				i++;
			}
		}
		return best;
	}

	/**
	 * Prints the time taken between start and finish in a human readable format for all times under a year.
	 * @param start The start time.
	 * @param finish The end time.
	 */
	private static String getTimeElapsed(long start, long finish) {
		double seconds = (finish - start) * 0.000000001;
		String time = seconds + " seconds";
		if (seconds > 60) {
			int minutes = (int) (seconds / 60);
			seconds %= 60;
			time = minutes + " minutes " + seconds + " seconds"; 
			if (minutes == 1) {
				time = minutes + " minute " + seconds + " seconds";
			}
			if (minutes > 60) {
				int hours = (int) (minutes / 60);
				minutes %= 60;
				time = hours + " hours " + minutes + " minutes " + seconds + " seconds";
				if (hours == 1) {
					time = hours + " hour " + minutes + " minutes " + seconds + " seconds";					
				}
				if (hours > 24) {
					int days = (int) (hours / 24);
					hours %= 24;
					time = days + " days " + hours + " hours " + minutes + " minutes " + seconds + " seconds ";
					if (days == 1) {
						time = days + " day " + hours + " hours " + minutes + " minutes " + seconds + " seconds ";
					}
				}
			}
		}
		return time;
	}
}