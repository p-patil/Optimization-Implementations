import java.lang.Math;

/**
 * Maximizes a given function (from R^n to R^m) over given constraints. The function to be maximized should be implemented in the 
 * "fitness_function" method. By default, the function f(x, y) = sin(x + y) * cos(xy) is used, with global optimum at (~1.5708, 0).
 */
public class MaximizeFunction {
	int swarm_size = 100; // Size of the swarm
	double inertia = 0.7; // Weighting of the component of the velocity towards the old velocity
	double cognitive_coefficient = 2.0; // Weighting of the component of the velocity towards a particle's local best
	double social_coefficient = 2.0; // Weighting of the component of the velocity towards the global best
	double velocity_clamping_coefficient = 0.55; // Coefficient for velocity clamping of the particle
	int max_iterations = 1000; // How many discrete time steps to run the simulation for when searching for the solution
	Particle[] swarm; // List of all particles in the swarm
	double[] global_best; // The best position in the swarm thus far

	/**
	 * Basic constructor.
	 */
	public MaximizeFunction() {
		swarm = new Particle[swarm_size];
		global_best = new double[2];
	}

	/**
	 * Class representing the particles in the swarm.
	 */
	public class Particle {
		double[] position; // Current position vector in the search space
		double[] velocity; // Current velocity vector in the search space
		double[] personal_best; // Position vector with the highest fitness that this particle has achieved
		double[] velocity_clamping; // Vector that imposes limitations on the velocity.
		double[][] constraints;

		/**
		 * Constructor.
		 * Constructs a particle with a random position and velocity in a search space of the inputted dimensions subject to the inputted constraints.
		 * @param dimensions The number of dimensions in the search space.
		 * @param constraints An array of size dimensions filled with arrays of size 2 representing the region in the search space to optimize the function over.
		 */
		public Particle(int dimensions, double[][] constraints) {
			this.constraints = constraints;
			this.position = new double[dimensions];
			this.velocity = new double[dimensions];
			this.personal_best = new double[dimensions];
			for (int i = 0; i < position.length; i++) {
				this.position[i] = constraints[i][0] + ( Math.random() * (constraints[i][1] - constraints[i][0]));
				this.velocity[i] = (Math.random() * 2.0) - 1.0;
			}
			this.personal_best = this.position;
			this.velocity_clamping = new double[constraints.length];
			for (int i = 0; i < this.velocity_clamping.length; i++) {
				this.velocity_clamping[i] = 0.5 * (constraints[i][1] - constraints[i][0]);
			}
		}

		/**
		 * Constructor.
		 * Constructs a particle with a random position and velocity in a search space of the inputted dimensions.
		 * @param dimensions The number of dimensions in the search space.
		 */
		public Particle(int dimensions) {
			this.constraints = null;
			this.position = new double[dimensions];
			this.velocity = new double[dimensions];
			this.personal_best = new double[dimensions];
			for (int i = 0; i < position.length; i++) {
				this.position[i] = Math.random();
				this.velocity[i] = (Math.random() * 2.0) - 1.0;
			}
			this.personal_best = this.position;
			this.velocity_clamping = new double[dimensions];
			for (int i = 0; i < this.velocity_clamping.length; i++) {
				this.velocity_clamping[i] = Double.POSITIVE_INFINITY;
			}
		}

		/**
		 * Updates the position and velocity vectors of this particle.
		 * Update equations: v(t + 1) = inertia * v(t) + cognitive_coefficient * (local_best - x(t)) + social_coefficient * (global_best - x(t))
		 *					 x(t + 1) = x(t) + v(t + 1)
		 * If constraints are imposed and the updated position vector exceeds the constraints, it is reset to the maximum value it may attain, ie the constraint
		 * itself. Velocities are bounded by the velocity clamping coefficient, unless there are no constraints (in which case the velocity may be unbounded).
		 */
		public void update() {
			for (int i = 0; i < this.velocity.length; i++) {
				double inertial_component = MaximizeFunction.this.inertia * this.velocity[i];
				double cognitive_component = Math.random() * MaximizeFunction.this.cognitive_coefficient * (this.personal_best[i] - this.position[i]);
				double social_component = Math.random() * MaximizeFunction.this.social_coefficient * (MaximizeFunction.this.global_best[i] - this.position[i]);
				this.velocity[i] = inertial_component + cognitive_component + social_component;
				if (this.velocity[i] > MaximizeFunction.this.velocity_clamping_coefficient * this.velocity_clamping[i]) {
					this.velocity[i] = MaximizeFunction.this.velocity_clamping_coefficient * this.velocity_clamping[i];
				} else if (this.velocity[i] < -(MaximizeFunction.this.velocity_clamping_coefficient * this.velocity_clamping[i])) {
					this.velocity[i] = -(MaximizeFunction.this.velocity_clamping_coefficient * this.velocity_clamping[i]);
				}
				this.position[i] += this.velocity[i];
			}
			if (MaximizeFunction.this.fitness_function(this.position) > MaximizeFunction.this.fitness_function(this.personal_best)) {
				this.personal_best = this.position;
			}
			if (this.constraints != null) {
				for (int i = 0; i < this.position.length; i++) {
					if (this.position[i] > this.constraints[i][1]) {
						this.position[i] = this.constraints[i][1];
					} else if (this.position[i] < this.constraints[i][0]) {
						this.position[i] = this.constraints[i][0];
					}
				}
			}
		}
	}

	/**
	 * Initializes the swarm.
	 * @param constraints Uses the given constraints when initializing the swarm.
	 */
	public void initializeSwarm(double[][] constraints) {
		for (int i = 0; i < this.swarm.length; i++) {
			this.swarm[i] = new Particle(2, constraints);
		}
		this.global_best = swarm[0].position;
		double max_fitness = fitness_function(swarm[0].position);
		for (Particle p : swarm) {
			double curr_fitness = fitness_function(p.position);
			if (curr_fitness > max_fitness) {
				max_fitness = curr_fitness;
				this.global_best = p.position;
			}
		}
	}

	/**
	 * Initializes the swarm.
	 */
	public void initializeSwarm() {
		for (int i = 0; i < this.swarm.length; i++) {
			this.swarm[i] = new Particle(2);
		}
		this.global_best = swarm[0].position;
		double max_fitness = fitness_function(swarm[0].position);
		for (Particle p : swarm) {
			double curr_fitness = fitness_function(p.position);
			if (curr_fitness > max_fitness) {
				max_fitness = curr_fitness;
				this.global_best = p.position;
			}
		}
	}

	/**
	 * Computes the fitness at a given position vector.
	 * @param position The position vector to evaluate at.
	 * @return The fitness.
	 */
	public double fitness_function(double[] position) {
		// return Math.sin(position[0] + position[1]) * Math.cos(position[0] * position[1]);
		// return Math.log(position[0]) + Math.log(position[1]) + 2 * Math.log(0.5 * Math.sqrt(1 - position[0] - position[1]));
		return (Math.exp(-((position[0] * position[0]) + (position[1] * position[1]))) + (Math.sqrt(5) * Math.pow(Math.sin(position[0] * position[0] * position[0]), 2)) + (2 * Math.pow(Math.cos((2 * position[0]) + (3 * position[1])), 2))) / (1 + (position[0] * position[0]) + (position[1] * position[1]));
	}

	/**
	 * Computes the next iteration in the swarm.
	 */
	public void updateParticles() {
		for (Particle p : this.swarm) {
			p.update();
		}
		double[] curr_best = swarm[0].position;
		double max_fitness = fitness_function(swarm[0].position);
		for (Particle p : swarm) {
			double curr_fitness = fitness_function(p.position);
			if (curr_fitness > max_fitness) {
				max_fitness = curr_fitness;
				curr_best = p.position;
			}
		}
		if (max_fitness > fitness_function(this.global_best)) {
			this.global_best = curr_best;
		}
	}

	/**
	 * Runs the swarm for the max_iterations number of iterations of discrete time steps and returns the best position vector found, subject to constraints.
	 * @param constraints The k-cell to optimize over.
	 * @return The point at which the function is maximized.
	 */
	public double[] fly(double[][] constraints) {
		initializeSwarm(constraints);
		for (int i = 0; i < this.max_iterations; i++) {
			updateParticles();
		}
		return this.global_best;
	}

	/**
	 * Runs the swarm for the max_iterations number of iterations of discrete time steps and returns the best position vector found.
	 * @return The point at which the function is maximized.
	 */
	public double[] fly() {
		initializeSwarm();
		for (int i = 0; i < this.max_iterations; i++) {
			updateParticles();
		}
		return this.global_best;
	}
}