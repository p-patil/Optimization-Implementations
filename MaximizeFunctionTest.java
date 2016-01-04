import java.lang.Math;

public class MaximizeFunctionTest {
	public static void main(String[] args) {
		MaximizeFunction m = new MaximizeFunction();
		double[][] constraints = {
			{-3.0, 3.0}, 
			{-3.0, 3.0}
		};
		long start = System.nanoTime();
		double[] test = m.fly(constraints);
		long finish = System.nanoTime();
		System.out.println("Time: " + getTimeElapsed(start, finish));
		System.out.println("Maximum at (x = " + test[0] + ", y = " + test[1] + ")");
		System.out.println("f(x, y) = " + m.fitness_function(test));
	}

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