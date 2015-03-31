/**
	03/30/2015 

	Mixture.java: implement the unified model with two multinomial mixtures

	The first mixture is background
	The second mixture is ideal point


**/

import java.util.*;

public class Main
{
	// Configuration
	public static final int K = 5;					// Number of Latent Features
	public static int N;						// Number of Users
	public static final int MAX_ITER = 100;				// Maximum Number of Iterations 
	public static SparseMatrix<Integer> data;
	public static Map<String, Integer> dict;

	// Model Parameters
	public static double[] alpha;					// K * 1
	public static double[][] beta;					// K * N
	public static double[] pi;					// K * 1
	public static double[] p, q, b;					// N * 1

	// Variational Parameters
	public static double[] gamma;					// K * 1
	public static double[][] phi;					// K * N
	public static double[] varphi;					// N * 1

	// latent Variables
	public static double[] theta;					// K * 1
	public static double[] lambda;					// N * 1
	public static double[] z;					// K * 1


	// Initialization and Read Data 
	public static void
	init(String[] args) {
		// Init
		data = new SparseMatrix<Integer>();
		dict = new HashMap<String, Integer>();
		// Read Data
		String num = "3k", rel = args[0];
		String dictDir = "../../data/" + num + "_" + rel + "/" + rel + "_dict_" + num;
		String trainDataDir = "../../data/" + num + "_" + rel + "/" + rel + "_list_" + num + ".train";
		dict = FileParser.readVocabulary(dictDir);
		FileParser.readData(data, trainDataDir, dict);
		N = dict.size();
		// Init
		alpha = new double[K]; beta = new double[K][N]; pi = new double[K];
		p = new double[N]; q = new double[N]; b = new double[N];
		gamma = new double[K]; phi = new double[K][N]; varphi = new double[N];
		theta = new double[K]; lambda = new double[N]; z = new double[K];
	}

	// Train
	public static void
	train() {
		for (int iter = 0; iter < MAX_ITER; iter++) {
			MixtureUpdate.update();
		}

		return;
	}

	// Test
	public static void
	test() {
		// todo
	}

	// Entry
	public static void
	main(String[] args) {
		if (args.length != 3) {
			System.out.println("Usage: java Main <relation> <#blocks> <reg>");
			System.out.println("Example: java Main friend 10 0.01");
			System.exit(0);
		}

		init(args);
		train();
		test();
	}
}
