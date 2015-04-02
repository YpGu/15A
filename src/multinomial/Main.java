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
	public static final int MAX_ITER = 10;				// Maximum Number of Iterations 
	public static SparseMatrix<Integer> trainData, testData;
	public static Map<String, Integer> dict;
	public static Map<Integer, String> invDict;

	// Model Parameters
	public static double[] alpha;					// K * 1
	public static double[][] beta;					// K * N
	public static double[] pi;					// N * 1
	public static double[] p, q, b;					// N * 1

	// Variational Parameters
	public static double[] gamma;					// K * 1
	public static double[][] phi;					// N * K
	public static double[] varphi;					// N * 1

	// Initialization and Read Data 
	public static void
	init(String[] args) {
		// Init
		trainData = new SparseMatrix<Integer>();
		dict = new HashMap<String, Integer>();
		invDict = new HashMap<Integer, String>();
		Random rand = new Random(0);
		// Read Data
		String num = "3k", rel = args[0];
		String dictDir = "../../data/" + num + "_" + rel + "/" + rel + "_dict_" + num;
		String trainDataDir = "../../data/" + num + "_" + rel + "/" + rel + "_list_" + num + ".train";
		dict = FileParser.readVocabulary(dictDir);
		invDict = FileParser.readInverseVocabulary(dictDir);
		FileParser.readData(trainData, trainDataDir, dict);
		N = dict.size();
		// Initialize Parameters 
		alpha = new double[K]; beta = new double[K][N]; pi = new double[N];
		p = new double[N]; q = new double[N]; b = new double[N];
		gamma = new double[K]; phi = new double[N][K]; varphi = new double[N];
		for (int k = 0; k < K; k++) alpha[k] = 2;
		for (int k = 0; k < K; k++) for (int j = 0; j < N; j++) beta[k][j] = 1/(double)N;
		for (int i = 0; i < N; i++) {
			pi[i] = 0.4 + 0.2 * rand.nextDouble();
			varphi[i] = pi[i];
		}
		for (int i = 0; i < N; i++) {
			double pqRange = 5;
			p[i] = -0.5 * pqRange + pqRange * rand.nextDouble();
			q[i] = -0.5 * pqRange + pqRange * rand.nextDouble();
			b[i] = -0.5 * pqRange + pqRange * rand.nextDouble();
		}
		for (int k = 0; k < k; k++) {
			gamma[k] = alpha[k] + (double)N/K;
		}
	}

	// Train
	public static void
	train() {
		for (int iter = 0; iter < MAX_ITER; iter++) {
			System.out.println("----- Iteration " + iter + " -----");
			MixtureUpdate.update(trainData, alpha, beta, pi, p, q, b, gamma, phi, varphi);
		}
		FileParser.output("./p", p, invDict);
		FileParser.output("./q", q, invDict);
		FileParser.output("./b", b, invDict);

		return;
	}

	// Test
	public static void
	test() {
		// todo
		Evaluation.auroc(trainData, alpha, beta, pi, p, q, b, gamma, phi, varphi, 1);
	}

	// Entry
	public static void
	main(String[] args) {
		if (args.length != 1) {
			System.out.println("Usage: java Main <relation>");
			System.out.println("Example: java Main friend");
			System.exit(0);
		}

		init(args);
		train();
		test();
	}
}
