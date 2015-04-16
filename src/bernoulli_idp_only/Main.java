/**
	03/30/2015 

	Main.java: implement the unified model with two multinomial mixtures

	The first mixture is background;
	The second mixture is ideal point;

	If we want to build the model purely based on background, set USE_BKG to true and USE_IPM to false;
	If we want to build the model purely based on ideal point, set USE_BKG to false and USE_IPM to true

	If we want to use the bias term in ideal point model, set USEB to true; otherwise, set USEB to false

**/

import java.util.*;

public class Main
{
	// Configuration
	public static int N;						// Number of Users
	public static final int MAX_ITER = 100;				// Maximum Number of Iterations 
	public static final double THRESHOLD = Math.pow(10,-5);
	public static final boolean USE_IPM = true;
	public static final boolean USEB = true;			// true if we use p[i]*q[j]+b[j]; fase if we use p[i]*q[j] 

	public static SparseMatrix<Integer> trainData, testData;
	public static SparseMatrix<Integer> trainDataNeg, testDataNeg;
	public static Map<String, Integer> dict;
	public static Map<Integer, String> invDict;

	// Model Parameters
	public static double[] p, q, b;					// N * 1

	// Initialization and Read Data 
	public static void
	init(String[] args) {
		// Init
		trainData = new SparseMatrix<Integer>(); testData = new SparseMatrix<Integer>();
		trainDataNeg = new SparseMatrix<Integer>(); testDataNeg = new SparseMatrix<Integer>();
		dict = new HashMap<String, Integer>(); invDict = new HashMap<Integer, String>();
		// Read Data
		String num = "3k", rel = args[0];
		String dictDir = "../../data/" + num + "_" + rel + "/" + rel + "_dict_" + num;
		String trainDataDir = "../../data/" + num + "_" + rel + "/" + rel + "_list_" + num + ".train";
		String testDataDir = "../../data/" + num + "_" + rel + "/" + rel + "_list_" + num + ".test";
		String trainDataDirNeg = "../../data/" + num + "_" + rel + "/n_" + rel + "_list_" + num + ".train";
		String testDataDirNeg = "../../data/" + num + "_" + rel + "/n_" + rel + "_list_" + num + ".test";
		dict = FileParser.readVocabulary(dictDir);
		invDict = FileParser.readInverseVocabulary(dictDir);
		FileParser.readData(trainData, trainDataDir, dict); FileParser.readData(testData, testDataDir, dict);
		FileParser.readData(trainDataNeg, trainDataDirNeg, dict); FileParser.readData(testDataNeg, testDataDirNeg, dict);
		N = dict.size();
		// Initialize Parameters 
		p = new double[N]; q = new double[N]; b = new double[N];
		FileParser.readInit(p, "./init/P");
		FileParser.readInit(q, "./init/Q");
		FileParser.readInit(b, "./init/B");
/*
		Random rand = new Random(10);
		for (int i = 0; i < N; i++) {
			double pqRange = 1;
			p[i] = -0.5 * pqRange + pqRange * rand.nextDouble();
			q[i] = -0.5 * pqRange + pqRange * rand.nextDouble();
			b[i] = -0.5 + rand.nextDouble();
		}
*/
	}

	// Train
	public static void
	train() {
		double oldLikelihood = -1;
		for (int iter = 0; iter < MAX_ITER; iter++) {
			System.out.println("----- Iteration " + iter + " -----");
			double likelihood = Update.update(trainData, trainDataNeg, p, q, b, iter);

			FileParser.output("./param/p", p, invDict);
			FileParser.output("./param/q", q, invDict);
			FileParser.output("./param/b", b, invDict);

			double rate = Math.abs((likelihood-oldLikelihood)/oldLikelihood);
			if (rate < THRESHOLD && iter != 0)
				break;
		}

		return;
	}

	// Test
	public static void
	test() {
		// todo
		System.out.println("Training (all):");
		Evaluation.auroc(trainData, trainDataNeg, p, q, b, 1);
		Evaluation.auprc(trainData, trainDataNeg, p, q, b, 1);

		System.out.println("Testing (all):");
		Evaluation.auroc(testData, testDataNeg, p, q, b, 1);
		Evaluation.auprc(testData, testDataNeg, p, q, b, 1);

/*		System.out.println("Training (aver):");
		Evaluation.auroc2(trainData, trainDataNeg, alpha, beta, pi, p, q, b, gamma, phi, varphi, 1);
		Evaluation.auprc2(trainData, trainDataNeg, alpha, beta, pi, p, q, b, gamma, phi, varphi, 1);
		System.out.println("Testing (aver):");
		Evaluation.auroc2(testData, testDataNeg, alpha, beta, pi, p, q, b, gamma, phi, varphi, 2);
		Evaluation.auprc2(testData, testDataNeg, alpha, beta, pi, p, q, b, gamma, phi, varphi, 2);
*/
		System.out.println("Classification:");
		Evaluation.partyClassify(p, q, b, invDict);
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
