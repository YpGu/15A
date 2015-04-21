/**
	04/20/2015 
	Main.java: implement the unified model with two multinomial mixtures

	The first mixture is background;
	The second mixture is ideal point;

	If we want to build the model purely based on background, set USE_BKG to true and USE_IPM to false;
	If we want to build the model purely based on ideal point, set USE_BKG to false and USE_IPM to true

	If we want to use the bias term in ideal point model, set USEB to true; otherwise, set USEB to false

	04/20/2015
	Problem exists in how to evaluate p(i--j), since it's unfair to directly compare the two mixtures 
**/

import java.util.*;
import java.io.*;

public class Main
{
	// Configuration
	public static final int K = 5;					// Number of Latent Features
	public static int N;						// Number of Users
	public static final int MAX_ITER = 10;				// Maximum Number of Iterations 
	public static final double THRESHOLD = Math.pow(10,-5);
	public static final boolean USE_BKG = true;
	public static final boolean USE_IPM = true;
	public static final boolean USEB = false;			// true if we use p[i]*q[j]+b[j]; fase if we use p[i]*q[j] 

	public static SparseMatrix<Integer> trainData, testData;
	public static SparseMatrix<Integer> trainDataNeg, testDataNeg;
	public static Map<String, Integer> dict;
	public static Map<Integer, String> invDict;

	// Model Parameters
	public static double[] pi;					// 2 * 1
	public static double[] gamma;					// N * 1
	public static double[] p, q, b;					// N * 1
	public static double[][] theta, beta;				// N * K; K * N 

	public static void
	init(String[] args, int option) {
		// Initialize Data 
		trainData = new SparseMatrix<Integer>(); testData = new SparseMatrix<Integer>();
		trainDataNeg = new SparseMatrix<Integer>(); testDataNeg = new SparseMatrix<Integer>();
		dict = new HashMap<String, Integer>(); invDict = new HashMap<Integer, String>();

		// Initialize Parameters 
		N = InitReader.init(args, trainData, testData, trainDataNeg, testDataNeg, dict, invDict, option);

		pi = new double[2]; gamma = new double[N];
		p = new double[N]; q = new double[N]; b = new double[N];
		theta = new double[N][K]; beta = new double[K][N];
		InitReader.init(pi, gamma, p, q, b, theta, beta, USE_BKG, USE_IPM, K, option);

		System.out.println("Size of invDict = " + invDict.size());

		return;
	}

	// Train
	public static void
	train() {
		double oldLikelihood = -1;
		for (int iter = 0; iter < MAX_ITER; iter++) {
			System.out.println("----- Iteration " + iter + " -----");
			double likelihood = Update.update(trainData, pi, gamma, theta, beta, p, q, b, USE_BKG, USE_IPM, iter);

			FileParser.output("./param/p", p, invDict);
			FileParser.output("./param/q", q, invDict);
			FileParser.output("./param/b", b, invDict);
			FileParser.output("./param/pi", pi);

			double rate = Math.abs((likelihood-oldLikelihood)/oldLikelihood);
			if (rate < THRESHOLD && iter != 0)
				break;
		}

		return;
	}

	// Test
	public static void
	test() {
		System.out.println("pi[0] = " + pi[0] + " pi[1] = " + pi[1]);
		// todo
		System.out.println("Training (all):");
		Evaluation.auroc(trainData, trainDataNeg, pi, theta, beta, p, q, b, 1);
		Evaluation.auprc(trainData, trainDataNeg, pi, theta, beta, p, q, b, 1);

		System.out.println("Testing (all):");
		Evaluation.auroc(testData, testDataNeg, pi, theta, beta, p, q, b, 1);
		Evaluation.auprc(testData, testDataNeg, pi, theta, beta, p, q, b, 1);

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
		if (args.length >= 2) {
			System.out.println("Usage: java Main <relation>");
			System.out.println("Example: java Main friend");
			System.out.println("If using java Main, settings.ini will be used as input");
			System.exit(0);
		}

		if (args.length == 0) {
			File f = new File("settings.ini");
			if (f.exists()) {
				System.out.println("Reading settings.ini ...");
			}
			else {
				System.out.println("settings.ini does not exist!");
				System.out.println("Usage: java Main <relation>");
				System.out.println("Example: java Main friend");
				System.exit(0);
			}	
		}

		init(args, args.length);

		train();
		test();
	}
}
