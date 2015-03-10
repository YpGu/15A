/**
	Implementation of the model in Discovering Latent Classes in Relational Data (Kemp 2004);
	Every node has a hard class label and the probability that two nodes link depends solely on the blocks they belong to;
	Here we do not use sampling to update parameters; instead we use EM;

	Time complexity may probably be an issue when the number of nodes is big;

	Todo: dictionary of the training data and all.
**/

import java.util.*;

public class BlockModel
{
	public static int NUM_BLOCKS;							// Number of Blocks (pre-defined) 
	public final static int MAX_ITER = 100;						// Maximum number of iterations until convergence 
	public final static int NUM_INITS = 5;						// init the configuration multiple times, and keep the one with largest likelihood 
	public final static boolean WRITE = false;
	public static String fileDir, dictDir;

	public static SparseMatrix trainData, testData;					// Graph data (training/testing) 

	public static Map<String, Double> vOut;						// idp - out Ideal Point 
	public static Map<String, Double> vIn;
	public static Map<String, Double> vBias;
	public static Map<String, Integer> z;						// bkg - Block Assignment
	public static double[][] eta;							// bkg - Block Parameters (K*K) 
	public static Map<String, Double> pi;						// weight of idp mixture 

	public static Map<String, Double> optOut;
	public static Map<String, Double> optIn;
	public static Map<String, Double> optBias;
	public static Map<String, Double> optPi;
	public static Map<String, Integer> optZ;					// [Optimal] Block Assignment 
	public static double[][] optEta;						// [Optimal] Block Parameters 

	public static Scanner sc;


	/// constructor 
	public BlockModel() {
		optOut = new HashMap<String, Double>();
		optIn = new HashMap<String, Double>();
		optBias = new HashMap<String, Double>();
		optPi = new HashMap<String, Double>();
		optZ = new HashMap<String, Integer>();
	}


	/// initialize parameters 
	public static void
	init(String[] args, int init) {

		fileDir = args[0];
		dictDir = args[1];
		NUM_BLOCKS = Integer.parseInt(args[2]);

		trainData = new SparseMatrix();
		testData = new SparseMatrix();
		FileParser.readData(dictDir, fileDir + ".train", trainData);
		FileParser.readData(dictDir, fileDir + ".test", testData);

		Random rand = new Random(10*init+1);

		// random initialization - in/out/bias & pi 
		vOut = new HashMap<String, Double>();
		vIn = new HashMap<String, Double>();
		vBias = new HashMap<String, Double>();
		pi = new HashMap<String, Double>();
		for (String s: trainData.getDict()) {
			vOut.put(s, (rand.nextDouble()-0.5)*1);
			vIn.put(s, (rand.nextDouble()-0.5)*1);
			vBias.put(s, (rand.nextDouble()-0.5)*1);
			pi.put(s, rand.nextDouble() * 0.2 + 0.4);
		}

		// random initialization - z/eta  
		z = new HashMap<String, Integer>();
		for (String s: trainData.getDict()) {
			z.put(s, rand.nextInt(NUM_BLOCKS));
		}
		eta = new double[NUM_BLOCKS][NUM_BLOCKS];
		UpdateBM.updateParamHard(trainData, z, eta);

		System.out.println("Objective at init = " + Evaluation.calcObj(trainData, eta, z));

		return;
	}


	/// estimate parameters using two-step EM-like algorithm 
	public static double
	train(SparseMatrix data) {
		for (int iter = 0; iter < MAX_ITER; iter++) {
			System.out.println("\t---- Iter = " + iter + " ----");
			if (UpdateBM.updateHard(data, z, eta)) {
				break;
			}
		}

		// output z
		System.out.println("    Final Block Assignments:");
		UpdateBM.checkBlocks(z, NUM_BLOCKS);

		return Evaluation.calcObj(trainData, eta, z);
	}


	/// main entry
	public static void
	main(String[] args) {
		if (args.length != 3) {
			System.out.println("Usage: java BlockModel <fileDir> <dictDir> <#blocks>");
			System.out.println("Example: java BlockModel ../../data/3k_mention/mention_list_3k ../../data/3k_mention/mention_dict_3k 10");
			System.exit(0);
		}

		double maxObj = -Double.MAX_VALUE, curObj;
		int init = 1;
		while (true) {
			System.out.println("------ Initialization #" + init + " ------");
			init(args, init);
			curObj = train(trainData);
			System.out.println("Objective function = " + curObj);
			if (curObj > maxObj) {
				maxObj = curObj; optEta = eta;
				System.out.println("optZ updated!");
				optZ = new HashMap<String, Integer>(z);
			}

			if (!UpdateBM.existEmptyBlock(z, NUM_BLOCKS, init) && init >= NUM_INITS) {
				break;
			}

			init += 1;
		}

		// output z
		UpdateBM.checkBlocks(optZ, NUM_BLOCKS);

		if (WRITE) {
			FileParser.output("./res/z_"+NUM_BLOCKS, optZ);
			FileParser.output("./res/eta_"+NUM_BLOCKS, optEta);
		}

		return;
	}
}
