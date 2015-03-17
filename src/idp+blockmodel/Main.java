/**
	Implementation of the model in Discovering Latent Classes in Relational Data (Kemp 2004);
	Every node has a hard class label and the probability that two nodes link depends solely on the blocks they belong to;
	Here we do not use sampling to update parameters; instead we use EM;

	Time complexity may probably be an issue when the number of nodes is big;

	Todo: dictionary of the training data and all.
**/

import java.util.*;

public class Main
{
	public static int NUM_BLOCKS;							// Number of Blocks (pre-defined) 
	public final static int MAX_ITER = 30;						// Maximum number of iterations until convergence 
	public final static int NUM_INITS = 5;						// init the configuration multiple times, and keep the one with largest likelihood 
	public final static boolean WRITE = true;

	public static double sw;							// sample weight or 1 
	public static double reg;							// regularization coefficient 

	public static SparseMatrix trainPositiveData, trainNegativeData;		// Graph data (training: existing/non-existing) 
	public static SparseMatrix testPositiveData, testNegativeData;			// Graph data (testing: existing/non-existing) 

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
	public Main() {
		optOut = new HashMap<String, Double>();
		optIn = new HashMap<String, Double>();
		optBias = new HashMap<String, Double>();
		optPi = new HashMap<String, Double>();
		optZ = new HashMap<String, Integer>();
	}


	/// initialize parameters 
	public static void
	init(String[] args, int init) {

		String rel = args[0];
		String num = "3k";
		NUM_BLOCKS = Integer.parseInt(args[1]);
		reg = Double.parseDouble(args[2]);

		trainPositiveData = new SparseMatrix(); trainNegativeData = new SparseMatrix();
		testPositiveData = new SparseMatrix(); testNegativeData = new SparseMatrix();
		FileParser.readData("../../data/" + num + "_" + rel + "/" + rel + "_dict_" + num, "../../data/" + num + "_" + rel + "/" + rel + "_list_" + num + ".train", trainPositiveData);
		FileParser.readData("../../data/" + num + "_" + rel + "/" + rel + "_dict_" + num, "../../data/" + num + "_" + rel + "/n_" + rel + "_list_" + num + ".train", trainNegativeData);
		FileParser.readData("../../data/" + num + "_" + rel + "/" + rel + "_dict_" + num, "../../data/" + num + "_" + rel + "/" + rel + "_list_" + num + ".test", testPositiveData);
		FileParser.readData("../../data/" + num + "_" + rel + "/" + rel + "_dict_" + num, "../../data/" + num + "_" + rel + "/n_" + rel + "_list_" + num + ".test", trainNegativeData);

/*
		for (String x: trainPositiveData.getDict()) {
			System.out.println(x);
			System.out.println(trainPositiveData.getRow(x).size());
			System.out.println(trainNegativeData.getRow(x).size());

			Scanner sc = new Scanner(System.in);
			int gu = sc.nextInt();
		}
*/
		// sample weight 
		sw = (double)trainPositiveData.getDictSize() * (trainPositiveData.getDictSize()-1) / (double)trainPositiveData.getSize() - 1;
//		sw = 1;

		Random rand = new Random(10*init+1);

		// random initialization - in/out/bias & pi 
		vOut = new HashMap<String, Double>();
		vIn = new HashMap<String, Double>();
		vBias = new HashMap<String, Double>();
		pi = new HashMap<String, Double>();
		for (String s: trainPositiveData.getDict()) {
			vOut.put(s, (rand.nextDouble()-0.5)*1);
			vIn.put(s, (rand.nextDouble()-0.5)*1);
			vBias.put(s, (rand.nextDouble()-0.5)*1);
			pi.put(s, rand.nextDouble() * 0.2 + 0.4);
		}

		// random initialization - z/eta  
		z = new HashMap<String, Integer>();
		for (String s: trainPositiveData.getDict()) {
			z.put(s, rand.nextInt(NUM_BLOCKS));
		}
		eta = new double[NUM_BLOCKS][NUM_BLOCKS];
		UpdateBM.updateParamHardAtInit(trainPositiveData, trainNegativeData, z, eta, sw);

		System.out.println("Objective at init = " + Evaluation.calcObj(trainPositiveData, trainNegativeData, eta, z, vOut, vIn, vBias, pi, sw, reg));

		return;
	}


	/// update parameters: for multiple runs, we keep the best set of parameters 
	public static void
	saveParam() {
		optPi = pi;
		optEta = eta;
		optOut = vOut; optIn = vIn; optBias = vBias;
		optZ = z;

		return;
	}


	/// estimate parameters using two-step EM-like algorithm 
	public static double
	train() {
		for (int iter = 0; iter < MAX_ITER; iter++) {
			System.out.println("\n\t---- Iter = " + iter + " ----");
			double lr = 0.00001;
			if (Update.update(trainPositiveData, trainNegativeData, vOut, vIn, vBias, z, eta, pi, sw, reg, lr)) {
				break;
			}
		}

		// output z
		System.out.println("    Final Block Assignments:");
		UpdateBM.checkBlocks(z, NUM_BLOCKS);

		return Evaluation.calcObj(trainPositiveData, trainNegativeData, eta, z, vOut, vIn, vBias, pi, sw, reg);
//		return 1;
	}


	/// main entry
	public static void
	main(String[] args) {
		if (args.length != 3) {
			System.out.println("Usage: java Main <relation> <#blocks> <reg>");
			System.out.println("Example: java Main friend 10 0.01");
			System.exit(0);
		}

		double maxObj = -Double.MAX_VALUE, curObj;
		int init = 1;
		while (true) {
			System.out.println("\n------ Initialization #" + init + " ------");
			init(args, init);
			curObj = train();
			System.out.println("Objective function = " + curObj);
			if (curObj > maxObj) {
				maxObj = curObj;
				System.out.println("Optimal parameter updated!");
				saveParam();
			}

			if (!UpdateBM.existEmptyBlock(optZ, NUM_BLOCKS) && init >= NUM_INITS) {
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

		Evaluation.auroc(trainPositiveData, optPi, optZ, optEta, optOut, optIn, optBias);

		return;
	}
}
