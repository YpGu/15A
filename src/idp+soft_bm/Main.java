/**
	First we implement the model in Discovering Latent Classes in Relational Data (Kemp 2004);
	Every node has a hard class label and the probability that two nodes link depends solely on the blocks they belong to;
	Here we do not use sampling to update parameters; instead we use EM;
	Then we combine the blockmodel with ideal point model;
	The individual \pi value is also trained by EM;

	Note:
	Time complexity may probably be an issue when the number of nodes is big;
	Learning rate should be set manually;
**/

import java.util.*;

public class Main
{
	public static int NUM_BLOCKS;							// Number of Blocks (pre-defined) 
	public final static int MAX_ITER = 30;						// Maximum number of iterations until convergence 
	public final static int NUM_INITS = 1;						// init the configuration multiple times, and keep the one with largest likelihood 
	public final static boolean WRITE = true;					// whether save to file

	public static double sw;							// sample weight or 1 
	public static double reg;							// regularization coefficient 

	public static SparseMatrix trainPositiveData, trainNegativeData;		// Graph data (training: existing/non-existing) 
	public static SparseMatrix testPositiveData, testNegativeData;			// Graph data (testing: existing/non-existing) 

	// Parameter Set 
	public static Map<String, Double> vOut;						// idp - outgoing Ideal Point 
	public static Map<String, Double> vIn;						// idp - incoming Ideal Point
	public static Map<String, Double> vBias;					// bias 
	public static Map<String, double[]> theta;
	public static double[][] eta;							// bkg - Block Parameters (K*K) 
	public static Map<String, Double> pi;						// global weight of [idp mixture]

	// [Optimal] Parameter Set
	public static Map<String, Double> optOut;
	public static Map<String, Double> optIn;
	public static Map<String, Double> optBias;
	public static Map<String, double[]> optTheta;
	public static double[][] optEta;
	public static Map<String, Double> optPi;


	/// constructor 
	public Main() {
		optOut = new HashMap<String, Double>();
		optIn = new HashMap<String, Double>();
		optBias = new HashMap<String, Double>();
		optPi = new HashMap<String, Double>();
		optTheta = new HashMap<String, double[]>();
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
		FileParser.readData("../../data/" + num + "_" + rel + "/" + rel + "_dict_" + num, "../../data/" + num + "_" + rel + "/n_" + rel + "_list_" + num + ".test", testNegativeData);

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
//			pi.put(s, 0.0);
		}

		// random initialization - z/eta  
		theta = new HashMap<String, double[]>();
		eta = new double[NUM_BLOCKS][NUM_BLOCKS];
		UpdateBM.updateParamAtInit(trainPositiveData, trainNegativeData, theta, eta, sw);

		System.out.println("Objective at init = " + Evaluation.calcObj(trainPositiveData, trainNegativeData, theta, eta, vOut, vIn, vBias, pi, sw, reg));

		return;
	}


	/// save the optimal parameters: for multiple runs, we keep the best set of parameters 
	public static void
	saveParam() {
		optPi = pi;
		optEta = eta;
		optOut = vOut; optIn = vIn; optBias = vBias;
		optTheta = theta;

		return;
	}


	/// save to file
	public static void 
	writeToFile() {
		FileParser.outputArray("./res/theta_" + NUM_BLOCKS, optTheta);
		FileParser.output("./res/eta_" + NUM_BLOCKS, optEta);
		FileParser.output("./res/pi_" + NUM_BLOCKS, optPi);
		FileParser.output("./res/out_" + NUM_BLOCKS, optOut);
		FileParser.output("./res/in_" + NUM_BLOCKS, optIn);
		FileParser.output("./res/bias_" + NUM_BLOCKS, optBias);
	}


	/// estimate parameters using two-step EM-like algorithm 
	public static double
	train() {
		double preObj = -1, curObj;
		for (int iter = 0; iter < MAX_ITER; iter++) {
			System.out.println("\n\t---- Iter = " + iter + "/" + MAX_ITER + " ----");
			double lr = 0.00005;
			curObj = Update.update(trainPositiveData, trainNegativeData, vOut, vIn, vBias, theta, eta, pi, sw, reg, lr);
			if (iter != 0) {
				double rate = -(curObj-preObj)/preObj;
				System.out.println("\tObjective function = " + curObj + " rate = " + rate);
//				if (rate < Math.pow(10,-6)) {
//					break;
//				}
			}
			preObj = curObj;
		}

		return Evaluation.calcObj(trainPositiveData, trainNegativeData, theta, eta, vOut, vIn, vBias, pi, sw, reg);
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
			System.out.println("\n------ Initialization #" + init + "/" + NUM_INITS + " ------");
			init(args, init);
			curObj = train();
			System.out.println("Objective function = " + curObj);
			if (curObj > maxObj) {
				maxObj = curObj;
				System.out.println("Optimal parameter updated!");
				saveParam();
			}

			if (init >= NUM_INITS) {
				break;
			}

//			if (init >= 2*NUM_INITS) {
//				break;
//			}

			init += 1;
		}

		if (WRITE) {
			writeToFile();
		}

		System.out.println("\n------ Evaluation ------");
		System.out.println("\t------ Training ------");
		Evaluation.auroc(trainPositiveData, trainNegativeData, optPi, optTheta, optEta, optOut, optIn, optBias, 1);
		System.out.println("\t------ Testing ------");
		Evaluation.auroc(testPositiveData, testNegativeData, optPi, optTheta, optEta, optOut, optIn, optBias, 2);

		return;
	}
}
