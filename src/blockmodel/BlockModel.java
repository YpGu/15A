/**
	Implementation of the model in Discovering Latent Classes in Relational Data (Kemp 2004)
	Every node has a hard class label and the probability that two nodes link depends solely on the blocks they belong to 
	Here we do not use sampling to update parameters; instead we use EM.
**/

import java.util.*;

public class BlockModel
{
	public final static int NUM_BLOCKS = 10;
	public final static int MAX_ITER = 100;
	public final static int NUM_INITS = 5;							// init the configuration multiple times, and keep the one with largest likelihood 
	public static String fileDir, dictDir;

	public static Map<String, Integer> idMap;						// RawID -> newID 
	public static Map<String, Integer> z;							// Block Assignment
	public static double[][] eta;								// Block Parameters (K*K) 

	public static SparseMatrix trainData, testData;

	public static Map<String, Integer> optZ;
	public static double[][] optEta;

	public static Scanner sc;

	/// initialize parameters 
	public static void
	init(String[] args, int init) {

		fileDir = args[0];
		dictDir = args[1];

		eta = new double[NUM_BLOCKS][NUM_BLOCKS];

		trainData = new SparseMatrix();
		testData = new SparseMatrix();
		FileParser.readData(dictDir, fileDir + ".train", trainData);
		FileParser.readData(dictDir, fileDir + ".test", testData);

		z = new HashMap<String, Integer>();
		optZ = new HashMap<String, Integer>();

		// random initialization - z 
		Random rand = new Random(10*init+1);
		for (String s: trainData.getDict()) {
			z.put(s, rand.nextInt(NUM_BLOCKS));
		}

		// random initialization - eta 
		double[][] m = new double[NUM_BLOCKS][NUM_BLOCKS];
		for (Map.Entry<Tuple<String, String>, Double> e: trainData.getMat().entrySet()) {
			String x = e.getKey().getX();
			String y = e.getKey().getY();
			double v = e.getValue();
			int zx = z.get(x);
			int zy = z.get(y);
			m[zx][zy] += v;
		}
		double[] counter = new double[NUM_BLOCKS];				// counter (K*1): #Block -> num of nodes 
		for (String i: trainData.getDict()) {					// TODO: if we iterate through trainData.dict, some nodes will never appear in the actual training dataset (in test only) 
			int block = z.get(i);
			counter[block] += 1;
		}
		for (int i = 0; i < NUM_BLOCKS; i++) {
			for (int j = 0; j < NUM_BLOCKS; j++) {
				double pb = counter[i] * counter[j];			// might be too large for almost all the block pairs 
				if (i == j) {
					pb -= counter[i];
				}
		//		System.out.println("m = " + m[i][j] + ", neg = " + neg);
		//		int gu = sc.nextInt();
				if (pb != 0) {
					eta[i][j] = (m[i][j]+1)/(pb+2);
				}
				else {
					eta[i][j] = 0;
				}
			}
		}

		sc = new Scanner(System.in);

		System.out.println("Objective at init = " + Evaluation.calcObj(trainData, eta, z));

		return;
	}


	/// estimate parameters using EM algorithm 
	public static double
	train(SparseMatrix data) {
		for (int iter = 0; iter < MAX_ITER; iter++) {
			System.out.println("\t---- Iter = " + iter);

			// convergence checker 
			boolean flag = true;

			// E-step 
			long sTime = System.currentTimeMillis();

			for (String s: trainData.getDict()) {
				int bestK = z.get(s), preK = z.get(s);
				double maxChange = -Double.MIN_VALUE, curChange;

				long aTime = System.currentTimeMillis();

				for (int k = 0; k < NUM_BLOCKS; k++) {
					if (k != preK) {
						curChange = Evaluation.changeInObj(data, eta, z, s, k);
						if (curChange > maxChange) {
							bestK = k;
							maxChange = curChange;
						}
					}
				}

				long bTime = System.currentTimeMillis();
//				System.out.println("Time: " + (bTime-aTime));
//				Scanner sc = new Scanner(System.in);
//				int gu = sc.nextInt();

				if (bestK == -1) {
					assert false;
				}
					
				z.put(s, bestK);

				flag = flag && (bestK == preK);
			}

			if (flag) {
				System.out.println("\tz has converged.");
				break;
			}

//// check blocks
	double[] counter = new double[NUM_BLOCKS];				// counter (K*1): #Block -> num of nodes 
	for (Map.Entry<String, Integer> i: z.entrySet()) {
		int block = i.getValue();
		counter[block] += 1;
	}
	for (int i = 0; i < counter.length; i++) {
		System.out.printf("\t%d", (int)counter[i]);
	}
	System.out.printf("\n");
///////

			long fTime = System.currentTimeMillis();
//			System.out.println("Time = " + (fTime-sTime));

			System.out.println("\tObj = " + Evaluation.calcObj(trainData, eta, z));

			// M-step 
			double[][] m = new double[NUM_BLOCKS][NUM_BLOCKS];

			for (Map.Entry<Tuple<String, String>, Double> e: trainData.getMat().entrySet()) {
				String x = e.getKey().getX();
				String y = e.getKey().getY();
				double v = e.getValue();
				int zx = z.get(x);
				int zy = z.get(y);
				m[zx][zy] += 1;
			}

			counter = null;
			counter = new double[NUM_BLOCKS];				// counter (K*1): #Block -> num of nodes 
			for (Map.Entry<String, Integer> i: z.entrySet()) {
				int block = i.getValue();
				counter[block] += 1;
			}


			for (int i = 0; i < NUM_BLOCKS; i++) {
				for (int j = 0; j < NUM_BLOCKS; j++) {
					double pb = counter[i] * counter[j];			// might be too large for almost all the block pairs 
					if (i == j) {
						pb -= counter[i];
					}
					if (pb != 0) {
						eta[i][j] = (m[i][j]+1)/(pb+2);
					}
					else {
						eta[i][j] = 0;
					}

//					System.out.println("i = " + i + ", j = " + j + ", eta = " + eta[i][j]);
				}
			}

			System.out.println("\tObj = " + Evaluation.calcObj(trainData, eta, z));
		}

		// output z
		double[] counter = new double[NUM_BLOCKS];				// counter (K*1): #Block -> num of nodes 
		for (Map.Entry<String, Integer> i: z.entrySet()) {
			int block = i.getValue();
			counter[block] += 1;
		}
		for (int i = 0; i < counter.length; i++) {
			System.out.printf("\t%d", (int)counter[i]);
		}
		System.out.printf("\n");

		double res = Evaluation.calcObj(trainData, eta, z);
		return res;
	}


	/// main entry
	public static void
	main(String[] args) {
		if (args.length != 2) {
			System.out.println("Usage: java BlockModel <fileDir> <dictDir>");
			System.out.println("Example: java BlockModel ../../data/3k_mention/mention_list_3k ../../data/3k_mention/mention_dict_3k");
			System.exit(0);
		}

		double maxObj = -Double.MAX_VALUE, curObj;
		int init = 0;
		double[] counter = new double[NUM_BLOCKS];				// counter (K*1): #Block -> num of nodes 
		boolean zeroFlag = true;
		while (true) {
			if (init > NUM_INITS) {
				break;
			}

			System.out.println("------ Initialization #" + init);
			init(args, init);
			curObj = train(trainData);
			System.out.println("Objective function = " + curObj);
			if (curObj > maxObj) {
				maxObj = curObj;
				optEta = eta;
				optZ = z;
			}
			init += 1;

			// check empty blocks (for [current] block): break if all blocks are non-empty 
			for (int k = 0; k < NUM_BLOCKS; k++) {
				counter[k] = 0;
			}
			for (Map.Entry<String, Integer> i: z.entrySet()) {
				int block = i.getValue();
				counter[block] += 1;
			}
			for (int k = 0; k < NUM_BLOCKS; k++) {
				zeroFlag = zeroFlag && (counter[k] != 0);
			}
			if (zeroFlag) {
				break;
			}
		}

		// output z
		for (int k = 0; k < NUM_BLOCKS; k++) {
			counter[k] = 0;
		}
		for (Map.Entry<String, Integer> i: z.entrySet()) {
			int block = i.getValue();
			counter[block] += 1;
		}
		for (int i = 0; i < counter.length; i++) {
			System.out.printf("\t%d", (int)counter[i]);
		}
		System.out.printf("\n");

//		FileParser.output("./res/z_"+NUM_BLOCKS, optZ);
//		FileParser.output("./res/eta_"+NUM_BLOCKS, optEta);

		return;
	}
}
