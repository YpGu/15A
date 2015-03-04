/**
	Implementation of the model in Discovering Latent Classes in Relational Data (Kemp 2004)
	Every node has a hard class label and the probability that two nodes link depends solely on the blocks they belong to 
	Here we do not use sampling to update parameters; instead we use EM.
**/

import java.util.*;

public class BlockModel
{
	public final static int NUM_BLOCKS = 3;
	public final static int MAX_ITER = 10;
	public final static int NUM_INITS = 5;							// init the configuration multiple times, and keep the one with largest likelihood 
	public static int NUM_NODES;
	public static String fileDir, dictDir;

	public static Map<String, Integer> idMap;						// RawID -> newID 
	public static double[][] data;								// network (N*N) 
	public static double[][] testData;
	public static double[][] eta;								// block*block (K*K)
	public static int[] z;									// block assignment (N*1): #Node -> #Block 
												// For eta and z: we need to init multiple times ...
												// , therefore we put an additional dimension
	public static double[][] optEta;
	public static int[] optZ;

	public static Scanner sc;

	/// initialize parameters 
	public static void
	init(
		String[] args,
		int init
	) {
		fileDir = args[0];
		dictDir = args[1];

		idMap = new HashMap<String, Integer>();
		eta = new double[NUM_BLOCKS][NUM_BLOCKS];
		NUM_NODES = FileParser.readDict(dictDir, idMap, "\t");
		data = new double[NUM_NODES][NUM_NODES];
		FileParser.readData(fileDir+".train", idMap, data);
		testData = new double[NUM_NODES][NUM_NODES];
		FileParser.readData(fileDir+".test", idMap, testData);

		z = new int[NUM_NODES];

		sc = new Scanner(System.in);

		// random initialization - z 
		Random rand = new Random(10*init);
		for (int i = 0; i < NUM_NODES; i++) {
			z[i] = rand.nextInt(NUM_BLOCKS);
		}
		// random initialization - eta 
		double[][] m = new double[NUM_BLOCKS][NUM_BLOCKS];
		for (int i = 0; i < NUM_NODES; i++) {
			for (int j = 0; j < NUM_NODES; j++) {
				m[z[i]][z[j]] = 0;
			}
		}
		for (int i = 0; i < NUM_NODES; i++) {
			for (int j = 0; j < NUM_NODES; j++) {
				m[z[i]][z[j]] += data[i][j];
			}
		}
		double[] counter = new double[NUM_BLOCKS];				// counter (K*1): #Block -> num of nodes 
		for (int i = 0; i < NUM_NODES; i++) {
			int block = z[i];
			counter[block] += 1;
		}
		for (int i = 0; i < NUM_BLOCKS; i++) {
			for (int j = 0; j < NUM_BLOCKS; j++) {
				double neg = counter[i] * counter[j] - m[i][j];		// might be too large for almost all the block pairs 
				if (i == j) {
					neg -= counter[i];
				}
		//		System.out.println("m = " + m[i][j] + ", neg = " + neg);
		//		int gu = sc.nextInt();

				eta[i][j] = (m[i][j]+1)/(m[i][j]+neg+2);
			}
		}

		// read data
		FileParser.readData(fileDir, idMap, data);

		System.out.println(Evaluation.calcObj(eta, z));

		return;
	}


	/// estimate parameters using EM algorithm 
	public static double
	train(
	) {
		for (int iter = 0; iter < MAX_ITER; iter++) {

			System.out.println("\tIter = " + iter);

			// convergence checker 
			boolean flag = true;

			// E-step 
			for (int n = 0; n < NUM_NODES; n++) {

				int bestK = -1, preK = z[n];

				double maxChange = -Double.MIN_VALUE, curChange;
	//			System.out.printf("Node %d\t", n);
				for (int k = 0; k < NUM_BLOCKS; k++) {
					curChange = Evaluation.changeInObj(data, eta, z, n, k);
	//				System.out.printf("K = %d, Change = %f\t", k, curChange);
					if (curChange > maxChange) {
						bestK = k;
						maxChange = curChange;
					}						
				}
	//			System.out.printf("Node %d: z changed from %d to %d\t", n, preK, bestK);
	//			int gu = sc.nextInt();

				z[n] = bestK;
				flag = flag && (bestK == preK);
			}

			if (flag) {
				System.out.println("\tz has converged.");
				break;
			}

			System.out.println(Evaluation.calcObj(eta, z));

			// M-step 
			double[][] m = new double[NUM_BLOCKS][NUM_BLOCKS];
			for (int i = 0; i < NUM_NODES; i++) {
				for (int j = 0; j < NUM_NODES; j++) {
					m[z[i]][z[j]] += data[i][j];
				}
			}
			double[] counter = new double[NUM_BLOCKS];				// counter (K*1): #Block -> num of nodes 
			for (int i = 0; i < NUM_NODES; i++) {
				int block = z[i];
				counter[block] += 1;
			}
			for (int i = 0; i < NUM_BLOCKS; i++) {
				for (int j = 0; j < NUM_BLOCKS; j++) {
					double neg = counter[i] * counter[j] - m[i][j];		// might be too large for almost all the block pairs 
					if (i == j) {
						neg -= counter[i];
					}
					eta[i][j] = (m[i][j]+1)/(m[i][j]+neg+2);
					if (m[i][j] == 0) {
						eta[i][j] = 0;
					}

	//				System.out.println("i = " + i + ", j = " + j + ", m = " + m[i][j] + ", neg = " + neg);
				}
			}
			System.out.println(Evaluation.calcObj(eta, z));
		}

		// output z
		double[] counter = new double[NUM_BLOCKS];				// counter (K*1): #Block -> num of nodes 
		for (int i = 0; i < NUM_NODES; i++) {
			int block = z[i];
			counter[block] += 1;
		}
		for (int i = 0; i < counter.length; i++) {
			System.out.printf("\t%d", (int)counter[i]);
		}
		System.out.printf("\n");

		double res = Evaluation.calcObj(eta, z);
		return res;
	}


	/// main entry
	public static void
	main(
		String[] args
	) {
		if (args.length != 2) {
			System.out.println("Usage: java BlockModel <fileDir> <dictDir>");
			System.out.println("Example: java BlockModel ../../data/3k_mention/mention_list_3k.train ../../data/3k_mention/mention_dict_3k");
			System.exit(0);
		}

		double maxObj = -Double.MAX_VALUE, curObj;
		int init = 0;
		double[] counter = new double[NUM_BLOCKS];				// counter (K*1): #Block -> num of nodes 
		boolean zeroFlag = true;
		while (true) {
			// TODO: check if any block contains no nodes 
			if (init > NUM_INITS) {
				break;
			}

			System.out.println("Initialization #" + init);
			init(args, init);
			curObj = train();
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
			for (int i = 0; i < NUM_NODES; i++) {
				int block = z[i];
				counter[block] += 1;
			}
			for (int k = 0; k < NUM_BLOCKS; k++) {
				zeroFlag = zeroFlag && counter[k];
			}
			if (zeroFlag) {
				break;
			}
		}

		// output z
		for (int i = 0; i < NUM_NODES; i++) {
			int block = optZ[i];
			counter[block] += 1;
		}
		for (int i = 0; i < counter.length; i++) {
			System.out.printf("\t%d", (int)counter[i]);
		}
		System.out.printf("\n");

		FileParser.output("./res/z_"+NUM_BLOCKS, optZ);
		FileParser.output("./res/eta_"+NUM_BLOCKS, optEta);

		return;
	}
}
