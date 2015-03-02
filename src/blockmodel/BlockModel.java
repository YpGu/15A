/**
	Implementation of the model in Discovering Latent Classes in Relational Data (Kemp 2004)
	Every node has a hard class label and the probability that two nodes link depends solely on the blocks they belong to
	Here we do not use sampling to update parameters; instead we use EM.
**/

import java.util.*;

public class BlockModel
{
//	public static ArrayList<
	public final static int NUM_BLOCKS = 10;
	public final static int MAX_ITER = 1000;
	public static int NUM_NODES;
	public static String fileDir, dictDir;

	public static Map<String, Integer> idMap;						// RawID -> newID 
	public static double[][] data;								// network (N*N) 
	public static double[][] eta;								// block*block (K*K)
	public static int[] z;									// block assignment (N*1): #Node -> #Block 


	/// initialize parameters 
	public static void
	init(
		String[] args
	) {
		fileDir = args[0];
		dictDir = args[1];

		idMap = new HashMap<String, Integer>();
		eta = new double[NUM_BLOCKS][NUM_BLOCKS];
		NUM_NODES = DataReader.readDict(dictDir, idMap, "\t");
		data = new double[NUM_NODES][NUM_NODES];
		z = new int[NUM_NODES];


		// random initialization - z 
		Random rand = new Random(0);
		for (int i = 0; i < NUM_NODES; i++) {
			z[i] = rand.nextInt(NUM_BLOCKS);
		}
		// random initialization - eta 
		double[][] m = new double[NUM_BLOCKS][NUM_BLOCKS];
		for (int i = 0; i < NUM_NODES; i++) {
			for (int j = 0; j < NUM_NODES; j++) {
				m[z[i]][z[j]] += 1;
			}
		}
		double[] counter = new double[NUM_BLOCKS];				// counter (K*1): #Block -> num of nodes 
		for (int i = 0; i < NUM_NODES; i++) {
			int block = z[i];
			counter[block] += 1;
		}
		for (int i = 0; i < NUM_BLOCKS; i++) {
			for (int j = 0; j < NUM_BLOCKS; j++) {
				double neg = counter[i] * counter[j];			// might be too large for almost all the block pairs 
				if (i == j) {
					neg -= counter[i];
				}
				eta[i][j] = (m[i][j]+1)/(m[i][j]+neg+2);
			}
		}

		// read data
		DataReader.readData(fileDir, idMap, data);

		return;
	}


	/// estimate parameters using EM algorithm 
	public static void
	train(
	) {
		for (int iter = 0; iter < MAX_ITER; iter++) {

			System.out.println("Iter = " + iter);

			// check convergence
			boolean flag = true;

			// E-step 
			for (int n = 0; n < NUM_NODES; n++) {
//				System.out.println("===N = " + n);
				double maxObj = -Double.MAX_VALUE, curObj;
				double maxChange = -Double.MAX_VALUE, curChange;
				int bestK = -1, preK = z[n];
				for (int k = 0; k < NUM_BLOCKS; k++) {
//					z[n] = k;
//					curObj = Evaluation.calcObj(eta, z);			// need 6-7 ms 
//					if (curObj > maxObj) {
//						bestK = k;
//						maxObj = curObj;
//					}
					curChange = Evaluation.changeInObj(data, eta, z, n, k);
					if (curChange > maxChange) {
						bestK = k;
						maxChange = curChange;
					}						
				}
				z[n] = bestK;
				flag = flag || (bestK == preK);
			}

			if (flag) {
				System.out.println("z has converged.");
				break;
			}

			// M-step 
			double[][] m = new double[NUM_BLOCKS][NUM_BLOCKS];
			for (int i = 0; i < NUM_NODES; i++) {
				for (int j = 0; j < NUM_NODES; j++) {
					m[z[i]][z[j]] += 1;
				}
			}
			double[] counter = new double[NUM_BLOCKS];				// counter (K*1): #Block -> num of nodes 
			for (int i = 0; i < NUM_NODES; i++) {
				int block = z[i];
				counter[block] += 1;
			}
			for (int i = 0; i < NUM_BLOCKS; i++) {
				for (int j = 0; j < NUM_BLOCKS; j++) {
					double neg = counter[i] * counter[j];			// might be too large for almost all the block pairs 
					if (i == j) {
						neg -= counter[i];
					}
					eta[i][j] = (m[i][j]+1)/(m[i][j]+neg+2);
				}
			}

		}
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

		init(args);
		train();

		return;
	}
}
