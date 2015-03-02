/**
	Implementation of the model in Discovering Latent Classes in Relational Data (Kemp 2004)
	Every node has a hard class label and the probability that two nodes link depends solely on the blocks they belong to
	Here we do not use sampling to update parameters; instead we use EM.
**/

import java.util.*;

public class BlockModel
{
//	public static ArrayList<
	public final static int NUM_BLOCKS = 2;
	public final static int MAX_ITER = 10;
	public static int NUM_NODES;
	public static String fileDir, dictDir;

	public static Map<String, Integer> idMap;						// RawID -> newID 
	public static double[][] data;								// network (N*N) 
	public static double[][] eta;								// block*block (K*K)
	public static int[] z;									// block assignment (N*1): #Node -> #Block 

	public static Scanner sc;

	/// initialize parameters 
	public static void
	init(
		String[] args
	) {
		fileDir = args[0];
		dictDir = args[1];

		idMap = new HashMap<String, Integer>();
		eta = new double[NUM_BLOCKS][NUM_BLOCKS];
		NUM_NODES = FileParser.readDict(dictDir, idMap, "\t");
		data = new double[NUM_NODES][NUM_NODES];
		FileParser.readData(fileDir, idMap, data);
		z = new int[NUM_NODES];

		sc = new Scanner(System.in);

		// random initialization - z 
		Random rand = new Random(0);
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
				double maxChange = -Double.MAX_VALUE, curChange;
				int bestK = -1, preK = z[n];
				for (int k = 0; k < NUM_BLOCKS; k++) {
					curChange = Evaluation.changeInObj(data, eta, z, n, k);
					if (curChange > maxChange) {
						bestK = k;
						maxChange = curChange;
					}						
				}
				z[n] = bestK;
				flag = flag && (bestK == preK);
			}

			if (flag) {
				System.out.println("z has converged.");
				break;
			}

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

					System.out.println("i = " + i + ", j = " + j + ", m = " + m[i][j] + ", neg = " + neg);
				}
			}
			int gu = sc.nextInt();

			// output z
			for (int i = 0; i < counter.length; i++) {
				System.out.printf("%d\t", (int)counter[i]);
			}
			System.out.printf("\n");

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

		FileParser.output("./res_z", z);
		FileParser.output("./res_eta", eta);

		return;
	}
}
