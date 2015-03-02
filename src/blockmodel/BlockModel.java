/**
	Implementation of the model in Discovering Latent Classes in Relational Data (Kemp 2004)
	Here we do not use sampling to update parameters; instead we use EM.
**/

import java.util.*;

public class BlockModel
{
//	public static ArrayList<
	public final static int NUM_BLOCKS = 10;
	public final static int MAX_ITER = 1000;
	public static int NUM_NODES;

	public static double[][] data;								// network (N*N) 
	public static double[][] eta;								// block*block (K*K)
	public static int[] z;									// block assignment (N*1): #Node -> #Block 


	public static void
	init(
	) {
		eta = new double[NUM_BLOCKS][NUM_BLOCKS];

		return;
	}


	public static void
	train(
	) {
		for (int iter = 0; iter < MAX_ITER; iter++) {

			// check convergence
			boolean flag = true;

			// E-step 
			for (int n = 0; n < NUM_NODES; n++) {
				double maxObj = 0, curObj;
				int bestK = -1, preK = z[n];
				for (int k = 0; k < NUM_BLOCKS; k++) {
					z[n] = k;
					curObj = Evaluation.calcObj(data, eta, z);
					if (curObj > maxObj) {
						bestK = k;
						maxObj = curObj;
					}
				}
				z[n] = bestK;
				flag = flag || (bestK != preK);
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

	public static void
	main(
		String[] args
	) {
		return;
	}
}
