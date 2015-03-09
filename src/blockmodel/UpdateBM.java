/**
	UpdateBM.java: update one round of parameters of blockmodel
**/

import java.util.*;

public class UpdateBM
{
	// check blocks
	public static void
	checkBlocks(Map<String, Integer> z, double[] counter) {
		for (Map.Entry<String, Integer> i: z.entrySet()) {
			int block = i.getValue();
			counter[block] += 1;
		}
		for (int i = 0; i < counter.length; i++) {
			System.out.printf("\t%d", (int)counter[i]);
		}
		System.out.printf("\n");
	}


	// update: 1. estimate posterior z; 2. update model parameters eta 
	public static boolean
	update(SparseMatrix trainData, Map<String, Integer> z, double[][] eta) {

		int NUM_BLOCKS = eta.length;
		boolean flag = true;						// convergence checker 

		// Step 1: estimate posterior z 
		long sTime = System.currentTimeMillis();

		for (String s: trainData.getDict()) {
			int bestK = z.get(s), preK = z.get(s);
			double maxChange = -Double.MIN_VALUE, curChange;

			for (int k = 0; k < NUM_BLOCKS; k++) {
				if (k != preK) {
					curChange = Evaluation.changeInObj(trainData, eta, z, s, k);
					if (curChange > maxChange) {
						bestK = k;
						maxChange = curChange;
					}
				}
			}
			z.put(s, bestK);
			flag = flag && (bestK == preK);
		}

		if (flag) {
			System.out.println("\tz has converged.");
			return true;
		}

		long fTime = System.currentTimeMillis();
//		System.out.println("Time = " + (fTime-sTime));

		double[] counter = new double[NUM_BLOCKS];			// counter (K*1): #Block -> num of nodes 
		checkBlocks(z, counter);
		System.out.println("\tObj = " + Evaluation.calcObj(trainData, eta, z));

		// Step 2: update model parameters eta 
		double[][] m = new double[NUM_BLOCKS][NUM_BLOCKS];

		for (Map.Entry<Tuple<String, String>, Double> e: trainData.getMat().entrySet()) {
			String x = e.getKey().getX();
			String y = e.getKey().getY();
			double v = e.getValue();
			int zx = z.get(x);
			int zy = z.get(y);
			m[zx][zy] += 1;
		}

		for (int i = 0; i < NUM_BLOCKS; i++) {
			for (int j = 0; j < NUM_BLOCKS; j++) {
				double pb = counter[i] * counter[j];
				if (i == j) {
					pb -= counter[i];
				}
				if (pb != 0) {
					eta[i][j] = (m[i][j]+1)/(pb+2);
				}
				else {
					eta[i][j] = 0;
				}

//				System.out.println("i = " + i + ", j = " + j + ", eta = " + eta[i][j]);
			}
		}

		System.out.println("\tObj = " + Evaluation.calcObj(trainData, eta, z));

		return false;
	}
}

