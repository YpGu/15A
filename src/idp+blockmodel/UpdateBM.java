/**
	UpdateBM.java: update one round of parameters of blockmodel

	Input:
		Data Matrix
		(Current) Block Assignment for Nodes
		(Current) Block Matrix
	Output:
		(Updated) Block Assignment for Nodes
		(Updated) Block Matrix
**/

import java.util.*;

public class UpdateBM
{
	// check blocks
	public static void
	checkBlocks(Map<String, Integer> z, int NUM_BLOCKS) {
		int[] counter = new int[NUM_BLOCKS];
		for (int k = 0; k < counter.length; k++) {
			counter[k] = 0;
		}
		for (Map.Entry<String, Integer> i: z.entrySet()) {
			int block = i.getValue();
			counter[block] += 1;
		}
		for (int i = 0; i < counter.length; i++) {
			System.out.printf("\t%d", counter[i]);
		}
		System.out.printf("\n");
	}


	// check empty blocks (for [current] block): break if all blocks are non-empty 
	public static boolean
	existEmptyBlock(Map<String, Integer> z, int NUM_BLOCKS, int init) {
		double[] tmpCounter = new double[NUM_BLOCKS];
		for (int k = 0; k < NUM_BLOCKS; k++) {
			tmpCounter[k] = 0;
		}
		for (Map.Entry<String, Integer> i: z.entrySet()) {
			int block = i.getValue();
			tmpCounter[block] += 1;
		}
		double zeroFlag = 1;
		for (int k = 0; k < NUM_BLOCKS; k++) {
			zeroFlag *= tmpCounter[k];
		}
		if (zeroFlag != 0) {
			return false;
		}

		return true;
	}


	/// Hard BM step 1 - estimate posterior z 
	public static boolean
	updateLatentHard(SparseMatrix trainData, Map<String, Integer> z, double[][] eta) {

		int NUM_BLOCKS = eta.length;
		boolean flag = true;						// convergence checker 

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

		checkBlocks(z, NUM_BLOCKS);
		System.out.println("\tObj = " + Evaluation.calcObj(trainData, eta, z));

		return false;
	}


	/// Hard BM step 2 - update mode parameters eta 
	public static void 
	updateParamHard(SparseMatrix trainData, Map<String, Integer> z, double[][] eta) {

		int NUM_BLOCKS = eta.length;
		double[][] m = new double[NUM_BLOCKS][NUM_BLOCKS];
		for (Map.Entry<Tuple<String, String>, Double> e: trainData.getMat().entrySet()) {
			String x = e.getKey().getX();
			String y = e.getKey().getY();
			double v = e.getValue();
			int zx = z.get(x);
			int zy = z.get(y);
			m[zx][zy] += v;
		}

		double[] counter = new double[NUM_BLOCKS];
		for (int k = 0; k < counter.length; k++) {
			counter[k] = 0;
		}
		for (Map.Entry<String, Integer> i: z.entrySet()) {
			int block = i.getValue();
			counter[block] += 1;
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
			}
		}

		System.out.println("\tObj = " + Evaluation.calcObj(trainData, eta, z));
	}

	
	/// Hard BM update: merge step 1 and 2 together 
	public static boolean
	updateHard(SparseMatrix trainData, Map<String, Integer> z, double[][] eta) {
		boolean res = updateLatentHard(trainData, z, eta);
		updateParamHard(trainData, z, eta);
		return res;
	}
}

