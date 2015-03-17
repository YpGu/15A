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
	public static boolean
	checkBlocks(Map<String, Integer> z, int NUM_BLOCKS) {
		boolean noZero = true;
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
			noZero = noZero && (counter[i] != 0);
		}
		System.out.printf("\n");
		return noZero;
	}


	// check empty blocks (for [current] block): break if all blocks are non-empty 
	public static boolean
	existEmptyBlock(Map<String, Integer> z, int NUM_BLOCKS) {
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
	updateLatentHard(
		SparseMatrix posData, 
		SparseMatrix negData, 
		Map<String, Integer> z, 
		double[][] eta, 
		Map<String, Map<String, Double>> gamma,
		double sw
	) {
//		long sTime = System.currentTimeMillis();
		int NUM_BLOCKS = eta.length;
		boolean flag = true;						// convergence checker 

		for (String s: posData.getDict()) {
			int bestK = z.get(s), preK = z.get(s);
			double maxChange = -Double.MIN_VALUE, curChange;
			double maxObj = -Double.MAX_VALUE, curObj;

			for (int k = 0; k < NUM_BLOCKS; k++) {
				if (k != preK) {
					curChange = Evaluation.changeInObj(posData, eta, gamma, z, s, k, sw);
//					curChange = Evaluation.changeInObj(posData, eta, z, s, k);
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

//		long fTime = System.currentTimeMillis();
//		System.out.println("Time = " + (fTime-sTime));

		if (!checkBlocks(z, NUM_BLOCKS)) {
			return true;
		}

		return false;
	}


	/// Hard BM step 2 - update model parameters eta 
	public static void 
	updateParamHardAtInit(
		SparseMatrix posData, 
		SparseMatrix negData, 
		Map<String, Integer> z, 
		double[][] eta,
		double sw
	) {
		int NUM_BLOCKS = eta.length;

		double[][] m = new double[NUM_BLOCKS][NUM_BLOCKS];
		double[][] mBar = new double[NUM_BLOCKS][NUM_BLOCKS];

		for (String x: posData.getDict()) {
			Set<String> s1 = posData.getRow(x);
			for (String y: s1) {
				int zx = z.get(x);
				int zy = z.get(y);
				m[zx][zy] += 1;
			}
			Set<String> s2 = posData.getRowComplement(x);
			for (String y: s2) {
				int zx = z.get(x);
				int zy = z.get(y);
				mBar[zx][zy] += 1;
			}
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
				if (m[i][j] + mBar[i][j] != 0) {
					eta[i][j] = (m[i][j]+1)/(m[i][j]+mBar[i][j]+2);
				}
				else {
					eta[i][j] = 0;
				}
			}
		}

		return;
	}


	/// Hard BM step 2 - update model parameters eta 
	public static void 
	updateParamHard(
		SparseMatrix posData, 
		SparseMatrix negData, 
		Map<String, Integer> z, 
		double[][] eta,
		Map<String, Map<String, Double>> gamma,
		double sw
	) {
		int NUM_BLOCKS = eta.length;

		double[][] m = new double[NUM_BLOCKS][NUM_BLOCKS];
		double[][] mBar = new double[NUM_BLOCKS][NUM_BLOCKS];
		for (String x: posData.getDict()) {
			Set<String> s1 = posData.getRow(x);
			for (String y: s1) {
				double v = (1-gamma.get(x).get(y));
				int zx = z.get(x);
				int zy = z.get(y);
				m[zx][zy] += v;
			}
			Set<String> s2 = posData.getRowComplement(x);
			for (String y: s2) {
				double v = (1-gamma.get(x).get(y));
				int zx = z.get(x);
				int zy = z.get(y);
				mBar[zx][zy] += v;
			}
		}

		for (int i = 0; i < NUM_BLOCKS; i++) {
			for (int j = 0; j < NUM_BLOCKS; j++) {
				if (m[i][j] + mBar[i][j] != 0) {
					eta[i][j] = (m[i][j]+1)/(m[i][j]+mBar[i][j]+2);
				}
				else {
					eta[i][j] = 0;
				}
			}
		}

		return;
	}

	
	/// Hard BM update: merge step 1 and 2 together 
	public static boolean
	updateHard(SparseMatrix posData, SparseMatrix negData, Map<String, Integer> z, double[][] eta, Map<String, Map<String, Double>> gamma, double sw,
			Map<String, Double> vOut, Map<String, Double> vIn, Map<String, Double> vBias, Map<String, Double> pi, double reg, boolean calc) {
		boolean res = updateLatentHard(posData, negData, z, eta, gamma, sw);
		if (calc) {
			double obj = Evaluation.calcObj(posData, negData, eta, z, vOut, vIn, vBias, pi, sw, reg);
			System.out.println("\t\tObjective function (after optimizing z) = " + obj);
		}
		updateParamHard(posData, negData, z, eta, gamma, sw);
		return res;
	}
}

