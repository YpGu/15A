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
	bkgUpdateLatentHard(SparseMatrix trainData, Map<String, Integer> z, double[][] eta) {

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
	bkgUpdateParamHard(SparseMatrix trainData, Map<String, Integer> z, double[][] eta) {

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
	bkgUpdateHard(SparseMatrix trainData, Map<String, Integer> z, double[][] eta) {
		boolean res = updateLatentHard(trainData, z, eta);
		bkgUpdateParamHard(trainData, z, eta);
		return res;
	}


	/// IDP update: using gradient descent 
	public static boolean
	idpUpdate(
		SparseMatrix data, SparseMatrix nData,
		Map<String, Double> vOut, Map<String, Double> vIn, Map<String, Double> vBias
		Map<String, Double> pi, double[][] gamma,						// gamma is the estimated weight before IDP mixture 
		double c,										// sample weight 
		double reg										// regularization coefficient 
	) {
		Map<String, Double> gradOut = new HashMap<String, Double>();
		Map<String, Double> gradIn = new HashMap<String, Double>();
		Map<String, Double> gradBias = new HashMap<String, Double>();

		for (String x: data.getDict()) {
			Set<String> s1 = data.getRow(x);
			for (String y: s1) {								// x -> y
				double p = logis(vOut.get(zx) * vIn.get(zy) + vBias.get(zy));
				double grad = gamma[x][y] * (1-p);
				gradOut.put(x, gradOut.get(x) + vIn.get(y) * grad);
				gradIn.put(y, gradIn.get(y) + vOut.get(x) * grad);
				gradBias.put(y, gradBias.get(y) + grad);
			}
		}
		for (String x: data.getDict()) {
			Set<String> s2 = nData.getRow(x);
			for (String y: s2) {								// x !-> y
				double p = logis(vOut.get(zx) * vIn.get(zy) + vBias.get(zy));
				double grad = gamma[x][y] * p;
				gradOut.put(x, gradOut.get(x) - vIn.get(y) * grad * c);
				gradIn.put(y, gradIn.get(y) - vOut.get(x) * grad * c);
				gradBias.put(y, gradBias.get(y) - grad * c);
			}
		}

		// regularizations
		for (String x: data.getDict()) {
			gradOut.put(x, gradOut.get(x) - reg * vOut.get(x));
			gradIn.put(x, gradIn.get(x) - reg * vIn.get(x));
		}

/*		// line search
		Map<String, Double> tmpOut = new HashMap<String, Double>();
		Map<String, Double> tmpIn = new HashMap<String, Double>();
		Map<String, Double> tmpBias = new HashMap<String, Double>();
		int lsIter = 0;
		double tmplr = lr;
		do {
			for (int x = 0; x < N; x++) {
				tmpOut.put(x, vOut.get(x) + tmplr * gradOut.get(x));
				tmpIn.put(x, vIn.get(x) + tmplr * gradIn.get(x));
				tmpBias.put(x, vBias.get(x) + tmplr * gradBias.get(x));
			}
			newObj = calObj(tmpAlpha, tmpBeta, tmpOut, tmpIn, tmpBias);
			tmplr *= 0.5;
			lsIter++;
		}
		while (newObj <= oldObj && lsIter < 10 && numIter > 0);
*/
		for (int x = 0; x < N; x++) {
			vOut.put(x, vOut.get(x) + tmplr * gradOut.get(x));
			vIn.put(x, vIn.get(x) + tmplr * gradIn.get(x));
			vBias.put(x, vBias.get(x) + tmplr * gradBias.get(x));
		}
}

