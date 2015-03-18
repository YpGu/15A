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
	// check blocks: return True if no block is empty 
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

	/// start working TODO 
	/// update \phi 
	public static void
	update1(
		SparseMatrix posData,
		double[][] matB,
		Map<String, Double[]> gamma,
		Map<String, Map<String, Double[]>> phiP2Q,
		Map<String, Map<String, Double[]>> phiQ2P,
	) {
		for (String p: posData.getDict()) {
			Map<String, Double[]> qs = new HashMap<String, Double[]>();

			Set<String> s1 = posData.getRow(p);
			for (String q: s1) {
				double norm = 0;
				double[] tmpV = new double[NUM_BLOCKS];
				for (int g = 0; g < NUM_BLOCKS; g++) {
					double v = 1;
					for (int h = 0; h < NUM_BLOCKS; h++) {
						double power = phiQ2P.get(q).get(p)[h];
						v *= Math.pow(matB[g][h], power);
					}
					v *= Math.exp(Evaluation.expt(gamma, p, g));

					tmpV[g] = v;
					norm += v;
				}
				if (norm != 0) {
					for (int g = 0; g < NUM_BLOCKS; g++) {
						tmpV[g] /= norm;
					}
				}
				qs.set(q, tmpV);
			}
			Set<String> s2 = posData.getRowComplement(p);
			for (String q: s2) {
				double norm = 0;
				double[] tmpV = new double[NUM_BLOCKS];
				for (int g = 0; g < NUM_BLOCKS; g++) {
					double v = 1;
					for (int h = 0; h < NUM_BLOCKS; h++) {
						double power = phiQ2P.get(q).get(p)[h];
						v *= Math.pow(1-matB[g][h], power);
					}
					v *= Math.exp(Evaluation.expt(gamma, p, g));

					tmpV[g] = v;
					norm += v;
				}
				if (norm != 0) {
					for (int g = 0; g < NUM_BLOCKS; g++) {
						tmpV[g] /= norm;
					}
				}
				qs.set(q, tmpV);
			}

			phiP2Q.set(p, qs);
		}

		for (String q: posData.getDict()) {
			Map<String, Double[]> ps = new HashMap<String, Double[]>();

			Set<String> s3 = posData.getColumn(q);
			for (String p: s3) {
				double norm = 0;
				double[] tmpV = new double[NUM_BLOCKS];
				for (int h = 0; h < NUM_BLOCKS; h++) {
					double v = 1;
					for (int g = 0; g < NUM_BLOCKS; g++) {
						double power = phiP2Q.get(p).get(q)[g];
						v *= Math.pow(matB[g][h], power);
					}
					v *= Math.exp(Evaluation.expt(gamma, q, h));

					tmpV[h] = v;
					norm += v;
				}
				if (norm != 0) {
					for (int g = 0; g < NUM_BLOCKS; g++) {
						tmpV[g] /= norm;
					}
				}
				ps.set(p, tmpV);
			}

			Set<String> s4 = posData.getColumnComplement(q);
			for (String p: s4) {
				double norm = 0;
				double[] tmpV = new double[NUM_BLOCKS];
				for (int h = 0; h < NUM_BLOCKS; h++) {
					double v = 1;
					for (int g = 0; g < NUM_BLOCKS; g++) {
						double power = phiP2Q.get(p).get(q)[g];
						v *= Math.pow(1-matB[g][h], power);
					}
					v *= Math.exp(Evaluation.expt(gamma, q, h));

					tmpV[h] = v;
					norm += v;
				}
				if (norm != 0) {
					for (int g = 0; g < NUM_BLOCKS; g++) {
						tmpV[g] /= norm;
					}
				}
				ps.set(p, tmpV);
			}

			phiQ2P.set(q, ps);
		}

		return;
	}


	/// update B and gamma 
	public static void
	update2(
		Map<String, Map<String, Double[]>> phiP2Q,
		Map<String, Map<String, Double[]>> phiQ2P,
		Map<String, Double[]> gamma,
		double[][] matB
	) {
		int l = matB.length;
		for (int g = 0; g < l; g++) {
			for (int h = 0; h < l; h++) {
				matB[g][h] = 0;
			}
		}
		double[][] matBDeno = new double[l][l];

		for (String p: posData.getDict()) {
			for (String q: posData.getRow(p)) {
				double v = posData.getElement(p,q);
				for (int g = 0; g < l; g++) {
					for (int h = 0; h < l; h++) {
						double c = phiP2Q.get(p).get(q)[g] * phiQ2P.get(q).get(p)[h];
						matB[g][h] += v * c;
						matBDeno[g][h] += c;
					}
				}
			}
		}

		for (int g = 0; g < l; g++) {
			for (int h = 0; h < l; h++) {
				if (matBDeno[g][h] != 0) {
					matB[g][h] /= matBDeno[g][h];
				}
			}
		}

		return;
	}


	/// Soft BM step 1 - estimate posterior z 
	public static boolean
	updateLatentSoft(
		SparseMatrix posData, 
		SparseMatrix negData, 
		Map<String, Map<String, Double[]>> phiP2Q,
		Map<String, Map<String, Double[]>> phiQ2P,
		double[][] matB,
		Map<String, Double[]> gamma,
		double sw
	) {
		update1(phiP2Q, phiQ2P);
		update2(gamma, matB);

		return false;
	}


	/// Hard BM step 2 - update model parameters eta 
	public static void 
	updateParamAtInit(
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

