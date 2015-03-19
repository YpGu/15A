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

	/// init
	public static void 
	updateParamAtInit(
		SparseMatrix posData, 
		Map<String, double[]> gamma,
		double[][] matB
	) {
		int NUM_BLOCKS = matB.length;
		int N = posData.getDict().size();

		// init gamma
		Random rand = new Random(0);
		for (String x: posData.getDict()) {
			double[] tmpV = new double[NUM_BLOCKS+1];
			for (int k = 0; k < NUM_BLOCKS; k++) {
	//			tmpV[k] = 1.0/(double)NUM_BLOCKS;
				tmpV[k] = rand.nextDouble();
				tmpV[NUM_BLOCKS] += tmpV[k];
			}
			if (tmpV[NUM_BLOCKS] != 0) {
				for (int k = 0; k < NUM_BLOCKS; k++) {
					tmpV[k] /= tmpV[NUM_BLOCKS];
				}
			}
			gamma.put(x, tmpV);
		}
	// output some gamma
		int ng = 0;
		for (Map.Entry<String, double[]> e: gamma.entrySet()) {
			if (ng == 3) break;
			double[] v = e.getValue();
			for (int i = 0; i < v.length-1; i++) {
				System.out.printf("%f\t", v[i]);
			}
			System.out.printf("\n");
			ng += 1;
		}

		// init B
		double[][] m = new double[NUM_BLOCKS][NUM_BLOCKS];
		double[][] mBar = new double[NUM_BLOCKS][NUM_BLOCKS];
		for (String x: posData.getDict()) {
			Set<String> s1 = posData.getRow(x);
			for (String y: s1) {
				for (int g = 0; g < NUM_BLOCKS; g++) {
					for (int h = 0; h < NUM_BLOCKS; h++) {
						m[g][h] += 1;
					}
				}
			}
			Set<String> s2 = posData.getRowComplement(x);
			for (String y: s2) {
				for (int g = 0; g < NUM_BLOCKS; g++) {
					for (int h = 0; h < NUM_BLOCKS; h++) {
						mBar[g][h] += 1;
					}
				}
			}
		}
		double c = m[0][0]/(m[0][0]+mBar[0][0]);
		for (int i = 0; i < NUM_BLOCKS; i++) {
			for (int j = 0; j < NUM_BLOCKS; j++) {
				matB[i][j] = c;
			}
		}

		return;
	}


	/// update \phi 
	public static void
	update1(
		SparseMatrix posData,
		double[][] matB,
		Map<String, double[]> gamma,
		Map<String, Map<String, double[]>> phiP2Q,
		Map<String, Map<String, double[]>> phiQ2P
	) {
		int NUM_BLOCKS = matB.length;

		for (String p: posData.getDict()) {
			Map<String, double[]> qs = new HashMap<String, double[]>();

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
				qs.put(q, tmpV);
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
				qs.put(q, tmpV);
			}

			phiP2Q.put(p, qs);
		}

		for (String q: posData.getDict()) {
			Map<String, double[]> ps = new HashMap<String, double[]>();

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
				ps.put(p, tmpV);
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
				ps.put(p, tmpV);
			}

			phiQ2P.put(q, ps);
		}

		return;
	}


	/// update B and gamma 
	public static void
	update2(
		SparseMatrix posData,
		Map<String, Map<String, double[]>> phiP2Q,
		Map<String, Map<String, double[]>> phiQ2P,
		Map<String, double[]> gamma,
		double[][] matB
	) {
		// update B
		int l = matB.length;
		for (int g = 0; g < l; g++) {
			for (int h = 0; h < l; h++) {
				matB[g][h] = 0;
			}
		}
		double[][] matBDeno = new double[l][l];
		for (String p: posData.getDict()) {
			for (String q: posData.getRow(p)) {
//				double v = posData.getElement(p,q);
				double v = 1;
				for (int g = 0; g < l; g++) {
					for (int h = 0; h < l; h++) {
						double c = phiP2Q.get(p).get(q)[g] * phiQ2P.get(q).get(p)[h];
						matB[g][h] += v * c;
						matBDeno[g][h] += c;
					}
				}
			}
			for (String q: posData.getRowComplement(p)) {
				for (int g = 0; g < l; g++) {
					for (int h = 0; h < l; h++) {
						double c = phiP2Q.get(p).get(q)[g] * phiQ2P.get(q).get(p)[h];
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

	// try outputing B
		for (int g = 0; g < l; g++) {
			for (int h = 0; h < l; h++) {
				System.out.printf("%f\t", matB[g][h]);
			}
			System.out.printf("\n");
		}

		// update gamma
		for (String p: posData.getDict()) {
			double norm = 0;
			double[] tmpV = new double[l];
			for (int k = 0; k < l; k++) {
				for (String q: posData.getRow(p)) {
					tmpV[k] += phiP2Q.get(p).get(q)[k];
					tmpV[k] += phiQ2P.get(q).get(p)[k];
				}
				norm += tmpV[k];
			}
			if (norm != 0) {
				for (int k = 0; k < l; k++) {
					tmpV[k] /= norm;
				}
			}
		}

		return;
	}


	/// Soft BM update: merge step 1 and 2 together 
	public static boolean
	updateLatentSoft(
		SparseMatrix posData, 
		Map<String, Map<String, double[]>> phiP2Q,
		Map<String, Map<String, double[]>> phiQ2P,
		double[][] matB,
		Map<String, double[]> gamma
	) {
		System.out.println("\tupdate1");
		update1(posData, matB, gamma, phiP2Q, phiQ2P);
		System.out.println("\tupdate2");
		update2(posData, phiP2Q, phiQ2P, gamma, matB);

		return false;
	}
	

}

