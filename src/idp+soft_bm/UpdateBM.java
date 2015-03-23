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
	// init parameters 
	public static void 
	updateParamAtInit(
		SparseMatrix posData, 
		SparseMatrix negData, 
		Map<String, double[]> theta, 
		double[][] eta,
		double sw
	) {
		int K = eta.length;

		// init theta 
		double[] thetaInit = new double[K];
		for (int k = 0; k < K; k++) {
			thetaInit[k] = 1/(double)K;
		}
		for (String x: posData.getDict()) {
			theta.put(x, thetaInit);
		}

		// update eta 
		double[][] posM = new double[K][K];
		double[][] negM = new double[K][K];
		for (String i: posData.getDict()) {
			Set<String> s1 = posData.getRow(i);
			for (String j: s1) {
				for (int g = 0; g < K; g++) 
					for (int h = 0; h < K; h++) 
						posM[g][h] += theta.get(i)[g] * theta.get(j)[h];
			}
			Set<String> s2 = posData.getRowComplement(i);
			for (String j: s2) {
				for (int g = 0; g < K; g++) 
					for (int h = 0; h < K; h++) 
						negM[g][h] += theta.get(i)[g] * theta.get(j)[h];
			}
		}

		for (int g = 0; g < K; g++) {
			for (int h = 0; h < K; h++) {
				if (posM[g][h] + negM[g][h] != 0) {
					eta[g][h] = posM[g][h] / (posM[g][h] + negM[g][h]);
				}
				else {
					eta[g][h] = 0;
				}
			}
		}

		return;
	}


	/// Soft BM step 1 - estimate posterior gamma and update theta 
	public static void
	updateLatentSoft(
		SparseMatrix posData, 
		SparseMatrix negData, 
		Map<String, double[]> theta, 
		Map<String, Map<String, Double>> piGamma,
		double[][] eta, 
		double sw
	) {
//		long sTime = System.currentTimeMillis();
		int K = eta.length;
		// this gamma stores the posterior for node pair (i,j) belongs to block pair (g,h)
		Map<String, Map<String, double[][]>> gamma = new HashMap<String, Map<String, double[][]>>();

		// calculate gamma (E-step)
		for (String i: posData.getDict()) {
			Map<String, double[][]> js = new HashMap<String, double[][]>();
			for (String j: posData.getRow(i)) {
				double normGamma = 0;
				double[][] localGamma = new double[K][K];
				for (int g = 0; g < K; g++) {
					for (int h = 0; h < K; h++) {
						localGamma[g][h] = theta.get(i)[g] * theta.get(j)[h] * eta[g][h];
						normGamma += localGamma[g][h];
					}
				}
				if (normGamma != 0) 
					for (int g = 0; g < K; g++) 
						for (int h = 0; h < K; h++) 
							localGamma[g][h] /= normGamma;
				js.put(j, localGamma);
			}
			for (String j: posData.getRowComplement(i)) {
				double normGamma = 0;
				double[][] localGamma = new double[K][K];
				for (int g = 0; g < K; g++) {
					for (int h = 0; h < K; h++) {
						localGamma[g][h] = theta.get(i)[g] * theta.get(j)[h] * (1-eta[g][h]);
						normGamma += localGamma[g][h];
					}
				}
				if (normGamma != 0) 
					for (int g = 0; g < K; g++) 
						for (int h = 0; h < K; h++) 
							localGamma[g][h] /= normGamma;
				js.put(j, localGamma);
			}
			gamma.put(i,js);
		}

		// update theta^{new} (one of M-step)
		for (String i: posData.getDict()) {
			double normTheta = 0;
			double[] localTheta = new double[K];
			for (int g = 0; g < K; g++) {
				for (String j: posData.getRow(i)) {
					double w = 1 - piGamma.get(i).get(j);
					for (int h = 0; h < K; h++) {
						localTheta[g] += w * gamma.get(i).get(j)[g][h];
						normTheta += localTheta[g];
					}
				}
				for (String j: posData.getRowComplement(i)) {
					double w = 1 - piGamma.get(i).get(j);
					for (int h = 0; h < K; h++) {
						localTheta[g] += w * gamma.get(i).get(j)[g][h];
						normTheta += localTheta[g];
					}
				}
	
				for (String j: posData.getColumn(i)) {
					double w = 1 - piGamma.get(i).get(j);
					for (int h = 0; h < K; h++) {
						localTheta[g] += w * gamma.get(j).get(i)[h][g];
						normTheta += localTheta[g];
					}
				}
				for (String j: posData.getColumnComplement(i)) {
					double w = 1 - piGamma.get(i).get(j);
					for (int h = 0; h < K; h++) {
						localTheta[g] += w * gamma.get(j).get(i)[h][g];
						normTheta += localTheta[g];
					}
				}
			}
			if (normTheta != 0) {
				for (int g = 0; g < K; g++) {
					localTheta[g] /= normTheta;
				}
			}
			theta.put(i, localTheta);
		}


//		long fTime = System.currentTimeMillis();
//		System.out.println("Time = " + (fTime-sTime));

		return;
	}


	/// Soft BM step 2 - update model parameters eta 
	public static void 
	updateParamSoft(
		SparseMatrix posData, 
		SparseMatrix negData, 
		Map<String, double[]> theta,
		Map<String, Map<String, Double>> piGamma,
		double[][] eta,
		double sw
	) {
		int K = eta.length;

		// update eta 
		double[][] posM = new double[K][K];
		double[][] negM = new double[K][K];
		for (String i: posData.getDict()) {
			Set<String> s1 = posData.getRow(i);
			for (String j: s1) {
				double w = 1 - piGamma.get(i).get(j);
				for (int g = 0; g < K; g++) 
					for (int h = 0; h < K; h++) 
						posM[g][h] += w * theta.get(i)[g] * theta.get(j)[h];
			}
			Set<String> s2 = posData.getRowComplement(i);
			for (String j: s2) {
				double w = 1 - piGamma.get(i).get(j);
				for (int g = 0; g < K; g++) 
					for (int h = 0; h < K; h++) 
						negM[g][h] += w * theta.get(i)[g] * theta.get(j)[h];
			}
		}

		for (int g = 0; g < K; g++) {
			for (int h = 0; h < K; h++) {
				if (posM[g][h] + negM[g][h] != 0) {
					eta[g][h] = posM[g][h] / (posM[g][h] + negM[g][h]);
				}
				else {
					eta[g][h] = 0;
				}
			}
		}

		return;
	}

	
	/// Hard BM update: merge step 1 and 2 together 
	public static void
	updateSoft(SparseMatrix posData, SparseMatrix negData, Map<String, double[]> theta, double[][] eta, Map<String, Map<String, Double>> piGamma, double sw,
			Map<String, Double> vOut, Map<String, Double> vIn, Map<String, Double> vBias, Map<String, Double> pi, double reg, boolean calc) {

		updateLatentSoft(posData, negData, theta, piGamma, eta, sw);
		if (calc) {
			double obj = Evaluation.calcObj(posData, negData, theta, eta, vOut, vIn, vBias, pi, sw, reg);
			System.out.println("\t\tObjective function (after optimizing z) = " + obj);
		}
		updateParamSoft(posData, negData, theta, piGamma, eta, sw);

		return;
	}
}

