/**
	MixtureUpdate.java: provide the variational E/M step for multinomial mixture training
**/

import java.util.*;

public class MixtureUpdate
{
	// Learning Rate
	public static final double lr = 0.001;

	// Estimate Variational Parameters using EM 
	public static void
	variationEM(
		SparseMatrix<Integer> trainData,
		int userIndex,								// maximize the likelihood for v_{userIndex} 
		double[] alpha,	double[][] beta, double[] pi,
		double[] p, double[] q, double[] b,
		double[] gamma, Map<Integer, double[][]> allPhi, double[] varphi
	) {
		int K = alpha.length, N = p.length, M = trainData.getRow(userIndex).size();
		double[][] phi = allPhi.get(userIndex);

//		double l = Evaluation.calcLikelihood(trainData, alpha, beta, pi, p, q, b, gamma, allPhi, varphi);
//		System.out.println("=== Likelihood: " + l + " ===");
	
		// Update Bernoulli Parameter \varphi |Time: O(KM)|
		while (true) {
			int i = userIndex;
			if (pi[i] == 1) {
				varphi[i] = 1;
				break;
			}
			else if (pi[i] == 0) {
				varphi[i] = 0;
				break;
			}
			double di = pi[i] / (1-pi[i] + Double.MIN_VALUE) + Double.MIN_VALUE;
			double ss = Evaluation.sumSigma(i, p, q, b);
			int m = 0;
			for (int j: trainData.getRow(i)) {
				double v = Math.exp(p[i] * q[j] + b[j]) / ss;
				for (int k = 0; k < K; k++) 
					v -= phi[m][k] * Math.log(beta[k][j] + Double.MIN_VALUE);
				di *= trainData.get(i, j);
				m += 1;
			}
			di *= Math.log(pi[i] / (1-pi[i]));
			di = Evaluation.logis(di);
			varphi[i] = di;

			break;
		}

		// Update Dirichlet Parameter \gamma |Time: O(KM)|
		if (true) {
			double[] resGamma = new double[K];
			double normGamma = 0;
			for (int k = 0; k < K; k++) {
				resGamma[k] = alpha[k];
				normGamma += alpha[k];
				for (int m = 0; m < M; m++) {
					resGamma[k] += phi[m][k];
					normGamma += phi[m][k];
				}
			}
			if (normGamma == 0) normGamma = 1;
			for (int k = 0; k < K; k++)
				resGamma[k] /= normGamma;

			for (int k = 0; k < K; k++) {
				gamma[k] = resGamma[k];
			}
 		}

		// Update Multinomial Parameter \phi |Time: O(KM)|
		if (true) {
			double[][] resPhi = new double[M][K];
			int i = userIndex, m = 0;
			double normPhi = 0;
			for (int j: trainData.getRow(i)) {
				for (int k = 0; k < K; k++) {				// note: phi is in log-scale 
					resPhi[m][k] = (1-varphi[i]) * trainData.get(i, j) * Math.log(beta[k][j] + Double.MIN_VALUE);
					if (k > 0) 
						normPhi = Evaluation.logSum(normPhi, resPhi[m][k]);
					else 
						normPhi = resPhi[m][k];
				}
				for (int k = 0; k < K; k++)				// Normalize
					resPhi[m][k] = Math.exp(resPhi[m][k] - normPhi);
				m += 1;
			}
			allPhi.put(i, resPhi);
		}

		return;
	}

	// Update Model Parameters using Gradient Ascent (after Variational Inference for all users) 
	public static void
	modelParamEstimate(
		SparseMatrix<Integer> trainData,
		double[] alpha,	double[][] beta, double[] pi,
		double[] p, double[] q, double[] b,
		double[] gamma, Map<Integer, double[][]> allPhi, double[] varphi
	) {
		int K = alpha.length, N = p.length;

		// Update \beta - O(KE) 
		double[][] tmpBeta = new double[K][N];
		for (int k = 0; k < K; k++) {
			double normBeta = 0;
			for (int i = 0; i < N; i++) {
				double[][] phi = allPhi.get(i);
				int m = 0;
				for (int j: trainData.getRow(i)) {
					double v = trainData.get(i, j) * phi[m][k] * (1-varphi[i]);
					tmpBeta[k][j] += v;
					normBeta += v;
					m += 1;
				}
			}
			if (normBeta == 0) normBeta = 1;
			for (int j = 0; j < N; j++) 
				tmpBeta[k][j] /= normBeta;
		}

		for (int k = 0; k < K; k++) 
			for (int j = 0; j < N; j++)
				beta[k][j] = tmpBeta[k][j];

		// Update \pi
		for (int k = 0; k < K; k++)
			pi[k] = varphi[k];

		// Update p, q, b - O(E) 
		for (int iter = 0; iter < 10; iter++) {
			System.out.println("*** Updating pqb " + iter + " ***");
			double[] gradP = new double[N], gradQ = new double[N], gradB = new double[N];
			for (int i = 0; i < N; i++) {
				if (varphi[i] != 0) {
					double ssw = Evaluation.sumSigmaWeighted(i, p, q, b);
					double ss = Evaluation.sumSigma(i, p, q, b);
					for (int j: trainData.getRow(i)) {
						double sg = Math.exp(p[i] * q[j] + b[j]) / ss;		// \sigma_{ij} 
						double v = trainData.get(i, j) * varphi[i];
						gradP[i] += v * (q[j] - ssw/ss);			// \sum_j { \varphi_i * n_{ij} * (q_j - weightedSum)
						gradQ[j] += v * p[i] * (1-sg);				// \sum_i { \varphi_i * n_{ij} * (1-\sigma_{ij}) * p_i 
//						gradB[j] += v * (1-sg);					// \sum_i { \varphi_i * n_{ij} * (1-\sigma_{ij}) 
					}
				}
			}

			// Update All
			for (int i = 0; i < N; i++) {
				p[i] += lr * gradP[i];
				q[i] += lr * gradQ[i];
//				b[i] += lr * gradB[i];
			}
			double l = 0;
			for (int i = 0; i < trainData.getXDictSize(); i++) {
				l += Evaluation.calcLikelihood(trainData, i, alpha, beta, pi, p, q, b, gamma, allPhi, varphi);
			}
			System.out.println("=== Likelihood: " + l + " ===");
		}

		return;
	}

	public static void
	update (
		SparseMatrix<Integer> trainData,
		double[] alpha,	double[][] beta, double[] pi,
		double[] p, double[] q, double[] b,
		double[] gamma, Map<Integer, double[][]> allPhi, double[] varphi
	) {
		for (int i = 0; i < trainData.getXDictSize(); i++) {
			if (i%500 == 0) System.out.println("  user " + i);
			variationEM(trainData, i, alpha, beta, pi, p, q, b, gamma, allPhi, varphi);
		}
		modelParamEstimate(trainData, alpha, beta, pi, p, q, b, gamma, allPhi, varphi);

		System.out.println("Calculating likelihood...");
		double l = 0;
		for (int i = 0; i < trainData.getXDictSize(); i++) {
			l += Evaluation.calcLikelihood(trainData, i, alpha, beta, pi, p, q, b, gamma, allPhi, varphi);
		}
		System.out.println("=== Likelihood: " + l + " ===");

		return;
	}
}
