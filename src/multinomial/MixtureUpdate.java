/**
	MixtureUpdate.java: provide the variational E/M step for multinomial mixture training
**/

import java.util.*;

public class MixtureUpdate
{
	// Learning Rate
	public static final double lr = 0.001;

	// Handle Large/Small Numbers 
	public static void
	saveToPower(double[][] resPhi, int[][] resPhiPosNegPower, int i, int k) {
		while (resPhi[i][k] >= 10) {
			resPhi[i][k] *= 0.1;
			resPhiPosNegPower[i][k] += 1;
		}
		while (resPhi[i][k] < 1) {
			resPhi[i][k] *= 10;
			resPhiPosNegPower[i][k] -= 1;
		}
		return;
	}

	// Normalize Vector in Special Forms 
	public static void
	normalize(double[][] resPhi, int[][] resPhiPosNegPower, int i) {
		int K = resPhi[0].length;
		double norm = 0;
		for (int k = 0; k < K; k++) 
			norm += resPhi[i][k] * Math.pow(10, resPhiPosNegPower[i][k]);
		for (int k = 0; k < K; k++)
			resPhi[i][k] /= norm;
		return;
	}

	// Estimate Variational Parameters using EM 
	public static void
	variationEM(
		SparseMatrix<Integer> trainData,
		double[] alpha,	double[][] beta, double[] pi,
		double[] p, double[] q, double[] b,
		double[] gamma, double[][] phi, double[] varphi
	) {
		int K = alpha.length, N = p.length;

		// Update Dirichlet Parameter \gamma |Time: O(KV)|
		double[] resGamma = new double[K];
		double normGamma = 0;
		for (int k = 0; k < K; k++) {
			resGamma[k] = alpha[k];
			normGamma += alpha[k];
			for (int i = 0; i < N; i++) {
				resGamma[k] += phi[i][k];
				normGamma += phi[i][k];
			}
		}
		if (normGamma == 0) normGamma = 1;
		for (int k = 0; k < K; k++)
			resGamma[k] /= normGamma;

		// Update Multinomial Parameter \phi |Time: O(KE)|
/*		double[][] resPhi = new double[N][K];
		int[][] resPhiPosNegPower = new int[N][K];				// pos-neg: 10 ^ (pos-neg) 
		for (int i = 0; i < N; i++)
			for (int k = 0; k < K; k++) 
				resPhi[i][k] = 1;
		for (int i = 0; i < N; i++) {
			double eBase = Math.exp(-varphi[i]);
			for (int k = 0; k < K; k++) {
				double a = dLogGamma(gamma[k]);
				resPhi[i][k] *= Math.exp(a);				// exp[ \Psi(\gamma_k) - \Psi(sumOfGamma) ] 
				saveToPower(resPhi, resPhiPosNegPower, i, k);
				double ePower = 0;
				// TODO: select one j or many j`s?
				// ----- if all j`s: \sum_j {\beta_{kj}} is tooooooo small 
				for (int j: trainData.getRow(i)) {
					resPhi[i][k] *= Math.pow(1 / (beta[k][j] + Double.MIN_VALUE), trainData.get(i, j) * varphi[i]);
					if (Double.isInfinite(resPhi[i][k])) {
//					if (true) {
						System.out.println("beta = " + beta[k][j] + " nij = " + trainData.get(i,j) + " varphi = " + varphi[i]);
						System.out.println("base = " + (1/(beta[k][j]+Double.MIN_VALUE)));
						System.out.println("res = " + Math.pow(1 / (beta[k][j] + Double.MIN_VALUE), trainData.get(i, j) * varphi[i]));
						Scanner sc = new Scanner(System.in);
						int gu = sc.nextInt();
					}
					saveToPower(resPhi, resPhiPosNegPower, i, k);
//					ePower += trainData.get(i, j) * Math.log(beta[k][j] + Double.MIN_VALUE);
				}
//				resPhi[i][k] *= Math.pow(eBase, ePower);		// exp[ -\varphi_i ] ^ [ \sum_j {n_{ij} log \beta_{kj} ] 
			}
			double sumPhi = 0;						// normalization 
			normalize(resPhi, resPhiPosNegPower, i);
//			for (int k = 0; k < K; k++) sumPhi += resPhi[i][k];
//			if (sumPhi == 0) sumPhi = 1;
//			for (int k = 0; k < K; k++)
//				resPhi[i][k] /= sumPhi;
		}
*/
		// Update Multinomial Parameter \phi |Time: O(KE)|
		double[][] resPhi = new double[N][K];
		for (int i = 0; i < N; i++)
			for (int k = 0; k < K; k++) 
				resPhi[i][k] = 0;
		for (int i = 0; i < N; i++) {
			double normPhi = 0;
			for (int k = 0; k < K; k++) {
				resPhi[i][k] = Evaluation.dLogGamma(gamma[k]);		// note: phi is in log-scale 
				for (int j: trainData.getRow(i)) {
					resPhi[i][k] -= varphi[i] * trainData.get(i, j) * Math.log(beta[k][j] + Double.MIN_VALUE);
				}
				if (k != 0) normPhi = Evaluation.logSum(normPhi, resPhi[i][k]);
				else normPhi = resPhi[i][k];
			}
			for (int k = 0; k < K; k++) 					// Normalize
				resPhi[i][k] = Math.exp(resPhi[i][k] - normPhi);
		}
		
		// Update Bernoulli Parameter \varphi |Time: O(KE)|
		double[] tmpVarphi = new double[N];					// !!! TODO: most of tmpVarphi's will be zero 
		for (int i = 0; i < N; i++) {
			double di = pi[i] / (1-pi[i] + Double.MIN_VALUE) + Double.MIN_VALUE;
			double ss = Evaluation.sumSigma(i, p, q, b);
			for (int j: trainData.getRow(i)) {
				double v = Math.exp(p[i] * q[j] + b[j]) / ss;

//				Scanner sc = new Scanner(System.in);
//				System.out.println("sigma = " + v);

				for (int k = 0; k < K; k++) {
//					System.out.println("beta = " + beta[k][j] + " phi = " + phi[i][k]);
					v /= Math.pow(beta[k][j], phi[i][k]);
				}

//				System.out.println("v = " + v + "\n");
//				int gu = sc.nextInt();

				di *= Math.pow(v, trainData.get(i, j));
			}
			if (Double.isInfinite(di))
				tmpVarphi[i] = 1;
			else
				tmpVarphi[i] = 1 - 1/(1+di);
		}

		// Update All
		System.out.printf("gamma: ");
		for (int k = 0; k < K; k++) {
			gamma[k] = resGamma[k];
			System.out.printf("%f\t", gamma[k]);
		} System.out.printf("\nvarphi: ");
		for (int i = 123; i < 133; i++) {
			System.out.printf("%f\t", tmpVarphi[i]);
		} System.out.printf("\n");
		for (int i = 0; i < N; i++)
			for (int k = 0; k < K; k++) 
				phi[i][k] = resPhi[i][k];
		for (int i = 0; i < N; i++) 
			varphi[i] = tmpVarphi[i];
//			varphi[i] = 0.9;

		return;
	}

	// Update Model Parameters using Gradient Ascent 
	public static void
	modelParamEstimate(
		SparseMatrix<Integer> trainData,
		double[] alpha,	double[][] beta, double[] pi,
		double[] p, double[] q, double[] b,
		double[] gamma, double[][] phi, double[] varphi
	) {
		int K = alpha.length, N = p.length;

		// Update \beta |Time: O(KE)| 
		double[][] tmpBeta = new double[K][N];
		for (int k = 0; k < K; k++) {
			double normBeta = 0;
			for (int j = 0; j < N; j++) {
				for (int i: trainData.getColumn(j)) {
					// beta_{kj} ~~ \sum_{i} { \phi_{ij} * (1-\varphi_i) * n_{ij} } 
					tmpBeta[k][j] += trainData.get(i, j) * phi[i][k] * (1-varphi[i]);
				}
				normBeta += tmpBeta[k][j];
			}
			if (normBeta == 0) normBeta = 1;
			for (int j = 0; j < N; j++) 
				tmpBeta[k][j] /= normBeta;
		}

		// Update \pi
		for (int k = 0; k < K; k++)
			pi[k] = varphi[k];

		// Update p, q, b 
		double[] gradP = new double[N], gradQ = new double[N], gradB = new double[N];
		for (int i = 0; i < N; i++) {
			double ssw = Evaluation.sumSigmaWeighted(i, p, q, b);
			double ss = Evaluation.sumSigma(i, p, q, b);
			for (int j: trainData.getRow(i)) {
				double sg = Math.exp(p[i] * q[j] + b[j]) / ss;		// \sigma_{ij} 
				double v = trainData.get(i, j) * varphi[i];
				gradP[i] += v * (q[j] - ssw/ss);			// \sum_j { \varphi_i * n_{ij} * (q_j - weightedSum)
				gradQ[j] += v * p[i] * (1-sg);				// \sum_i { \varphi_i * n_{ij} * (1-\sigma_{ij}) * p_i 
				gradB[j] += v * (1-sg);					// \sum_i { \varphi_i * n_{ij} * (1-\sigma_{ij}) 
			}
		}

		// Update All
		for (int k = 0; k < K; k++) 
			for (int j = 0; j < N; j++)
				beta[k][j] = tmpBeta[k][j];
		for (int i = 0; i < N; i++) {
			p[i] += lr * gradP[i];
			q[i] += lr * gradQ[i];
			b[i] += lr * gradB[i];
		}
//		FileParser.output("./p", p);
//		FileParser.output("./q", q);
//		FileParser.output("./b", b);

		return;
	}

	public static void
	update (
		SparseMatrix<Integer> trainData,
		double[] alpha,	double[][] beta, double[] pi,
		double[] p, double[] q, double[] b,
		double[] gamma, double[][] phi, double[] varphi
	) {
		variationEM(trainData, alpha, beta, pi, p, q, b, gamma, phi, varphi);
		modelParamEstimate(trainData, alpha, beta, pi, p, q, b, gamma, phi, varphi);
		double l = Evaluation.calcLikelihood(trainData, alpha, beta, pi, p, q, b, gamma, phi, varphi);
		System.out.println("=== Likelihood: " + l + " ===");

		return;
	}
}
