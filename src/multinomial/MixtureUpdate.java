/**
	MixtureUpdate.java: provide the variational E/M step for multinomial mixture training
**/

import java.util.*;

public class MixtureUpdate
{
	// Learning Rate
	public static final double lr = 0.001;

	// Logistic Function
	public static double logis(double x) {
		if (x > 100) return 1;
		else return Math.pow(Math.E, x) / (1 + Math.pow(Math.E, x));
	}

	// Calculate \sigma_{ij} 
	public static double sigma(int i, int j, double[] p, double[] q, double[] b) {
		return logis(p[i] * q[j] + b[j]);
	}

	// Magic. Do not touch. 
	// Calculate the Derivative of log-Gamma (Digamma) Function 
	public static double dLogGamma(double x) {
		if (x == 0) return Math.pow(10,-9);
		double dtmp = (x - 0.5) / (x + 4.5) + Math.log(x + 4.5) - 1;
		double ser = 1.0 + 76.18009173 / (x + 0) - 86.50532033 / (x + 1)
                       + 24.01409822 / (x + 2) - 1.231739516 / (x + 3)
                       +  0.00120858003 / (x + 4) -  0.00000536382 / (x + 5);
		double dser = -76.18009173 / (x + 0) / (x + 0)  + 86.50532033 / (x + 1) / (x + 1)
                       - 24.01409822 / (x + 2) / (x + 2) + 1.231739516 / (x + 3) / (x + 3)
                       -  0.00120858003 / (x + 4) / (x + 4) + 0.00000536382 / (x + 5) / (x + 5);
		double res = dtmp + dser / ser;
		if (res != res) {
			System.out.println("dLog error");
			System.out.println("x = " + x);
			Scanner sc = new Scanner(System.in);
			int gu = sc.nextInt(); 
		}
		return res;
	}

	// Handle Large/Small Numbers 
	public static void
	saveToPower(double[][] resPhi, int[][] resPhiPosNegPower, int i, int k) {
		if (Double.isInfinite(resPhi[i][k])) System.out.println("doge");
//		System.out.printf("%f ", resPhi[i][k]);
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
		double[][] resPhi = new double[N][K];
		int[][] resPhiPosNegPower = new int[N][K];				// pos-neg: 10 ^ (pos-neg) 
		for (int i = 0; i < N; i++)
			for (int k = 0; k < K; k++) 
				resPhi[i][k] = 1;
		double sumGamma = 0;
		for (int k = 0; k < K; k++) sumGamma += gamma[k];

		for (int i = 0; i < N; i++) {
			double eBase = Math.exp(-varphi[i]);
			for (int k = 0; k < K; k++) {
				double a = dLogGamma(gamma[k]) - dLogGamma(sumGamma);
				resPhi[i][k] *= Math.exp(a);				// exp[ \Psi(\gamma_k) - \Psi(sumOfGamma) ] 
				saveToPower(resPhi, resPhiPosNegPower, i, k);
				// TODO: a might be too large; the same with ePower 
				double ePower = 0;
				// TODO: select one j or many j`s?
				// ----- if all j`s: \sum_j {\beta_{kj}} is tooooooo small 
				for (int j: trainData.getRow(i)) {
					resPhi[i][k] *= Math.pow(1 / (beta[k][j] + Double.MIN_VALUE), trainData.getElement(i, j) * varphi[i]);
//					if (Double.isInfinite(resPhi[i][k])) {
					if (true) {
						System.out.println("beta = " + beta[k][j] + " nij = " + trainData.getElement(i,j) + " varphi = " + varphi[i]);
						System.out.println("base = " + (1/(beta[k][j]+Double.MIN_VALUE)));
						System.out.println("res = " + Math.pow(1 / (beta[k][j] + Double.MIN_VALUE), trainData.getElement(i, j) * varphi[i]));
						Scanner sc = new Scanner(System.in);
						int gu = sc.nextInt();
					}
					saveToPower(resPhi, resPhiPosNegPower, i, k);
//					ePower += trainData.getElement(i, j) * Math.log(beta[k][j] + Double.MIN_VALUE);
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

		// Update Bernoulli Parameter \varphi |Time: O(KE)|
		double[] tmpVarphi = new double[N];
		for (int i = 0; i < N; i++) {
			double di = pi[i] / (1-pi[i]);
			for (int j: trainData.getRow(i)) {
				double v = sigma(i, j, p, q, b);			// \sigma_{ij} / \prod_{k} {\beta_{kj} ^ \phi_{ik}} 
				for (int k = 0; k < K; k++)
					v /= Math.pow(beta[k][j], phi[i][k]);
				di *= Math.pow(v, trainData.getElement(i, j));
			}
			tmpVarphi[i] = 1 - 1/(1+di);
		}

		// Update All
		for (int k = 0; k < K; k++) {
			gamma[k] = resGamma[k];
			System.out.printf("%f\t", gamma[k]);
		} System.out.printf("\n");
		for (int i = 0; i < N; i++)
			for (int k = 0; k < K; k++) 
				phi[i][k] = resPhi[i][k];
		for (int i = 0; i < N; i++) 
			varphi[i] = tmpVarphi[i];

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
				for (int i: trainData.getColumn(j)) {		// TODO: check if getColumn exists
					// beta_{kj} ~~ \sum_{i} { \phi_{ij} * (1-\varphi_i) * n_{ij} } 
					tmpBeta[k][j] += trainData.getElement(i, j) * phi[i][k] * (1-varphi[i]);
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
			for (int j: trainData.getRow(i)) {
				double s = sigma(i, j, p, q, b);
				double v = trainData.getElement(i, j) * varphi[i] * (1-s);
				gradP[i] += v * q[j];					// \sum_j { \varphi_i * n_{ij} * (1-\sigma_{ij}) * q_j 
				gradQ[j] += v * p[i];					// \sum_i { \varphi_i * n_{ij} * (1-\sigma_{ij}) * p_i 
				gradB[j] += v;						// \sum_i { \varphi_i * n_{ij} * (1-\sigma_{ij}) 
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

		return;
	}
}
