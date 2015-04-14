/**
	Update.java: gives one round of parameter update.
**/

import java.util.*;

public class Update
{
	// Learning Rate
	public static final double reg = 0.001;

	public static final int MAX_ITER_IPM = 1;
	public static final int MAX_ITER_BKG = 1;

	// Estimate Variational Parameters using EM 
	public static double
	update(
		SparseMatrix<Integer> trainData,
		double[] pi, double[][] gamma,
		double[][] theta, double[][] beta,
		double[] p, double[] q, double[] b,
		boolean bkg, boolean ipm,
		int iterRecord
	) {
		int N = p.length, K = beta.length;

		// Update \pi 
		if (bkg && ipm) {							// update pi only if both BKG and IPM are used 
			double oldPi0 = pi[0], oldPi1 = pi[1];
			pi[0] = 0; pi[1] = 0;

			for (int i: trainData.getXDict()) {
				double ss = Evaluation.sumSigma(i, p, q, b);
				double newPi0 = 0, newPi1 = 0;
				for (int j: trainData.getRow(i)) {
					for (int k = 0; k < K; k++) 
						newPi0 += theta[i][k] * beta[k][j];
					double sg = 0;
					if (Main.USEB)
						sg = Math.exp(p[i] * q[j] + b[j]) / ss;	// \sigma_{ij} 
					else
						sg = Math.exp(p[i] * q[j]) / ss;	// \sigma_{ij} 
					newPi1 += sg;

					newPi0 *= oldPi0; newPi1 *= oldPi1;
					if (newPi0 + newPi1 != 0)			// update pi (next round) 	
						gamma[i][j] = newPi1 / (newPi0 + newPi1);
					else
						gamma[i][j] = 0.5;

					pi[0] += (1-gamma[i][j]);
					pi[1] += gamma[i][j];
				}
			}
			double nor = pi[0] + pi[1];
			pi[0] /= nor;
			pi[1] /= nor;
		}

		double l = 0;
		if (ipm) {
			for (int iter = 0; iter < MAX_ITER_IPM; iter++) {
				System.out.println("    *** Updating p,q,b " + iter + " ***");
				IdealPointInference.update(trainData, gamma, theta, beta, p, q, b, reg, iterRecord);
				l = Evaluation.calcLikelihood(trainData, gamma, theta, beta, p, q, b);
				System.out.println("\tlogL (lower bound) = " + l);
			}
		}
		if (bkg) {
			for (int iter = 0; iter < MAX_ITER_BKG; iter++) {
				System.out.println("    *** Updating background parameters " + iter + " ***");
				BackgroundInference.update(trainData, gamma, theta, beta);
				l = Evaluation.calcLikelihood(trainData, gamma, theta, beta, p, q, b);
				System.out.println("\tlogL (lower bound) = " + l);
			}
		}

		return l;
	}
}

