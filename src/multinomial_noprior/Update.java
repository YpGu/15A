/**
	Update.java: gives one round of parameter update.
**/

import java.util.*;

public class Update
{
	// Learning Rate
	public static final double reg = 0.01;

	public static final int MAX_ITER_IPM = 3;
	public static final int MAX_ITER_BKG = 3;

	// Estimate Variational Parameters using EM 
	public static void
	update(
		SparseMatrix<Integer> trainData,
		double[] pi,
		double[][] theta, double[][] beta,
		double[] p, double[] q, double[] b,
		boolean bkg, boolean ipm,
		int iterRecord
	) {
		int N = p.length, K = beta.length;

		// Update \pi 
		if (bkg && ipm) {							// update pi only if both BKG and IPM are used 
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
				}
				newPi0 *= (1-pi[i]);
				newPi1 *= pi[i];

				if (newPi0 + newPi1 != 0)				// update pi (next round) 	
					pi[i] = newPi1 / (newPi0 + newPi1);
				else
					pi[i] = 0.5;
			}
		}

		if (ipm) {
			for (int iter = 0; iter < MAX_ITER_IPM; iter++) {
				System.out.println("    *** Updating p,q,b " + iter + " ***");
				IdealPointInference.update(trainData, pi, theta, beta, p, q, b, reg, iterRecord);
				double l = Evaluation.calcLikelihood(trainData, pi, theta, beta, p, q, b);
				System.out.println("\tlogL (lower bound) = " + l);
			}
		}
		if (bkg) {
			for (int iter = 0; iter < MAX_ITER_BKG; iter++) {
				System.out.println("    *** Updating background parameters " + iter + " ***");
				BackgroundInference.update(trainData, pi, theta, beta);
				double l = Evaluation.calcLikelihood(trainData, pi, theta, beta, p, q, b);
				System.out.println("\tlogL (lower bound) = " + l);
			}
		}

		return;
	}
}

