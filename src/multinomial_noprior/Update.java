/**
	Update.java: gives one round of parameter update.
**/

import java.util.*;

public class Update
{
	// Learning Rate
	public static final double reg = 0.001;

	public static final int MAX_ITER_IPM = 3;
	public static final int MAX_ITER_BKG = 3;

	// Estimate Variational Parameters using EM 
	public static double
	update(
		SparseMatrix<Integer> trainData,
		double[] pi, double[] gamma,
		double[][] theta, double[][] beta,
		double[] p, double[] q, double[] b,
		boolean bkg, boolean ipm,
		int iterRecord
	) {
		int N = p.length, K = beta.length;

		// Update \pi 
		if (bkg && ipm) {							// update pi only if both BKG and IPM are used 
			pi[0] = 0; pi[1] = 0;

			for (int i: trainData.getXDict()) {
				if (false) {
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
					newPi0 *= (1-pi[i]); newPi1 *= pi[i];

					if (newPi0 + newPi1 != 0)					// update pi (next round) 	
						pi[i] = newPi1 / (newPi0 + newPi1);
					else
						pi[i] = 0.5;
				}

				if (true) {							// use multinomial distribution (REMEMBER to divide by x_ij!) 
					double ss = Evaluation.sumSigma(i, p, q, b);
					double newPi0 = 1, newPi1 = 1;
					for (int j: trainData.getRow(i)) {
						double delta = 0;				// mixture 0
						for (int k = 0; k < K; k++) 
							delta += theta[i][k] * beta[k][j];
						newPi0 *= delta;
						double sg = 0;					// mixture 1
						if (Main.USEB)
							sg = Math.exp(p[i] * q[j] + b[j]) / ss;	// \sigma_{ij} 
						else
							sg = Math.exp(p[i] * q[j]) / ss;	// \sigma_{ij} 
						newPi1 *= sg;
					}
					newPi0 *= pi[0]; newPi1 *= pi[1];

					if (newPi0 + newPi1 != 0)				// update pi (next round) 	
						gamma[i] = newPi1 / (newPi0 + newPi1);
					else
						gamma[i] = 0.5;
				}
				pi[0] += (1-gamma[i]);
				pi[1] += gamma[i];
			}
		}

		double l = 0;
		if (ipm) {
			for (int iter = 0; iter < MAX_ITER_IPM; iter++) {
				System.out.println("    *** Updating p,q,b " + iter + " ***");
				IdealPointInference.update(trainData, gamma, theta, beta, p, q, b, reg, iterRecord);
				l = Evaluation.calcLikelihood(trainData, pi, theta, beta, p, q, b);
				System.out.println("\tlogL (lower bound) = " + l);
			}
		}
		if (bkg) {
			for (int iter = 0; iter < MAX_ITER_BKG; iter++) {
				System.out.println("    *** Updating background parameters " + iter + " ***");
				BackgroundInference.update(trainData, gamma, theta, beta);
				l = Evaluation.calcLikelihood(trainData, pi, theta, beta, p, q, b);
				System.out.println("\tlogL (lower bound) = " + l);
			}
		}

		return l;
	}
}

