/**
	Update.java: gives one round of parameter update.
**/

import java.util.*;

public class Update
{
	// Learning Rate
//	public static final double reg = 0.1;
	public static final double reg = 0;

	public static final int MAX_ITER_IPM = 3;
	public static final int MAX_ITER_BKG = 1;

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
		if (bkg && ipm) {								// update pi only if both BKG and IPM are used 
			double oldPi0 = pi[0], oldPi1 = pi[1];
			pi[0] = 0; pi[1] = 0;

			double[][] sg = new double[N][N];
			Evaluation.sumSigma(p, q, b, sg);

			for (int i: trainData.getXDict()) {
				while (true) {							// use multinomial distribution (REMEMBER to divide by x_ij!) 
					int powerTen = 0;
					double newPi0 = 1, newPi1 = 1;
					for (int j: trainData.getRow(i)) {
						double delta = 0;				// mixture 0
						for (int k = 0; k < K; k++) 
							delta += theta[i][k] * beta[k][j];
						newPi0 *= (delta * (1-sg[i][j]) / sg[i][j]);
						if (Double.isInfinite(newPi0)) {
							gamma[i] = 1;
							break;
						}
						if (newPi0 == 0) {
							gamma[i] = 0;
							break;
						}
						while (newPi0 >= 10) {
							newPi0 *= 0.1; powerTen += 1;
						}
						while (newPi0 < 1) {
							newPi0 *= 10; powerTen -= 1;
						}
					}
					for (int j: trainData.getXDict()) {
						newPi0 /= (1-sg[i][j]);
						if (Double.isInfinite(newPi0)) {
							gamma[i] = 1;
							break;
						}
						if (newPi0 == 0) {
							gamma[i] = 0;
							break;
						}
						while (newPi0 >= 10) {
							newPi0 *= 0.1; powerTen += 1;
						}
						while (newPi0 < 1) {
							newPi0 *= 10; powerTen -= 1;
						}
					}

					newPi0 *= oldPi0; newPi1 *= oldPi1;

					if (newPi0 + newPi1 != 0)				// update pi (next round) 	
						if (powerTen < 100)
							gamma[i] = newPi1 / (newPi0*Math.pow(10,powerTen) + newPi1);
						else
							gamma[i] = 0;
					else
						gamma[i] = 0.5;

					break;
				}
				pi[0] += (1-gamma[i]);
				pi[1] += gamma[i];
			}
			double norm = pi[0] + pi[1];
			pi[0] /= norm;
			pi[1] /= norm;
		}

		double l = 0;
		if (ipm) {
			for (int iter = 0; iter < MAX_ITER_IPM; iter++) {
				System.out.println("    *** Updating p,q,b " + iter + " ***");
				IdealPointInference.update(trainData, gamma, theta, beta, p, q, b, reg, iterRecord);
				l = Evaluation.calcLikelihood(trainData, gamma, theta, beta, p, q, b, reg);
				System.out.println("\tlogL (lower bound) = " + l);
			}
		}
		if (bkg) {
			for (int iter = 0; iter < MAX_ITER_BKG; iter++) {
				System.out.println("    *** Updating background parameters " + iter + " ***");
				BackgroundInference.update(trainData, gamma, theta, beta);
				l = Evaluation.calcLikelihood(trainData, gamma, theta, beta, p, q, b, reg);
				System.out.println("\tlogL (lower bound) = " + l);
			}
		}

		return l;
	}
}

