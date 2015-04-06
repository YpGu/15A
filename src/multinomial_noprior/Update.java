/**
	Update.java: gives one round of parameter update.
**/

import java.util.*;

public class Update
{
	// Estimate Variational Parameters using EM 
	public static void
	update(
		SparseMatrix<Integer> trainData,
		double[] pi,
		double[][] theta, double[][] beta,
		double[] p, double[] q, double[] b
	) {
		int N = p.length, K = beta.length;

		// Update \pi 
		for (int i: trainData.getXDict()) {
			double ss = Evaluation.sumSigma(i, p, q, b);
			double newPi0 = 0, newPi1 = 0;

			for (int j: trainData.getRow(i)) {
				for (int k = 0; k < K; k++) {
					newPi0 += theta[i][k] * beta[k][j];
				}

				double sg = Math.exp(p[i] * q[j] + b[j]) / ss;		// \sigma_{ij} 
				newPi1 += sg;
			}
			newPi0 *= (1-pi[i]);
			newPi1 *= pi[i];

			pi[i] = newPi1 / (newPi0 + newPi1);				// update pi (next round) 
		}

		IdealPointInference.update(trainData, pi, p, q, b);
		BackgroundInference.update(trainData, pi, theta, beta);

		return;
	}
}

