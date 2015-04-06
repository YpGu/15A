/**
	IdealPointInference.java: gives one round of parameter update in the ideal point part.
**/

import java.util.*;

public class IdealPointInference
{
	// Learning Rate
	public static final double lr = 0.001;
	public static final int MAX_ITER = 1;

	// Update Model Parameters using Gradient Ascent (after Variational Inference for all users) 
	public static void
	update(
		SparseMatrix<Integer> trainData,
		double[] pi,
		double[] p, double[] q, double[] b
	) {
		int N = p.length;

		// Update p, q, b - O(E) 
		for (int iter = 0; iter < MAX_ITER; iter++) {
			System.out.println("*** Updating pqb " + iter + " ***");
			double[] gradP = new double[N], gradQ = new double[N], gradB = new double[N];
			for (int i = 0; i < N; i++) {
				if (pi[i] != 0) {
					double ssw = Evaluation.sumSigmaWeighted(i, p, q, b);
					double ss = Evaluation.sumSigma(i, p, q, b);
					for (int j: trainData.getRow(i)) {
						double sg = Math.exp(p[i] * q[j] + b[j]) / ss;		// \sigma_{ij} 
						double v = trainData.get(i, j) * pi[i];
						gradP[i] += v * (q[j] - ssw/ss);			// \sum_j { \pi_i * n_{ij} * (q_j - weightedSum)
						gradQ[j] += v * p[i] * (1-sg);				// \sum_i { \pi_i * n_{ij} * (1-\sigma_{ij}) * p_i 
//						gradB[j] += v * (1-sg);					// \sum_i { \pi_i * n_{ij} * (1-\sigma_{ij}) 
					}
				}
			}

			// Update All
			for (int i = 0; i < N; i++) {
				p[i] += lr * gradP[i];
				q[i] += lr * gradQ[i];
//				b[i] += lr * gradB[i];
			}
//			double l = 0;
//			for (int i = 0; i < trainData.getXDictSize(); i++) {
//				l += Evaluation.calcLikelihood(trainData, i, alpha, beta, pi, p, q, b, gamma, allPhi, pi);
//			}
//			System.out.println("=== Likelihood: " + l + " ===");
		}

		return;
	}
}
