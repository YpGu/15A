/**
	IdealPointInference.java: gives one round of parameter update in the ideal point part.
**/

import java.util.*;

public class IdealPointInference
{
	public static double LR = 0.0001;

	public static void
	update(
		SparseMatrix<Integer> trainData,
		double[] gamma,
		double[][] theta, double[][] beta,
		double[] p, double[] q, double[] b,
		double reg,									// coefficient of the regularization term 
		int iterRecord
	) {
		int N = p.length;
		double[] tmpP = new double[N], tmpQ = new double[N], tmpB = new double[N];

		if (iterRecord%10 == 0) LR *= 1.2;						// provide some chances for LR to increase 

		// Update p, q, b - O(E) 
		double[][] sg = new double[N][N];						// \sigma_{ij} 
		double[] ss = new double[N];							// \sum_{i} {\sigma_{ij}} 
		Evaluation.sumSigma(p, q, b, sg, ss);
		double[] ssw = new double[N];							// \sum_{i} {q_j \sigma_{ij}} 
		Evaluation.sumSigmaWeighted(p, q, b, ssw);

		double[] gradP = new double[N], gradQ = new double[N], gradB = new double[N];
		for (int i = 0; i < N; i++) {
			if (gamma[i] != 0) {
				for (int j: trainData.getRow(i)) {
					double v = trainData.get(i, j) * gamma[i];
					gradP[i] += v * (q[j] - ssw[i]/ss[i]);			// \sum_j { \gamma_i * n_{ij} * (q_j - weightedSum/sum)
					gradQ[j] += v * p[i] * (1-sg[i][j]);			// \sum_i { \gamma_i * n_{ij} * (1-\sigma_{ij}) * p_i 
					if (Main.USEB)
						gradB[j] += v * (1-sg[i][j]);			// \sum_i { \gamma_i * n_{ij} * (1-\sigma_{ij}) 
				}
			}
		}

		FileParser.output("./param/gradP", gradP);
		FileParser.output("./param/gradQ", gradQ);
		FileParser.output("./param/gradB", gradB);

		// Line Search 
		int count = 0;
		double lr = LR;
		while (true) {
			for (int i = 0; i < N; i++) {
				tmpP[i] = p[i] + lr * gradP[i] - reg * p[i];
				tmpQ[i] = q[i] + lr * gradQ[i] - reg * q[i];
				if (Main.USEB)
					tmpB[i] = b[i] + lr * gradB[i] - reg * b[i];
			}

			double increaseInObj = Evaluation.calcLikelihood(trainData, gamma, theta, beta, tmpP, tmpQ, tmpB, reg)
					- Evaluation.calcLikelihood(trainData, gamma, theta, beta, p, q, b, reg);
			System.out.println("\tIncrease = " + increaseInObj + ", learning rate = " + lr);
			if (increaseInObj > 0) break;

			lr *= 0.5;
			count += 1;
			if (count == 5) break;
		}
		if (count != 0) LR *= 0.5;

		// Update 
//		if (count != 5) {
		if (true) {
			for (int i = 0; i < N; i++) {
				p[i] = tmpP[i]; q[i] = tmpQ[i]; b[i] = tmpB[i];
			}
		}

		return;
	}
}
