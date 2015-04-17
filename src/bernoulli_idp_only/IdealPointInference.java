/**
	IdealPointInference.java: gives one round of parameter update in the ideal point part.
**/

import java.util.*;

public class IdealPointInference
{
	public static double LR = 0.0005;
	public static double LR_M = LR/100;

	public static double
	update(
		SparseMatrix<Integer> trainData, SparseMatrix<Integer> negData,
		double[] p, double[] q, double[] b,
		double reg,									// coefficient of the regularization term 
		int iterRecord
	) {
		int N = p.length;
		double[] tmpP = new double[N], tmpQ = new double[N], tmpB = new double[N];

//		if (iterRecord%10 == 9) LR *= 1.2;						// provide some chances for LR to increase 
		if (LR < LR_M) LR = LR_M;

		// Update p, q, b - O(E) 
		double c = N*(N-1) / (double)trainData.getSize() - 1;
//		double c = 1;
		System.out.println("c = " + c);

		double[] gradP = new double[N], gradQ = new double[N], gradB = new double[N];
		for (int i: trainData.getXDict()) {
			for (int j: trainData.getRow(i)) {
				double sigma = Evaluation.logis(p[i]*q[j]+b[j]);
				gradP[i] += q[j]*(1-sigma);
				gradQ[j] += p[i]*(1-sigma);
				gradB[j] += (1-sigma);
			}
		}
		for (int i: negData.getXDict()) {
			for (int j: negData.getRow(i)) {
				double sigma = Evaluation.logis(p[i]*q[j]+b[j]);
				gradP[i] -= q[j]*sigma*c;
				gradQ[j] -= p[i]*sigma*c;
				gradB[j] -= sigma*c;
			}
		}

		FileParser.output("./param/gradP", gradP);
		FileParser.output("./param/gradQ", gradQ);
		FileParser.output("./param/gradB", gradB);

		// Line Search 
		int count = 0;
		double lr = LR;
		double res = 0;
		while (true) {
			for (int i = 0; i < N; i++) {
				tmpP[i] = p[i] + lr * gradP[i] - reg * p[i];
				tmpQ[i] = q[i] + lr * gradQ[i] - reg * q[i];
				tmpB[i] = b[i] + lr * gradB[i] - reg * b[i];
			}
			res = Evaluation.calcLikelihood(trainData, negData, tmpP, tmpQ, tmpB);
			double increaseInObj = res - Evaluation.calcLikelihood(trainData, negData, p, q, b);
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

		return res;
	}
}
