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
		SparseMatrix<Integer> trainData, SparseMatrix<Integer> trainDataNeg,
		double[] p, double[] q, double[] b,
		int iterRecord
	) {
		double l = 0;
		for (int iter = 0; iter < MAX_ITER_IPM; iter++) {
			System.out.println("    *** Updating p,q,b " + iter + " ***");
			IdealPointInference.update(trainData, trainDataNeg, p, q, b, reg, iterRecord);
			l = Evaluation.calcLikelihood(trainData, trainDataNeg, p, q, b);
			System.out.println("\tlogL (lower bound) = " + l);
		}

		return l;
	}
}

