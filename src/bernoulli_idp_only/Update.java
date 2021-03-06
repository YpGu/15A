/**
	Update.java: gives one round of parameter update.
**/

import java.util.*;

public class Update
{
	// Learning Rate
//	public static final double reg = 0.001;
	public static final double reg = 0;

	public static final int MAX_ITER_IPM = 1;

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
			l = IdealPointInference.update(trainData, trainDataNeg, p, q, b, reg, iterRecord);
			System.out.println("\tlogL (lower bound) = " + l);
		}

		return l;
	}
}

