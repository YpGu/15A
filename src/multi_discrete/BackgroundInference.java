/**
	BackgroundInference.java: gives one round of parameter update in the background part.
**/

import java.util.*;

public class BackgroundInference
{
	public static void
	update(
		SparseMatrix<Integer> trainData,
		double[][] gamma,						// N * 1
		double[][] theta,						// N * K
		double[][] beta							// K * V
	) {
		int K = theta[0].length, N = theta.length;

		double[][] tmpTheta = new double[N][K];
		double[][] tmpBeta = new double[K][N];
		double[] sumTheta = new double[N];
		double[] sumBeta = new double[K];

		for (int d: trainData.getXDict()) {
			for (int w: trainData.getRow(d)) {			// It's not necessary to update those p(z|d,w) where n(d,w) = 0 
				if (gamma[d][w] == 1) continue;
				double[] p_z_dw = new double[K];		// update p(z|d,w) 
				double p_z_dw_norm = 0;
				for (int k = 0; k < K; k++) {
					p_z_dw[k] = theta[d][k] * beta[k][w];
					p_z_dw_norm += p_z_dw[k];
				}
				if (p_z_dw_norm == 0) p_z_dw_norm = 1;
				for (int k = 0; k < K; k++) {			// update theta & beta
					p_z_dw[k] /= p_z_dw_norm;
					double v = p_z_dw[k] * trainData.get(d, w) * (1-gamma[d][w]);	// n(i,j)*p(k|i,j)*gamma(ij)(0) 
					tmpTheta[d][k] += v;			// theta: p(z|i) (un-normalized) 
					sumTheta[d] += v;
					tmpBeta[k][w] += v;			// beta: p(j|z) (un-normalized) 
					sumBeta[k] += v;
				}
			}
		}

		for (int i = 0; i < N; i++) if (sumTheta[i] == 0) sumTheta[i] = 1;
		for (int k = 0; k < K; k++) if (sumBeta[k] == 0) sumBeta[k] = 1;
		for (int i = 0; i < N; i++) 					// normalize theta
			for (int k = 0; k < K; k++) 
				theta[i][k] = tmpTheta[i][k] / sumTheta[i];
		for (int k = 0; k < K; k++) 
			for (int j = 0; j < N; j++) 
				beta[k][j] = tmpBeta[k][j] / sumBeta[k];

		return;
	}

}
