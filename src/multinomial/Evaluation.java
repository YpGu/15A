/**
	Evaluation.java: Evaluation and Calculator.
**/

import java.util.*;
import java.lang.*;

public class Evaluation
{
	// Logistic Function
	public static double logis(double x) {
		double v = 1;
		if (x < 100)
			v = 1 - 1 / (1 + Math.exp(x));
		return v;
	}

	// Calculate Sum: \sum_{l} { exp(p_i * q_l + b_l) }
	public static double sumSigma(int i, double[] p, double[] q, double[] b) {
		double res = 0;
		for (int l = 0; l < p.length; l++) {
			double power = p[i] * q[l] + b[l];
			res += Math.exp(power);
		}
		return res;
	}

	// Calculate Weighted Sum: \sum_{l} { q_l * exp(p_i * q_l + b_l) }
	public static double sumSigmaWeighted(int i, double[] p, double[] q, double[] b) {
		double res = 0;
		for (int l = 0; l < p.length; l++) {
			double power = p[i] * q[l] + b[l];
			res += q[l] * Math.exp(power);
		}
		return res;
	}

	// Input: log(a) and log(b); Output: log(a+b)
	public static double logSum(double logA, double logB) {
 		if (logA < logB) return logB + Math.log(1 + Math.exp(logA-logB));
		else return logA + Math.log(1 + Math.exp(logB-logA));
	}

	// Magic. Do not touch. 
	// Calculate the Derivative of log-Gamma (Digamma) Function 
	public static double dLogGamma(double x) {
		if (x == 0) return Math.pow(10,-9);
		double dtmp = (x - 0.5) / (x + 4.5) + Math.log(x + 4.5) - 1;
		double ser = 1.0 + 76.18009173 / (x + 0) - 86.50532033 / (x + 1)
                       + 24.01409822 / (x + 2) - 1.231739516 / (x + 3)
                       +  0.00120858003 / (x + 4) -  0.00000536382 / (x + 5);
		double dser = -76.18009173 / (x + 0) / (x + 0)  + 86.50532033 / (x + 1) / (x + 1)
                       - 24.01409822 / (x + 2) / (x + 2) + 1.231739516 / (x + 3) / (x + 3)
                       -  0.00120858003 / (x + 4) / (x + 4) + 0.00000536382 / (x + 5) / (x + 5);
		double res = dtmp + dser / ser;
		if (res != res) {
			System.out.println("dLog error");
			System.out.println("x = " + x);
			Scanner sc = new Scanner(System.in);
			int gu = sc.nextInt(); 
		}
		return res;
	}

	public static double calcLikelihood(
		SparseMatrix<Integer> data,
		double[] alpha,	double[][] beta, double[] pi,
		double[] p, double[] q, double[] b,
		double[] gamma, double[][] phi, double[] varphi
	) {
		int K = alpha.length, N = p.length;
		double res = 0;
		for (int i = 0; i < N; i++) {
			double ss = Evaluation.sumSigma(i, p, q, b);
			for (int j: data.getRow(i)) {
				double p1 = 0;
				for (int k = 0; k < K; k++) 
					p1 += phi[i][k] * beta[k][j];
				double p2 = Math.exp(p[i] * q[j] + b[j]) / ss;
				// n(i,j) * log p(i,j) = n(i,j) * log{ (1-pi) * \sum_k{theta_{ik} * beta_{kj}} + pi * sigma(i,j) } 
				// \theta ~~ \phi (variational) 
				res += data.get(i, j) * Math.log( (1-pi[i]) * p1 + pi[i] * p2 + Double.MIN_VALUE );
			}
		}
		return res;
	}
}

