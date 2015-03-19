/**
	Evaluation.java: evaluate and objective function calculator.
**/

import java.util.*;
import java.lang.*;

public class Evaluation
{
	/// logistic (sigmond) function
	public static double logis(double x) {
		if (x > 100) {
			return 1;
		}
		else {
			return Math.pow(Math.E, x) / (1 + Math.pow(Math.E, x));
		}
	}


	/// calculate the derivative of log-gamma function 
	public static double dLogGamma(double x) {
		double dtmp = (x - 0.5) / (x + 4.5) + Math.log(x + 4.5) - 1;
		double ser = 1.0 + 76.18009173 / (x + 0) - 86.50532033 / (x + 1)
                       + 24.01409822 / (x + 2) - 1.231739516 / (x + 3)
                       +  0.00120858003 / (x + 4) -  0.00000536382 / (x + 5);
		double dser = -76.18009173 / (x + 0) / (x + 0)  + 86.50532033 / (x + 1) / (x + 1)
                       - 24.01409822 / (x + 2) / (x + 2) + 1.231739516 / (x + 3) / (x + 3)
                       -  0.00120858003 / (x + 4) / (x + 4) + 0.00000536382 / (x + 5) / (x + 5);
		double res = dtmp + dser / ser;
		if (res != res) System.out.println("dLog error");
		return res;
	}


	/// calculate the expectation of log(\pi) w.r.t. q
	public static double
	expt(Map<String, double[]> gamma, String p, int k) {
		double v1 = gamma.get(p)[k];
		double v2 = 0;
		for (double v: gamma.get(p)) {
			v2 += v;
		}

		return dLogGamma(v1)-dLogGamma(v2);
	}


	/// calculate the change in objective function, by changing the class of one node only 
	public static double
	changeInObj(
		SparseMatrix posData,
		double[][] eta,
		Map<String, Map<String, Double>> gamma, 
		Map<String, Integer> z,
		String x,									// node 
		int newClassLabelForX,
		double sw
	) {
//		long sTime = System.currentTimeMillis();

		if (newClassLabelForX < 0 || newClassLabelForX >= eta.length) {
			throw new ArrayIndexOutOfBoundsException();
		}

		double res = 0;
		int preX = z.get(x), curX = newClassLabelForX;

		// x -> y
		for (String y: posData.getRow(x)) {
			int curY = z.get(y);
			double gamma1 = 1-gamma.get(x).get(y);
			if (eta[preX][curY] != 0)
				res -= gamma1 * Math.log(eta[preX][curY] + Double.MIN_VALUE);
			if (eta[curX][curY] != 0)
				res += gamma1 * Math.log(eta[curX][curY] + Double.MIN_VALUE);
		}
		for (String y: posData.getRowComplement(x)) {
//		for (String y: negData.getRow(x)) {
			int curY = z.get(y);
			double gamma1 = 1-gamma.get(x).get(y);
			if (eta[preX][curY] != 0)
//				res -= sw * gamma1 * Math.log(1 - eta[preX][curY] + Double.MIN_VALUE);
				res -= gamma1 * Math.log(1 - eta[preX][curY] + Double.MIN_VALUE);
			if (eta[curX][curY] != 0)
//				res += sw * gamma1 * Math.log(1 - eta[curX][curY] + Double.MIN_VALUE);
				res += gamma1 * Math.log(1 - eta[curX][curY] + Double.MIN_VALUE);
		}

		// y -> x
		for (String y: posData.getColumn(x)) {
			int curY = z.get(y);
			double gamma1 = 1-gamma.get(y).get(x);
			if (eta[curY][preX] != 0)
				res -= gamma1 * Math.log(eta[curY][preX] + Double.MIN_VALUE);
			if (eta[curY][curX] != 0)
				res += gamma1 * Math.log(eta[curY][curX] + Double.MIN_VALUE);
		}
		for (String y: posData.getColumnComplement(x)) {
//		for (String y: negData.getColumn(x)) {
			int curY = z.get(y);
			double gamma1 = 1-gamma.get(y).get(x);
			if (eta[curY][preX] != 0)
//				res -= sw * gamma1 * Math.log(1 - eta[curY][preX] + Double.MIN_VALUE);
				res -= gamma1 * Math.log(1 - eta[curY][preX] + Double.MIN_VALUE);
			if (eta[curY][curX] != 0)
//				res += sw * gamma1 * Math.log(1 - eta[curY][curX] + Double.MIN_VALUE);
				res += gamma1 * Math.log(1 - eta[curY][curX] + Double.MIN_VALUE);
		}

//		long fTime = System.currentTimeMillis();
//		System.out.println("Time: " + (fTime-sTime));

		return res;
	}


	/// calculate the overall objective function (log-likelihood) 
	public static double 
	calcObj(
		SparseMatrix posData, SparseMatrix negData, 
		Map<String, Map<String, double[]>> phiP2Q, Map<String, Map<String, double[]>> phiQ2p, 
		double[][] matB, Map<String, double[]> gamma,
		Map<String, Double> vOut, Map<String, Double> vIn, Map<String, Double> vBias,
		Map<String, Double> pi,									// weight of ideology mixture
		double c,										// sample weight 
		double reg										// regularization coefficient 
	) {
//		long sTime = System.currentTimeMillis();
		double res = 0;
		Scanner sc = new Scanner(System.in);
		for (String x: posData.getDict()) {
			Set<String> s1 = posData.getRow(x);
			for (String y: s1) {								// x -> y
				double p1 = 0;
				for (int g = 0; g < matB.length; g++) {
					for (int h = 0; h < matB.length; h++) {
						double t1 = gamma.get(x)[g], t2 = gamma.get(y)[h];
						if (matB[g][h] != 0) {
							p1 += t1 * t2 * matB[g][h];
						}
					}
				}
				double p2 = logis(vOut.get(x) * vIn.get(y) + vBias.get(y));
				res += Math.log( (1-pi.get(x)) * p1 + pi.get(x) * p2 + Double.MIN_VALUE );
				if (res != res) {
					System.out.println("res error 1, res = " + res);
					System.out.printf("x = %s, y = %s\n", x, y);
					System.out.printf("p1 = %f, p2 = %f, pi = %f\n", p1, p2, pi.get(x));
					for (int i = 0; i < matB.length; i++) {
						for (int j = 0; j < matB.length; j++) {
							System.out.printf("%f\t", matB[i][j]);
						}
						System.out.printf("\n");
					}
					for (int i = 0; i < matB.length; i++) {
						System.out.printf("%f\t", gamma.get(x)[i]);
						System.out.printf("\n");
					}
					for (int i = 0; i < matB.length; i++) {
						System.out.printf("%f\t", gamma.get(y)[i]);
						System.out.printf("\n");
					}
					int gu = sc.nextInt();
				}
			}
			Set<String> s2 = posData.getRowComplement(x);
			for (String y: s2) {								// x !-> y
				double p1 = 0;
				for (int g = 0; g < matB.length; g++) {
					for (int h = 0; h < matB.length; h++) {
						double t1 = gamma.get(x)[g], t2 = gamma.get(y)[h];
						if (matB[g][h] != 1) {
							p1 += t1 * t2 * (1-matB[g][h]);
						}
					}
				}
				double p2 = 1 - logis(vOut.get(x) * vIn.get(y) + vBias.get(y));
				res += Math.log( (1-pi.get(x)) * p1 + pi.get(x) * p2 + Double.MIN_VALUE );
				if (res != res) {
					System.out.println("res error 2");
					System.out.printf("x = %s, y = %s\n", x, y);
					System.out.printf("p1 = %f, p2 = %f, pi = %f\n", p1, p2, pi.get(x));
					int gu = sc.nextInt();
				}
			}
		}

		// regularization
		if (reg != 0) {
			for (String x: posData.getDict()) {
				res -= 0.5 * reg * (vOut.get(x) * vOut.get(x) + vIn.get(x) * vIn.get(x));
			}
		}

		if (res != res) {
//			int NUM_BLOCKS = eta.length;
//			FileParser.output("./res/pi_" + NUM_BLOCKS + ".err", pi);
//			FileParser.output("./res/out_" + NUM_BLOCKS + ".err", vOut);
//			FileParser.output("./res/in_" + NUM_BLOCKS + ".err", vIn);
//			FileParser.output("./res/bias_" + NUM_BLOCKS + ".err", vBias);
			System.out.println("res NAN!");
		}
//		long fTime = System.currentTimeMillis();
//		System.out.println("\t\tTime = " + (fTime-sTime) + " ms");
	
		return res;
	}


	// auroc 
	public static void 
	auroc(
		SparseMatrix posData, SparseMatrix negData,
		Map<String, Double> pi,
		Map<String, Integer> z, double[][] eta,
		Map<String, Double> vOut, Map<String, Double> vIn, Map<String, Double> vBias
	) {
		Map<Integer, Double> recProbs = new HashMap<Integer, Double>();
		Set<Integer> posGroundTruth = new HashSet<Integer>();
		Set<Integer> negGroundTruth = new HashSet<Integer>();

		int tupleID = 0;
		for (String x: posData.getDict()) {
			Set<String> s1 = posData.getRow(x);
			for (String y: s1) {
				int zx = z.get(x), zy = z.get(y);
				double p1 = eta[zx][zy];
				double p2 = Evaluation.logis(vOut.get(x) * vIn.get(y) + vBias.get(y));
				double prob = (1-pi.get(x)) * p1 + pi.get(x) * p2;
				recProbs.put(tupleID, prob);
				posGroundTruth.add(tupleID);
				tupleID += 1;
			}
		}
		for (String x: negData.getDict()) {
			Set<String> s2 = negData.getRow(x);
			for (String y: s2) {
				int zx = z.get(x), zy = z.get(y);
				double p1 = eta[zx][zy];
				double p2 = Evaluation.logis(vOut.get(x) * vIn.get(y) + vBias.get(y));
				double prob = (1-pi.get(x)) * p1 + pi.get(x) * p2;
				recProbs.put(tupleID, prob);
				negGroundTruth.add(tupleID);
				tupleID += 1;
			}
		}

		double posSamples = posGroundTruth.size();
		double negSamples = negGroundTruth.size();
		System.out.println("\tSize of +'s = " + posSamples + "  Size of -'s = " + negSamples);

		// calculate AUC
		Map<Integer, Double> sortedProbs = ArrayTools.ValueComparator.sortByValue(recProbs);
		System.out.println("\tsortedProbs size = " + sortedProbs.size());
		double newX = 0, newY = 0, oldX = 0, oldY = 0;
		double upperAUC = 0, lowerAUC = 0;
		for (Map.Entry<Integer, Double> e: sortedProbs.entrySet()) {
			int curKey = e.getKey();

			if (posGroundTruth.contains(curKey)) {
				newY += 1.0/posSamples;
			}
			else if (negGroundTruth.contains(curKey)) {
				newX += 1.0/negSamples;
			}
			else {
				// check key? 
			}
			upperAUC += (newX - oldX) * newY;
			lowerAUC += (newX - oldX) * oldY;

			oldX = newX;
			oldY = newY;
		}

		System.out.println("\tAUC between " + lowerAUC + " and " + upperAUC);
		System.out.println("\tnewY = " + newY + " newX = " + newX);

		return;
	}
}
