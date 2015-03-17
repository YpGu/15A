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
		long sTime = System.currentTimeMillis();

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

		long fTime = System.currentTimeMillis();
//		System.out.println("Time: " + (fTime-sTime));

		return res;
	}


	/// calculate the overall objective function 
	public static double 
	calcObj(
		SparseMatrix posData, SparseMatrix negData, double[][] eta, Map<String, Integer> z,	
		Map<String, Double> vOut, Map<String, Double> vIn, Map<String, Double> vBias,
		Map<String, Double> pi,									// weight of ideology mixture
		double c,										// sample weight 
		double reg										// regularization coefficient 
	) {
		// log likelihood
		long sTime = System.currentTimeMillis();
		double res = 0;
		for (String x: posData.getDict()) {
			Set<String> s1 = posData.getRow(x);
			for (String y: s1) {								// x -> y
				int zx = z.get(x);
				int zy = z.get(y);
				double p1 = eta[zx][zy];
				double p2 = logis(vOut.get(x) * vIn.get(y) + vBias.get(y));
				res += Math.log( (1-pi.get(x)) * p1 + pi.get(x) * p2 + Double.MIN_VALUE );

				if (pi.get(x) > 1 || pi.get(x) < 0) {
					System.out.println("pi error");
				}
				if (pi.get(y) > 1 || pi.get(y) < 0) {
					System.out.println("pi error");
				}
			}
			Set<String> s2 = posData.getRowComplement(x);
			for (String y: s2) {								// x !-> y
				int zx = z.get(x);
				int zy = z.get(y);
				double p1 = 1 - eta[zx][zy];
				double p2 = 1 - logis(vOut.get(x) * vIn.get(y) + vBias.get(y));
				res += Math.log( (1-pi.get(x)) * p1 + pi.get(x) * p2 + Double.MIN_VALUE );

				if (pi.get(x) > 1 || pi.get(x) < 0) {
					System.out.println("pi error");
				}
				if (pi.get(y) > 1 || pi.get(y) < 0) {
					System.out.println("pi error");
				}
			}
		}
/*		for (String x: negData.getDict()) {
			Set<String> s2 = posData.getRow(x);
			for (String y: s2) {								// x !-> y
				int zx = z.get(x);
				int zy = z.get(y);
				double p1 = eta[zx][zy];
				double p2 = logis(vOut.get(x) * vIn.get(y) + vBias.get(y));
				res += Math.log( 1 - (1-pi.get(x)) * p1 - pi.get(y) * p2 + Double.MIN_VALUE ) * c;
			}
		}
*/

		// regularization
		if (reg != 0) {
			for (String x: posData.getDict()) {
				res -= 0.5 * reg * (vOut.get(x) * vOut.get(x) + vIn.get(x) * vIn.get(x));
			}
		}

		if (res != res) {
/*
			output("./" + rel + "Res/pi2", pi_2);
			output("./" + rel + "Res/alpha", alpha);
			output("./" + rel + "Res/beta", beta);
			output("./" + rel + "Res/vOut", vOut);
			output("./" + rel + "Res/vIn", vIn);
			output("./" + rel + "Res/vBias", vBias);
*/
			System.out.println("res NAN!");
			Scanner myInput = new Scanner(System.in);
			int s = myInput.nextInt();
		}
		long fTime = System.currentTimeMillis();
		System.out.println("time = " + (fTime-sTime));
		
		return res;
	}


	// auroc  (TODO) 
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
		System.out.println("Size of +'s = " + posSamples);
		double negSamples = negGroundTruth.size();
		System.out.println("Size of -'s = " + negSamples);

		// calculate AUC
		Map<Integer, Double> sortedProbs = ArrayTools.ValueComparator.sortByValue(recProbs);
		System.out.println("sortedProbs size = " + sortedProbs.size());
		double newX = 0, newY = 0, oldX = 0, oldY = 0;
		double upperAUC = 0, lowerAUC = 0;
		for (Map.Entry<Integer, Double> e: sortedProbs.entrySet()) {
			int curKey = e.getKey();

//			System.out.println(e.getValue());
//			Scanner sc = new Scanner(System.in);
//			int gu = sc.nextInt();

			if (posGroundTruth.contains(curKey)) {
				newY += 1.0/posSamples;
			}
			else {
				newX += 1.0/negSamples;
			}
			upperAUC += (newX - oldX) * newY;
			lowerAUC += (newX - oldX) * oldY;

			oldX = newX;
			oldY = newY;
		}

		System.out.println("AUC between " + lowerAUC + " and " + upperAUC);
		System.out.println("newY = " + newY + " newX = " + newX);

		return;
	}
}
