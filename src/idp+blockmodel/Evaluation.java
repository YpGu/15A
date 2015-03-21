/**
	Evaluation.java: evaluate and objective function calculator.
**/

import java.util.*;
import java.lang.*;
import java.io.*;

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
		SparseMatrix posData, SparseMatrix negData,
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
//		for (String y: posData.getRowComplement(x)) {
		for (String y: negData.getRow(x)) {
			int curY = z.get(y);
			double gamma1 = 1-gamma.get(x).get(y);
			if (eta[preX][curY] != 0)
				res -= sw * gamma1 * Math.log(1 - eta[preX][curY] + Double.MIN_VALUE);
//				res -= gamma1 * Math.log(1 - eta[preX][curY] + Double.MIN_VALUE);
			if (eta[curX][curY] != 0)
				res += sw * gamma1 * Math.log(1 - eta[curX][curY] + Double.MIN_VALUE);
//				res += gamma1 * Math.log(1 - eta[curX][curY] + Double.MIN_VALUE);
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
//		for (String y: posData.getColumnComplement(x)) {
		for (String y: negData.getColumn(x)) {
			int curY = z.get(y);
			double gamma1 = 1-gamma.get(y).get(x);
			if (eta[curY][preX] != 0)
				res -= sw * gamma1 * Math.log(1 - eta[curY][preX] + Double.MIN_VALUE);
//				res -= gamma1 * Math.log(1 - eta[curY][preX] + Double.MIN_VALUE);
			if (eta[curY][curX] != 0)
				res += sw * gamma1 * Math.log(1 - eta[curY][curX] + Double.MIN_VALUE);
//				res += gamma1 * Math.log(1 - eta[curY][curX] + Double.MIN_VALUE);
		}

//		long fTime = System.currentTimeMillis();
//		System.out.println("Time: " + (fTime-sTime));

		return res;
	}


	/// calculate the overall objective function (log-likelihood) 
	public static double 
	calcObj(
		SparseMatrix posData, SparseMatrix negData, double[][] eta, Map<String, Integer> z,	
		Map<String, Double> vOut, Map<String, Double> vIn, Map<String, Double> vBias,
		Map<String, Double> pi,									// weight of ideology mixture
		double c,										// sample weight 
		double reg										// regularization coefficient 
	) {
//		long sTime = System.currentTimeMillis();
		double res = 0;
		for (String x: posData.getDict()) {
			Set<String> s1 = posData.getRow(x);
			for (String y: s1) {								// x -> y
				int zx = z.get(x);
				int zy = z.get(y);
				double p1 = eta[zx][zy];
				double p2 = logis(vOut.get(x) * vIn.get(y) + vBias.get(y));
				res += Math.log( (1-pi.get(x)) * p1 + pi.get(x) * p2 + Double.MIN_VALUE );
			}
/*
			Set<String> s2 = posData.getRowComplement(x);
			for (String y: s2) {								// x !-> y
				int zx = z.get(x);
				int zy = z.get(y);
				double p1 = 1 - eta[zx][zy];
				double p2 = 1 - logis(vOut.get(x) * vIn.get(y) + vBias.get(y));
				res += Math.log( (1-pi.get(x)) * p1 + pi.get(x) * p2 + Double.MIN_VALUE );
			}
*/
		}
		for (String x: negData.getDict()) {
			Set<String> s1 = negData.getRow(x);
			for (String y: s1) {								// x -> y
				int zx = z.get(x);
				int zy = z.get(y);
				double p1 = 1 - eta[zx][zy];
				double p2 = 1 - logis(vOut.get(x) * vIn.get(y) + vBias.get(y));
				res += Math.log( (1-pi.get(x)) * p1 + pi.get(x) * p2 + Double.MIN_VALUE ) * c;
			}
		}

		// regularization
		if (reg != 0) {
			for (String x: posData.getDict()) {
				res -= 0.5 * reg * (vOut.get(x) * vOut.get(x) + vIn.get(x) * vIn.get(x));
			}
		}

		if (res != res) {
			int NUM_BLOCKS = eta.length;
			FileParser.output("./res/z_" + NUM_BLOCKS + ".err", z);
			FileParser.output("./res/eta_" + NUM_BLOCKS + ".err", eta);
			FileParser.output("./res/pi_" + NUM_BLOCKS + ".err", pi);
			FileParser.output("./res/out_" + NUM_BLOCKS + ".err", vOut);
			FileParser.output("./res/in_" + NUM_BLOCKS + ".err", vIn);
			FileParser.output("./res/bias_" + NUM_BLOCKS + ".err", vBias);
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
		Map<String, Double> vOut, Map<String, Double> vIn, Map<String, Double> vBias,
		int type
	) {
		Map<Integer, Double> recProbs = new HashMap<Integer, Double>();
		Map<Integer, Double> recProbs1 = new HashMap<Integer, Double>();
		Map<Integer, Double> recProbs2 = new HashMap<Integer, Double>();
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
				recProbs1.put(tupleID, p1);
				recProbs2.put(tupleID, p2);
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
				recProbs1.put(tupleID, p1);
				recProbs2.put(tupleID, p2);
				negGroundTruth.add(tupleID);
				tupleID += 1;
			}
		}

		double posSamples = posGroundTruth.size();
		double negSamples = negGroundTruth.size();
		System.out.println("\tSize of +'s = " + posSamples + "  Size of -'s = " + negSamples);

		// calculate AUC
		Map<Integer, Double> sortedProbs = ArrayTools.ValueComparator.sortByValue(recProbs);
		Map<Integer, Double> sortedProbs1 = ArrayTools.ValueComparator.sortByValue(recProbs1);
		Map<Integer, Double> sortedProbs2 = ArrayTools.ValueComparator.sortByValue(recProbs2);

		int aa = 0, bb = 0;
		if (type == 1) {
			try (PrintWriter writer = new PrintWriter(new BufferedWriter(new FileWriter("./record/secMixtureTrain")))) {
				for (Map.Entry<Integer, Double> e: sortedProbs.entrySet()) {
					if (e.getKey() < sortedProbs.size()/2.0) {
						writer.printf("%s\t%f\t1\n", e.getKey(), e.getValue());
						if (sortedProbs1.get(e.getKey()) > sortedProbs2.get(e.getKey())) {
							aa += 1;
						}
						else bb += 1;
					}
					else
					{
						writer.printf("%s\t%f\t-1\n", e.getKey(), e.getValue());
						if (sortedProbs1.get(e.getKey()) < sortedProbs2.get(e.getKey())) {
							aa += 1;
						}
						else bb += 1;
					}
				}
				System.out.println("aa = " + aa + " bb = " + bb);
			}
			catch (IOException e) {
				e.printStackTrace();
			}
			try (PrintWriter writer = new PrintWriter(new BufferedWriter(new FileWriter("./record/secMixtureTrain_p1")))) {
				for (Map.Entry<Integer, Double> e: sortedProbs1.entrySet()) {
					if (e.getKey() < sortedProbs1.size()/2.0) 
						writer.printf("%s\t%f\t1\n", e.getKey(), e.getValue());
					else
						writer.printf("%s\t%f\t-1\n", e.getKey(), e.getValue());
				}
			}
			catch (IOException e) {
				e.printStackTrace();
			}
			try (PrintWriter writer = new PrintWriter(new BufferedWriter(new FileWriter("./record/secMixtureTrain_p2")))) {
				for (Map.Entry<Integer, Double> e: sortedProbs2.entrySet()) {
					if (e.getKey() < sortedProbs2.size()/2.0) 
						writer.printf("%s\t%f\t1\n", e.getKey(), e.getValue());
					else
						writer.printf("%s\t%f\t-1\n", e.getKey(), e.getValue());
				}
			}
			catch (IOException e) {
				e.printStackTrace();
			}
		}
		if (type == 2) {
			try (PrintWriter writer = new PrintWriter(new BufferedWriter(new FileWriter("./record/secMixtureTest")))) {
				for (Map.Entry<Integer, Double> e: sortedProbs.entrySet()) {
					if (e.getKey() < sortedProbs.size()/2.0) 
						writer.printf("%s\t%f\t1\n", e.getKey(), e.getValue());
					else
						writer.printf("%s\t%f\t-1\n", e.getKey(), e.getValue());
				}
			}
			catch (IOException e) {
				e.printStackTrace();
			}
			try (PrintWriter writer = new PrintWriter(new BufferedWriter(new FileWriter("./record/secMixtureTest_p1")))) {
				for (Map.Entry<Integer, Double> e: sortedProbs1.entrySet()) {
					if (e.getKey() < sortedProbs1.size()/2.0) 
						writer.printf("%s\t%f\t1\n", e.getKey(), e.getValue());
					else
						writer.printf("%s\t%f\t-1\n", e.getKey(), e.getValue());
				}
			}
			catch (IOException e) {
				e.printStackTrace();
			}
			try (PrintWriter writer = new PrintWriter(new BufferedWriter(new FileWriter("./record/secMixtureTest_p2")))) {
				for (Map.Entry<Integer, Double> e: sortedProbs2.entrySet()) {
					if (e.getKey() < sortedProbs2.size()/2.0) 
						writer.printf("%s\t%f\t1\n", e.getKey(), e.getValue());
					else
						writer.printf("%s\t%f\t-1\n", e.getKey(), e.getValue());
				}
			}
			catch (IOException e) {
				e.printStackTrace();
			}
		}


		System.out.println("\tsortedProbs size = " + sortedProbs.size());
		if (true) {
			double newX = 0, newY = 0, oldX = 0, oldY = 0;
			double upperAUC = 0, lowerAUC = 0;
			for (Map.Entry<Integer, Double> e: sortedProbs.entrySet()) {
				Scanner sc = new Scanner(System.in);
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
		}
		if (true) {
			double newX = 0, newY = 0, oldX = 0, oldY = 0;
			double upperAUC = 0, lowerAUC = 0;
			for (Map.Entry<Integer, Double> e: sortedProbs1.entrySet()) {
				Scanner sc = new Scanner(System.in);
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
		}
		if (true) {
			double newX = 0, newY = 0, oldX = 0, oldY = 0;
			double upperAUC = 0, lowerAUC = 0;
			for (Map.Entry<Integer, Double> e: sortedProbs2.entrySet()) {
				Scanner sc = new Scanner(System.in);
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
		}

		return;
	}
}
