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


	/// calculate the overall objective function (log-likelihood) 
	public static double 
	calcObj(
		SparseMatrix posData, SparseMatrix negData, Map<String, double[]> theta, double[][] eta, Map<Integer, Double> rho,
		Map<String, Double> vOut, Map<String, Double> vIn, Map<String, Double> vBias,
		Map<String, Double> pi,									// weight of ideology mixture
		double c,										// sample weight 
		double reg										// regularization coefficient 
	) {
//		long sTime = System.currentTimeMillis();
		int K = eta.length;
		double res = 0;
		for (String x: posData.getDict()) {
			Set<String> s1 = posData.getRow(x);
			for (String y: s1) {								// x -> y
				double p1 = 0;
				for (int g = 0; g < K; g++) 
					for (int h = 0; h < K; h++) 
						p1 += theta.get(x)[g] * theta.get(y)[h] * eta[g][h];
				p1 *= rho.get(0);
				double p2 = logis(vOut.get(x) * vIn.get(y) + vBias.get(y));
				res += Math.log( (1-pi.get(x)) * p1 + pi.get(x) * p2 + Double.MIN_VALUE );
			}

			Set<String> s2 = posData.getRowComplement(x);
			for (String y: s2) {								// x !-> y
				double p1 = 1;
				for (int g = 0; g < K; g++) 
					for (int h = 0; h < K; h++) 
						p1 -= theta.get(x)[g] * theta.get(y)[h] * eta[g][h];
				p1 = (1-rho.get(0)) * p1 + rho.get(0);
				double p2 = 1 - logis(vOut.get(x) * vIn.get(y) + vBias.get(y));
				res += Math.log( (1-pi.get(x)) * p1 + pi.get(x) * p2 + Double.MIN_VALUE );
			}

		}
/*
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
*/
		// regularization
		if (reg != 0) {
			for (String x: posData.getDict()) {
				res -= 0.5 * reg * (vOut.get(x) * vOut.get(x) + vIn.get(x) * vIn.get(x));
			}
		}

		if (res != res) {
			int NUM_BLOCKS = eta.length;
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
		Map<String, double[]> theta, double[][] eta, Map<Integer, Double> rho,
		Map<String, Double> vOut, Map<String, Double> vIn, Map<String, Double> vBias,
		int type
	) {
		int K = eta.length;
		Map<Integer, Double> recProbs = new HashMap<Integer, Double>();
		Map<Integer, Double> recProbs1 = new HashMap<Integer, Double>();
		Map<Integer, Double> recProbs2 = new HashMap<Integer, Double>();
		Set<Integer> posGroundTruth = new HashSet<Integer>();
		Set<Integer> negGroundTruth = new HashSet<Integer>();

		int tupleID = 0;
		for (String x: posData.getDict()) {
			Set<String> s1 = posData.getRow(x);
			for (String y: s1) {
				double p1 = 0;
				for (int g = 0; g < K; g++) 
					for (int h = 0; h < K; h++) 
						p1 += theta.get(x)[g] * theta.get(y)[h] * eta[g][h];
				p1 *= rho.get(0);
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
				double p1 = 0;
				for (int g = 0; g < K; g++) 
					for (int h = 0; h < K; h++) 
						p1 += theta.get(x)[g] * theta.get(y)[h] * eta[g][h];
				p1 *= rho.get(0);
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
		System.out.printf("\tSize of +'s = %f Size of -'s = %f", posSamples, negSamples);

		// calculate AUC
		Map<Integer, Double> sortedProbs = ArrayTools.ValueComparator.sortByValue(recProbs);
		Map<Integer, Double> sortedProbs1 = ArrayTools.ValueComparator.sortByValue(recProbs1);
		Map<Integer, Double> sortedProbs2 = ArrayTools.ValueComparator.sortByValue(recProbs2);

		if (type == 1) {
			try (PrintWriter writer = new PrintWriter(new BufferedWriter(new FileWriter("./record/secMixtureTrain")))) {
				for (Map.Entry<Integer, Double> e: sortedProbs.entrySet()) {
					if (e.getKey() < sortedProbs.size()/2.0) {
						writer.printf("%s\t%f\t1\n", e.getKey(), e.getValue());
					}
					else
					{
						writer.printf("%s\t%f\t-1\n", e.getKey(), e.getValue());
					}
				}
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

		System.out.println(" sortedProbs size = " + sortedProbs.size());
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
			System.out.printf("\tUsing unified model: AUC between %f and %f", lowerAUC, upperAUC);
			System.out.println(" (newY = " + newY + " newX = " + newX + ")");
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
			System.out.printf("\tUsing block model: AUC between %f and %f", lowerAUC, upperAUC);
			System.out.println(" (newY = " + newY + " newX = " + newX + ")");
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
			System.out.printf("\tUsing ideal point model: AUC between %f and %f", lowerAUC, upperAUC);
			System.out.println(" (newY = " + newY + " newX = " + newX + ")");
		}

		return;
	}


	// au p-r c
	public static void 
	auprc(
		SparseMatrix posData, SparseMatrix negData,
		Map<String, Double> pi,
		Map<String, double[]> theta, double[][] eta, Map<Integer, Double> rho,
		Map<String, Double> vOut, Map<String, Double> vIn, Map<String, Double> vBias,
		int type
	) {
		int K = eta.length;
		Map<Integer, Double> recProbs = new HashMap<Integer, Double>();
		Map<Integer, Double> recProbs1 = new HashMap<Integer, Double>();
		Map<Integer, Double> recProbs2 = new HashMap<Integer, Double>();
		Set<Integer> posGroundTruth = new HashSet<Integer>();
		Set<Integer> negGroundTruth = new HashSet<Integer>();

		int tupleID = 0;
		for (String x: posData.getDict()) {
			Set<String> s1 = posData.getRow(x);
			for (String y: s1) {
				double p1 = 0;
				for (int g = 0; g < K; g++) 
					for (int h = 0; h < K; h++) 
						p1 += theta.get(x)[g] * theta.get(y)[h] * eta[g][h];
				p1 *= rho.get(0);
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
				double p1 = 0;
				for (int g = 0; g < K; g++) 
					for (int h = 0; h < K; h++) 
						p1 += theta.get(x)[g] * theta.get(y)[h] * eta[g][h];
				p1 *= rho.get(0);
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
		System.out.printf("\tSize of +'s = %f Size of -'s = %f", posSamples, negSamples);

		// calculate AUC
		Map<Integer, Double> sortedProbs = ArrayTools.ValueComparator.sortByValue(recProbs);
		Map<Integer, Double> sortedProbs1 = ArrayTools.ValueComparator.sortByValue(recProbs1);
		Map<Integer, Double> sortedProbs2 = ArrayTools.ValueComparator.sortByValue(recProbs2);

		if (type == 1) {
			try (PrintWriter writer = new PrintWriter(new BufferedWriter(new FileWriter("./record/prc/secMixtureTrain")))) {
				for (Map.Entry<Integer, Double> e: sortedProbs.entrySet()) {
					if (e.getKey() < sortedProbs.size()/2.0) {
						writer.printf("%s\t%f\t1\n", e.getKey(), e.getValue());
					}
					else
					{
						writer.printf("%s\t%f\t-1\n", e.getKey(), e.getValue());
					}
				}
			}
			catch (IOException e) {
				e.printStackTrace();
			}
			try (PrintWriter writer = new PrintWriter(new BufferedWriter(new FileWriter("./record/prc/secMixtureTrain_p1")))) {
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
			try (PrintWriter writer = new PrintWriter(new BufferedWriter(new FileWriter("./record/prc/secMixtureTrain_p2")))) {
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
			try (PrintWriter writer = new PrintWriter(new BufferedWriter(new FileWriter("./record/prc/secMixtureTest")))) {
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
			try (PrintWriter writer = new PrintWriter(new BufferedWriter(new FileWriter("./record/prc/secMixtureTest_p1")))) {
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
			try (PrintWriter writer = new PrintWriter(new BufferedWriter(new FileWriter("./record/prc/secMixtureTest_p2")))) {
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

		System.out.println(" sortedProbs size = " + sortedProbs.size());
		// x: recall; y: precision
		if (true) {
			double newX = 0, newY = 0, oldX = 0, oldY = 0;
			double upperAUC = 0, lowerAUC = 0;
			double numOfRetrieved = 0, numOfRelevant = 0;
			for (Map.Entry<Integer, Double> e: sortedProbs.entrySet()) {
				int curKey = e.getKey();
				numOfRetrieved += 1;

				if (posGroundTruth.contains(curKey)) {
					numOfRelevant += 1;
				}
				newY = numOfRelevant/numOfRetrieved;
				newX = numOfRelevant/(double)posSamples;

				upperAUC += (newX - oldX) * newY;
				lowerAUC += (newX - oldX) * oldY;

				oldX = newX;
				oldY = newY;
			}
			System.out.printf("\tUsing unified model: AUC between %f and %f", lowerAUC, upperAUC);
			System.out.println(" (newY = " + newY + " newX = " + newX + ")");
		}
		if (true) {
			double newX = 0, newY = 0, oldX = 0, oldY = 0;
			double upperAUC = 0, lowerAUC = 0;
			double numOfRetrieved = 0, numOfRelevant = 0;
			for (Map.Entry<Integer, Double> e: sortedProbs1.entrySet()) {
				int curKey = e.getKey();
				numOfRetrieved += 1;

				if (posGroundTruth.contains(curKey)) {
					numOfRelevant += 1;
				}
				newY = numOfRelevant/numOfRetrieved;
				newX = numOfRelevant/(double)posSamples;

				upperAUC += (newX - oldX) * newY;
				lowerAUC += (newX - oldX) * oldY;

				oldX = newX;
				oldY = newY;
			}
			System.out.printf("\tUsing block model: AUC between %f and %f", lowerAUC, upperAUC);
			System.out.println(" (newY = " + newY + " newX = " + newX + ")");
		}
		if (true) {
			double newX = 0, newY = 0, oldX = 0, oldY = 0;
			double upperAUC = 0, lowerAUC = 0;
			double numOfRetrieved = 0, numOfRelevant = 0;
			for (Map.Entry<Integer, Double> e: sortedProbs2.entrySet()) {
				int curKey = e.getKey();
				numOfRetrieved += 1;

				if (posGroundTruth.contains(curKey)) {
					numOfRelevant += 1;
				}
				newY = numOfRelevant/numOfRetrieved;
				newX = numOfRelevant/(double)posSamples;

				upperAUC += (newX - oldX) * newY;
				lowerAUC += (newX - oldX) * oldY;

				oldX = newX;
				oldY = newY;
			}
			System.out.printf("\tUsing ideal point model: AUC between %f and %f", lowerAUC, upperAUC);
			System.out.println(" (newY = " + newY + " newX = " + newX + ")");
		}

		return;
	}
}
